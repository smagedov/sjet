// Include necessary system headers
#include <iostream>
#include <limits>
#include <random>

// Relativistic kinematics header (for 4-vectors, etc)
#include "rk/rk.hh"

// Event structure for this analysis
#include "PythiaEvent.hh"

// Analysis sequencer
#include "FrameworkAnalysisSequence.hh"

// Event production and analysis modules
#include "PythiaEventMaker.hh"
#include "PythiaEventTTbarGen.hh"
#include "RunDiffusionClustering.hh"
#include "EventCounter.hh"
#include "GenParticleDump.hh"
#include "PythiaOutputPrinter.hh"
#include "PythiaJetMaker.hh"
#include "RandJaccardMMDCalculator.hh"
#include "ClusteringVisualisation.hh"
#include "FullJetHistoryPrinter.hh"

// Utilities, etc
#include "stringUtils.hh"
#include "tv_to_usec.hh"
#include "gnuc_demangle_type.hh"

// Command line procesing
#include "CmdLine.hh"

// Function for printing usage instructions
static void print_usage(const char* progname)
{
    std::cout << "\nUsage: " << progname
              << " [-m] [-r]"
              << " [--countPeriod value]"
              << " [--dumpFile value]"
              << " [--maxDistance value]"
              << " [--pileupEtaRange value]"
              << " [--ptMin value]"
              << " [--refJetWidth value]"
              << " [--seed value]"
              << " avePtRef muRef avePtProbe muProbe avePtPileup muPileup"
              << " probeJetWidth deltaR nEvents nJets distanceCutoff outputFile\n"
              << std::endl
              << "The following options can be specified to modify the program behavior:\n"
              << std::endl
              << " -m                This switch turns off Poisson randomization of the jet\n"
              << "                   and pileup multiplicities.\n"
              << std::endl
              << " -r                This switch turns off randomization of the eta and phi\n"
              << "                   location of the probe jet.\n"
              << std::endl
              << " --dumpFile        This option can be used to specify the name of a file into\n"
              << "                   which generated particles will be printed in a text form.\n"
              << "                   Default value of this option is an empty string which means\n"
              << "                   that no such file will be created.\n"
              << std::endl
              << " --countPeriod     If the value of this option is positive (0 is the default),\n"
              << "                   a message about processed event count will be printed to\n"
              << "                   the standard output every time \"value\" more events are\n"
              << "                   processed.\n"
              << std::endl
              << " --maxDistance     This option allows you to specify the maximum clustering\n"
              << "                   distance. If this option is not specified, clustering will\n"
              << "                   run until only one top-level cluster remains.\n"
              << std::endl
              << " --pileupEtaRange  This option allows you to specify the pseudorapidity range\n"
              << "                   for pileup particles. They will be generated inside the eta\n"
              << "                   interval [-value, value]. Default value of this option\n"
              << "                   is 5.0.\n"
              << std::endl
              << " --ptMin           This option allows you to specify the minimum pt cutoff\n"
              << "                   for the event particles (mimicking detector inefficiency).\n"
              << "                   Default value of this option is 0.0 (cluster all paricles).\n"
              << std::endl
              << " --refJetWidth     This option allows you to specify the width of the\n"
              << "                   reference jet in the eta-phi space. Default value of this\n"
              << "                   option is 0.2.\n"
              << std::endl
              << " --seed            This option allows you to provide the seed for the random\n"
              << "                   number generator. Default value of this option is 0 which\n"
              << "                   means that the seed will be chosen automatically (resulting\n"
              << "                   in a different set of events every time you run the\n"
              << "                   program). If you want to specify the seed, choose a large\n"
              << "                   prime number (but less than 2^64).\n"
              << std::endl
              << "The following arguments must be provided on the command line:\n"
              << std::endl
              << " avePtRef          Average pt (per particle) for the reference jet particles\n"
              << "                   (before minimum pt cutoff is applied).\n"
              << std::endl
              << " muRef             Average particle multiplicity for the reference jet (before\n"
              << "                   minimum pt cutoff is applied). The actual multiplicity will\n"
              << "                   be Poisson distributed with this average.\n"
              << std::endl
              << " avePtProbe        Average pt for the probe jet particles (before pt cutoff).\n"
              << std::endl
              << " muProbe           Average particle multiplicity for the probe jet (Poisson).\n"
              << std::endl
              << " avePtPileup       Average pt for the pileup particles. Their pt values will\n"
              << "                   be exponentially distributed with this average (before pt\n"
              << "                   cutoff).\n"
              << std::endl
              << " muPileup          Average particle multiplicity for the pileup (Poisson).\n"
              << std::endl
              << " probeJetWidth     Width of the probe jet in the eta-phi space.\n"
              << std::endl
              << " deltaR            The expected distance between reference and probe jet\n"
              << "                   directions in the eta-phi space (the actual distance\n"
              << "                   will fluctuate about this value).\n"
              << std::endl
              << " nEvents           Number of events to generate and analyze.\n"
              << std::endl
              << " distanceCutoff    Distance cutoff for the simple diffusion clustering.\n"
              << std::endl
              << " outputFile        The name of the file into which event information will be\n"
              << "                   written in a text form.\n"
              << std::endl
	      << " nJets             The number of jets to generate.\n"
	      << std::endl;
}

int main(int argc, char *argv[])
{
    typedef std::mt19937_64 MyRng;
    // typedef std::mt19937 MyRng;
    typedef PythiaEvent<rk::P4> MyEvent;
    typedef frw::FrameworkAnalysisSequence<MyEvent> MyAnalysis;

    // The convention is to print usage instructions and
    // exit with status 0 in case the program is run without
    // any command line arguments
    CmdLine cmdline(argc, argv);
    if (argc == 1)
    {
        print_usage(cmdline.progname());
        return 0;
    }

    // Parse command line arguments
    //
    // The seed for the random number generator
    unsigned long seed = 0;

    // How often we are going to print event counts?
    unsigned countPeriod = 0;

    // Number of events to generate
    unsigned nEvents;

    // Number of jets to generate
    unsigned nJets;

    // Maximum distance for the clustering process
    double maxDistance = std::numeric_limits<double>::max();

    // Distance for the "naive" jet definition based
    // only on the distance cutoff
    double distanceCutoff;

    std::string inputFile;

    // File names for dumping information
    std::string dumpFile, outputFile;

    try {
        // Process short options

        // Process long options
        cmdline.option(0, "--dumpFile") >> dumpFile;
        cmdline.option(0, "--countPeriod") >> countPeriod;
        cmdline.option(0, "--maxDistance") >> maxDistance;
        cmdline.option(0, "--seed") >> seed;

        // Process positional arguments
        cmdline.optend();
        if (cmdline.argc() != 5)
            throw CmdLineError("wrong number of command line arguments");
        cmdline >> nEvents >> nJets >> distanceCutoff >> outputFile >> inputFile;

        if (distanceCutoff <= 0.0 || maxDistance <= 0.0)
            throw CmdLineError("jet recombination distance cutoffs must be positive");

        if (distanceCutoff > maxDistance)
            throw CmdLineError("jet definition distance cutoff is too large");
    }
    catch (const CmdLineError& e) {
        std::cerr << "Error in " << cmdline.progname() << ": "
                  << e.str() << std::endl;
        return 1;
    }

    // Print the command line arguments
    std::cout << "Command line: " << cmdline.progname();
    for (int iarg=1; iarg<argc; ++iarg)
        std::cout << ' ' << argv[iarg];
    std::cout << std::endl;
    
    // Storage space for the event data
    MyEvent evt;
    std::cout << "\nEvent type is \"" << gnuc_demangle_type(evt) << '"' << std::endl;
    std::cout << "Event version is " << evt.version() << std::endl;

    // Initialize the random number generator
    unsigned long actualSeed = seed;
    if (!seed)
    {
        std::random_device rd;
        actualSeed = rd();
    }
    std::cout << "\nRng seed is " << actualSeed
              << (seed ? "" : " (auto)") << "\n\n";
    MyRng gen(actualSeed);

    // Create the event generation and analysis sequence.
    // The "add" function of the FrameworkAnalysisSequence class
    // clones the framework modules. After that the original
    // modules should immediately go out of scope. This is because
    // modules could have unique and shared pointers as their members,
    // while copy constructors could be auto-generated by the compiler.
    // In such cases calling the destructors at the end of the
    // program (rather than immediately) could potentially have
    // undesired effects.
    MyAnalysis mySequence;
    {
        const PythiaEventMaker<MyEvent, MyRng> PythiaJetGenerator(
            "PythiaEventMaker", gen, inputFile);
        mySequence.add(PythiaJetGenerator);
    }
    mySequence.add(PythiaEventTTbarGen<MyEvent>(
                            "PythiaEventTTbarGen"));
    if (countPeriod)
        mySequence.add(frw::EventCounter<MyEvent>("Counter_0", countPeriod));
    mySequence.add(RunDiffusionClustering<MyEvent>(
                       "DiffusionClustering", maxDistance));
    mySequence.add(PythiaJetMaker<MyEvent>(
                       "PythiaJetMaker", distanceCutoff));
    mySequence.add(PythiaOutputPrinter<MyEvent>(
    		       "PythiaOutputPrinter", "output.txt"));
    if (!dumpFile.empty())
        mySequence.add(GenParticleDump<MyEvent>(
                           "GenParticleDump", dumpFile));
    mySequence.add(FullJetHistoryPrinter<MyEvent>("0", "pythiatHistFile", 1000, distanceCutoff));
    mySequence.add(ClusteringVisualisation<MyEvent>(
                       "ClusteringVisualisation", outputFile));
    //mySequence.add(frw::EventCounter<MyEvent>("Last_counter"));

    // Print module labels in this analysis sequence
    std::cout << "Sequence of framework modules ("
              << mySequence.nModules() << " total):" << std::endl;
    for (const std::string& l : mySequence.labels())
        std::cout << l << std::endl;

    // Note the start time of the event loop
    struct timeval tv0;
    const int st0 = gettimeofday(&tv0, NULL);
    assert(st0 == 0);

    // Generate events and run the analysis code
    // (i.e., process the event loop)
    for (unsigned iev=0; iev<nEvents; ++iev)
    {
        evt.clear();
        evt.setNumber(iev);
        mySequence.run(evt);
    }


    // Measure the run time of the event loop
    struct timeval tv1;
    const int st1 = gettimeofday(&tv1, NULL);
    assert(st1 == 0);
    const double sec = (tv_to_usec(&tv1) - tv_to_usec(&tv0))/1.0e6;
    std::cout << "\nEvent cycle run time is " << sec << " sec / "
              << nEvents << " events" << std::endl;

    // Run the "endJob" methods of the framework modules
    std::cout << "\nReports from analysis modules:" << std::endl;
    mySequence.endJob();

    // We are done
    return 0;
}


