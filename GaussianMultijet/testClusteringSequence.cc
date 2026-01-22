#include <cmath>
#include <iostream>
#include <cstdlib>

#include "rk/rk.hh"

#include "sjet/ClusteringSequence.hh"
#include "sjet/DistanceCalculator.hh"
#include "sjet/NodeVisitors.hh"

#include "CmdLine.hh"

typedef sjet::ClusteringSequence<sjet::DistanceCalculator, rk::P4> MyClustSeq;

static void print_usage(const char* progname)
{
    std::cout << "\nUsage: " << progname << " nTimes nParticles maxDist\n"
              << std::endl;
}

rk::P4 random_particle()
{
    const double pt = drand48()*9.0 + 1.0;
    const double eta = drand48()*6.0 - 3.0;
    const double phi = drand48()*2.0*M_PI;
    const double m = 1.0 + drand48()*4.0;
    return rk::P4(pt*geom3::Vector3(cos(phi), sin(phi), sinh(eta)), m);
}

static void run_test(const unsigned nTries,
                     const unsigned nParticles,
                     const double maxDist)
{
    for (unsigned itry=0; itry<nTries; ++itry)
    {
        std::vector<rk::P4> particles;
        particles.reserve(nParticles);
        for (unsigned i=0; i<nParticles; ++i)
            particles.push_back(random_particle());

        MyClustSeq seq((sjet::DistanceCalculator()));
        seq.init(particles);
        if (maxDist <= 0.0)
            seq.run();
        else
            seq.run(maxDist);

        MyClustSeq seq2((sjet::DistanceCalculator()));
        seq2.init(particles);
        if (maxDist <= 0.0)
        {
            seq2.run(1.0);
            seq2.run();
        }
        else
        {
            seq2.run(maxDist/2.0);
            seq2.run(maxDist);
        }
        assert(seq.clustHist() == seq2.clustHist());

        if (maxDist > 0.0)
        {
            const std::vector<unsigned>& clusDist0 = seq2.clusterIndices();
            for (unsigned iclus : clusDist0)
            {
                sjet::ParticleAdder<rk::P4> adder;
                seq2.visitParentParticles(adder, iclus);
                const double m2 = (adder.result() - seq2.clustHist()[iclus].p()).squared();
                assert(fabs(m2) < 1.0e-8);
            }
            const std::vector<unsigned>& clusAssign0 = seq2.clusterAssignments(maxDist);

            seq2.run(maxDist*1.5);
            const std::vector<unsigned>& clusDist1 = seq2.clusterIndices(maxDist);
            assert(clusDist0 == clusDist1);
            const std::vector<unsigned>& clusAssign1 = seq2.clusterAssignments(maxDist);
            assert(clusAssign1 == clusAssign0);
            
            seq2.run();
            const std::vector<unsigned>& clusDist2 = seq2.clusterIndices(maxDist);
            assert(clusDist0 == clusDist2);
            const std::vector<unsigned>& clusAssign2 = seq2.clusterAssignments(maxDist);
            assert(clusAssign2 == clusAssign0);
        }

        std::cout << "\nEvent " << itry << std::endl;
        std::cout << "nParticles = " << seq.nParticles() << std::endl;
        std::cout << "nClusters = " << seq.nClusters() << std::endl;
        std::cout << "Last distance is " << seq.lastDistance() << std::endl;
        std::cout << "Next distance is " << seq.nextDistance() << std::endl;
        std::cout << "Recombination distances are:";
        const std::vector<double>& recoDist = seq.recombinationDistances();
        for (const double& d : recoDist)
            std::cout << ' ' << d;
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[])
{
    CmdLine cmdline(argc, argv);
    if (argc == 1)
    {
        print_usage(cmdline.progname());
        return 0;
    }

    unsigned nTries, nParticles;
    double maxDist;

    try {
        cmdline.optend();

        // We are expecting 3 arguments on the command line
        if (cmdline.argc() != 3)
            throw CmdLineError("wrong number of command line arguments");
        cmdline >> nTries >> nParticles >> maxDist;
    }
    catch (const CmdLineError& e) {
        std::cerr << "Error in " << cmdline.progname() << ": "
                  << e.str() << std::endl;
        return 1;
    }

    run_test(nTries, nParticles, maxDist);
    
    return 0;
}
