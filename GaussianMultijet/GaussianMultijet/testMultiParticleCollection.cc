#include <iostream>

#include "GaussianJetV1.hh"
#include "ExponentialPileup.hh"
#include "MultiParticleCollectionMaker.hh"

#include "CmdLine.hh"

typedef rk::P4 My4Vector;
typedef GaussianJetV1<My4Vector> MyJet;
typedef ExponentialPileup<My4Vector> MyPileup;

static void print_usage(const char* progname)
{
    std::cout << "\nUsage: " << progname << " mult1 mult2 mult3\n"
              << std::endl;
}

static void print_particle(const My4Vector& p, std::ostream& os)
{
    os << p.pt() << ' ' << p.eta() << ' ' << p.phi();
}

static void make_and_print(const MultiParticleCollectionMaker<My4Vector>& ev)
{
    std::random_device rd{};
    std::mt19937_64 gen{rd()};

    const std::pair<std::vector<My4Vector>, std::vector<unsigned> >& particles =
        ev.make(gen);

    std::cout.precision(10);
    for (const My4Vector& p : particles.first)
    {
        print_particle(p, std::cout);
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

    double mult1, mult2, mult3;

    try {
        cmdline.optend();

        if (cmdline.argc() != 3)
            throw CmdLineError("wrong number of command line arguments");
        cmdline >> mult1 >> mult2 >> mult3;

        if (mult1 < 0.0 || mult2 < 0.0 || mult2 < 0.0)
            throw CmdLineError("multiplicity parameters must be non-negative");
    }
    catch (const CmdLineError& e) {
        std::cerr << "Error in " << cmdline.progname() << ": "
                  << e.str() << std::endl;
        return 1;
    }
    
    MultiParticleCollectionMaker<My4Vector> eventBuilder;

    const MyJet jet1(-2.0, 0.0, 0.5, 20.0, mult1);
    eventBuilder.add(jet1);

    const MyJet jet2(2.0, 0.0, 1.0, 10.0, mult2);
    eventBuilder.add(jet2);

    const MyPileup pileup(5.0, 2.0, mult3);
    eventBuilder.add(pileup);

    make_and_print(eventBuilder);

    return 0;
}
