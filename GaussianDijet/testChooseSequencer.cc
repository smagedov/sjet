#include <iostream>

#include "CmdLine.hh"
#include "ChooseSequencer.hh"
#include "permutation.hh"

static void print_usage(const char* progname)
{
    std::cout << "\nUsage: " << progname << " [-c] [-p] m n\n"
              << std::endl
              << "  -c      this switch tells the program to count and print\n" 
              << "          the total number of choices/permutations\n" 
              << std::endl
              << "  -p      this switch tells the program to permute elements\n" 
              << "          in the set assignments\n" 
              << std::endl;
}

int main(int argc, char *argv[])
{
    CmdLine cmdline(argc, argv);
    if (argc == 1)
    {
        print_usage(cmdline.progname());
        return 0;
    }

    unsigned m, n;
    bool printCount, permute;

    try {
        printCount = cmdline.has("-c");
        permute = cmdline.has("-p");

        cmdline.optend();
        cmdline >> m >> n;

        if (permute && m > 20U)
            throw CmdLineError("Argument m is too large");
        if (m > n)
            throw CmdLineError("Must have m <= n");
    }
    catch (const CmdLineError& e) {
        std::cerr << "Error in " << cmdline.progname() << ": "
                  << e.str() << std::endl;
        return 1;
    }

    std::vector<unsigned> chosen(m);
    std::vector<unsigned> p(m);
    std::vector<unsigned> permuted(m);
    const unsigned long nPerms = factorial(m);

    unsigned long count = 0;
    for (ChooseSequencer seq(m, n); seq.isValid(); ++seq)
    {
        seq.getChoice(m ? &chosen[0] : nullptr, chosen.size());

        if (m && permute)
            for (unsigned long iperm=0; iperm<nPerms; ++iperm)
            {
                orderedPermutation(iperm, &p[0], p.size());
                for (unsigned i=0; i<m; ++i)
                    permuted[i] = chosen[p[i]];
                for (unsigned i=0; i<m; ++i)
                {
                    if (i) std::cout << ' ';
                    std::cout << permuted[i];
                }
                std::cout << '\n';
                ++count;
            }                
        else
        {
            for (unsigned i=0; i<m; ++i)
            {
                if (i) std::cout << ' ';
                std::cout << chosen[i];
            }
            std::cout << '\n';
            ++count;
        }
    }

    if (printCount)
        std::cout << count
                  << (permute ? " permutations" : " choices")
                  << " total\n";

    return 0;
}
