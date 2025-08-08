#ifndef PRINTPTETAPHI_HH_
#define PRINTPTETAPHI_HH_

#include <iostream>

template <class Particle>
inline void printPtEtaPhi(const Particle& p, std::ostream& os)
{
    os << p.pt() << ' ' << p.eta() << ' ' << p.phi();
}

template <class Particle>
inline void printPtEtaPhiM(const Particle& p, std::ostream& os)
{
    printPtEtaPhi(p, os);
    os << ' ' << p.m();
}

#endif // PRINTPTETAPHI_HH_
