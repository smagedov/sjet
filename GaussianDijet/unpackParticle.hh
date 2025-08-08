#ifndef UNPACKPARTICLE_HH_
#define UNPACKPARTICLE_HH_

template <class Particle>
inline void unpackParticle(const Particle& p,
                           double& pt, double& eta,
                           double& phi, double& m)
{
    pt = p.pt();
    eta = p.eta();
    phi = p.phi();
    m = p.m();
}

#endif // UNPACKPARTICLE_HH_
