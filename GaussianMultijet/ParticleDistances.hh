#ifndef PARTICLEDISTANCES_HH_
#define PARTICLEDISTANCES_HH_

#include <cmath>

inline double deltaPhi(const double phi1, const double phi2)
{
    double delta = phi1 - phi2;
    while (delta > M_PI)
        delta -= 2.0*M_PI;
    while (delta < -M_PI)
        delta += 2.0*M_PI;
    return delta;
}

struct ParticleDeltaR
{
    template <class Particle>
    inline double operator()(const Particle& left, const Particle& right) const
    {
        const double deta = left.eta() - right.eta();
        const double dphi = deltaPhi(left.phi(), right.phi());
        return std::sqrt(dphi*dphi + deta*deta);
    }
};

struct MassSquaredDist
{
    template <class Particle>
    inline double operator()(const Particle& left, const Particle& right) const
    {
        const double msq = (left - right).squared();
        return std::abs(msq);
    }
};

#endif // PARTICLEDISTANCES_HH_
