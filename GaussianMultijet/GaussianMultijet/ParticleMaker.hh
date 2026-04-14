#ifndef PARTICLEMAKER_HH_
#define PARTICLEMAKER_HH_

#include <cassert>

#include "rk/rk.hh"

#define MEAN_PION_MASS 0.138

template <class Particle>
struct ParticleMaker
{
    typedef Particle particle_type;

    static Particle make(double pt, double eta,
                         double phi, double m);
};

template <>
struct ParticleMaker<rk::P4>
{
    typedef rk::P4 particle_type;

    inline static rk::P4 make(const double pt, const double eta,
                              const double phi, const double m)
    {
        assert(pt >= 0.0);
        assert(m >= 0.0);
        return rk::P4(geom3::Vector3(pt*cos(phi), pt*sin(phi), pt*sinh(eta)), m);
    }
};

#endif // PARTICLEMAKER_HH_
