#ifndef PYTHIAJET_HH_
#define PYTHIAJET_HH_

#include <cmath>

#include "AbsParticleCollectionMaker.hh"
#include "ParticleMaker.hh"

// Simple distribution in which pt-weighted delta eta and delta phi
// will be jointly Gaussian distributed with the same predefined width
template <class Particle, class Rng=std::mt19937_64>
class PythiaJet : public PoissonParticleCollectionMaker<Particle,Rng>
{
public:
    inline PythiaJet(const bool randomizeMultiplicity = true,
		    const double m = MEAN_PION_MASS)
        : eta0_(0.0), phi0_(0.0), width_(0.0), averagePt_(0.0),
	  randomizeMult_(randomizeMultiplicity)
    {
    }

    inline virtual ~PythiaJet() override {}

    inline virtual PythiaJet* clone() const override
        {return new PythiaJet(*this);}

protected:
    virtual Particle makeOneParticle(Rng& gen) const override
    {
        return ParticleMaker<Particle>::make(pt, eta, phi, m);
    }

private:
    typedef PoissonParticleCollectionMaker<Particle,Rng> Base;
}

#endif // PYTHIAJET_HH_
