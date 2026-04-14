#ifndef GAUSSIANJETV1_HH_
#define GAUSSIANJETV1_HH_

#include <cmath>

#include "AbsParticleCollectionMaker.hh"
#include "ParticleMaker.hh"

// Simple distribution in which pt-weighted delta eta and delta phi
// will be jointly Gaussian distributed with the same predefined width
template <class Particle, class Rng=std::mt19937_64>
class GaussianJetV1 : public PoissonParticleCollectionMaker<Particle,Rng>
{
public:
    inline GaussianJetV1(const bool randomizeMultiplicity = true,
                         const double m = MEAN_PION_MASS)
        : eta0_(0.0), phi0_(0.0), width_(0.0), averagePt_(0.0),
          randomizeMult_(randomizeMultiplicity)
    {
        setParticleMass(m);
    }

    // "averagePt" is the average expected Pt value among all
    // particles in the jet
    inline GaussianJetV1(const double eta, const double phi,
                         const double width, const double averagePt,
                         const double particleMultiplicity,
                         const bool randomizeMultiplicity = true,
                         const double m = MEAN_PION_MASS)
        : eta0_(eta), phi0_(phi), randomizeMult_(randomizeMultiplicity)
    {
        setWidth(width);
        setMagnitude(averagePt);
        Base::setMultiplicity(particleMultiplicity);
        setParticleMass(m);
    }

    inline virtual ~GaussianJetV1() override {}

    inline virtual GaussianJetV1* clone() const override
        {return new GaussianJetV1(*this);}

    inline virtual bool multiplicityIsRandom() const override
        {return randomizeMult_;}

    inline virtual void setLocation(const double eta,
                                    const double phi) override
    {
        eta0_ = eta;
        phi0_ = phi;
    }

    inline virtual void setWidth(const double widthParameter) override
    {
        assert(widthParameter > 0.0);
        width_ = widthParameter;
    }

    inline virtual void setParticleMass(const double m)
    {
        assert(m >= 0.0);
        m_ = m;
    }

    inline virtual void setMagnitude(const double magnitudeParameter) override
    {
        assert(magnitudeParameter > 0.0);
        averagePt_ = magnitudeParameter;
    }

protected:
    virtual Particle makeOneParticle(Rng& gen) const override
    {
        assert(width_ > 0.0);
        assert(averagePt_ > 0.0);
        assert(m_ >= 0.0);

        // Generate delta eta and delta phi from the jet axis
        // assuming that they are Gaussian distributed
        const double s = width_*M_SQRT2;
        std::normal_distribution d(0.0, s);
        const double deta = d(gen);
        const double dphi = d(gen);

        // Average particle Pt for this shift from the jet axis,
        // also suppressed according to the same Gaussian. The overall
        // effect is that the width of jet Pt flow becomes "width_".
        const double avePt = averagePt_*2.0*exp(-(deta*deta + dphi*dphi)/2.0/s/s);

        // Generate particle Pt from the exponential distribution
        std::exponential_distribution<> expd(1.0/avePt);
        const double pt = expd(gen);

        return ParticleMaker<Particle>::make(pt, eta0_+deta, phi0_+dphi, m_);
    }

private:
    typedef PoissonParticleCollectionMaker<Particle,Rng> Base;

    double eta0_;
    double phi0_;
    double width_;
    double averagePt_;
    double m_;
    bool randomizeMult_;
};

#endif // GAUSSIANJETV1_HH_
