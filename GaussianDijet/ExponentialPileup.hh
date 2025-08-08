#ifndef EXPONENTIALPILEUP_HH_
#define EXPONENTIALPILEUP_HH_

#include <cmath>

#include "AbsParticleCollectionMaker.hh"
#include "ParticleMaker.hh"

template <class Particle, class Rng=std::mt19937_64>
class ExponentialPileup : public PoissonParticleCollectionMaker<Particle,Rng>
{
public:
    inline ExponentialPileup(const bool randomizeMultiplicity = true,
                             const double m = MEAN_PION_MASS)
        : etaWidth_(0.0), averagePt_(0.0),
          randomizeMult_(randomizeMultiplicity)
    {
        setParticleMass(m);
    }

    inline ExponentialPileup(const double etaWidth, const double averagePt,
                             const double particleMultiplicity,
                             const bool randomizeMultiplicity = true,
                             const double m = MEAN_PION_MASS)
        : randomizeMult_(randomizeMultiplicity)
    {
        setWidth(etaWidth);
        setMagnitude(averagePt);
        Base::setMultiplicity(particleMultiplicity);
        setParticleMass(m);
    }

    inline virtual ~ExponentialPileup() override {}

    inline virtual ExponentialPileup* clone() const override
        {return new ExponentialPileup(*this);}

    inline virtual bool multiplicityIsRandom() const override
        {return randomizeMult_;}

    inline virtual void setLocation(double /* eta */,
                                    double /* phi */) override
    {
        // This call is ignored by this class
    }

    inline virtual void setWidth(const double widthParameter) override
    {
        assert(widthParameter >= 0.0);
        etaWidth_ = widthParameter;
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
        assert(etaWidth_ >= 0.0);
        assert(averagePt_ > 0.0);
        assert(m_ >= 0.0);

        std::exponential_distribution<> expd(1.0/averagePt_);
        const double pt = expd(gen);

        std::uniform_real_distribution<> uni(-1.0, 1.0);
        const double eta = etaWidth_*uni(gen);
        const double phi = M_PI*uni(gen);
        
        return ParticleMaker<Particle>::make(pt, eta, phi, m_);
    }

private:
    typedef PoissonParticleCollectionMaker<Particle,Rng> Base;

    double etaWidth_;
    double averagePt_;
    double m_;
    bool randomizeMult_;
};

#endif // EXPONENTIALPILEUP_HH_
