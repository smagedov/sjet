#ifndef ABSPARTICLECOLLECTIONMAKER_HH_
#define ABSPARTICLECOLLECTIONMAKER_HH_

#include <cassert>
#include <vector>
#include <random>

// Abstract class for generating particle collections such as
// jets, unfderlying event, noise hits, etc
template <class Particle, class Rng=std::mt19937_64>
struct AbsParticleCollectionMaker
{
    typedef Particle particle_type;
    typedef Rng rng_type;

    inline AbsParticleCollectionMaker() {}
    inline virtual ~AbsParticleCollectionMaker() {}

    // "Virtual copy constructor"
    virtual AbsParticleCollectionMaker* clone() const = 0;

    // Accessors
    virtual double getMultiplicity() const = 0;
    virtual bool multiplicityIsRandom() const = 0;

    // Modifiers
    virtual void setLocation(double eta, double phi) = 0;
    virtual void setWidth(double widthParameter) = 0;
    virtual void setMagnitude(double magnitudeParameter) = 0;
    virtual void setMultiplicity(const double multiplicityParameter) = 0;

    // Main generation function
    virtual std::vector<Particle> make(Rng& eng) const = 0;
};

// Abstract class for generating particle collections which
// takes care of randomizing particle multiplicity
template <class Particle, class Rng=std::mt19937_64>
class PoissonParticleCollectionMaker : public AbsParticleCollectionMaker<Particle,Rng>
{
public:
    inline PoissonParticleCollectionMaker() : mult_(0.0) {}

    inline PoissonParticleCollectionMaker(const double multiplicityParameter)
        : mult_(multiplicityParameter) {assert(mult_ >= 0.0);}

    // "Virtual copy constructor"
    virtual PoissonParticleCollectionMaker* clone() const = 0;

    inline virtual ~PoissonParticleCollectionMaker() {}

    inline virtual double getMultiplicity() const
        {return mult_;}

    // You can disable multiplicity randomization by overriding
    // this method and returning "false"
    inline virtual bool multiplicityIsRandom() const
        {return true;}

    virtual void setLocation(double eta, double phi) = 0;
    virtual void setWidth(double widthParameter) = 0;
    virtual void setMagnitude(double magnitudeParameter) = 0;

    inline virtual void setMultiplicity(const double multiplicityParameter)
    {
        assert(multiplicityParameter >= 0.0);
        mult_ = multiplicityParameter;
    }

    // This implementation of "make" assumes that particles are independent
    // of each other and that their number is Poisson distributed.
    // Of course, conservation laws do not work in this simple scenario.
    virtual std::vector<Particle> make(Rng& eng) const
    {
        const double mult = getMultiplicity();
        assert(mult >= 0.0);
        std::vector<Particle> result;
        if (mult > 0.0)
        {
            unsigned n;
            if (multiplicityIsRandom())
            {
                std::poisson_distribution<unsigned> poiss(mult);
                n = poiss(eng);
            }
            else
                n = static_cast<unsigned>(mult);
            result.reserve(n);
            for (unsigned i=0; i<n; ++i)
                result.push_back(makeOneParticle(eng));
        }
        return result;
    }

protected:
    virtual Particle makeOneParticle(Rng& gen) const = 0;

private:
    double mult_;
};

#endif // ABSPARTICLECOLLECTIONMAKER_HH_
