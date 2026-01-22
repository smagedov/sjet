#ifndef MULTIPARTICLECOLLECTIONMAKER_HH_
#define MULTIPARTICLECOLLECTIONMAKER_HH_

#include <memory>
#include <utility>

#include "AbsParticleCollectionMaker.hh"

template <class Particle, class Rng=std::mt19937_64>
class MultiParticleCollectionMaker
{
public:
    inline MultiParticleCollectionMaker() {}

    inline unsigned nMakers() const {return makers_.size();}

    inline void add(const AbsParticleCollectionMaker<Particle,Rng>& i_maker)
        {makers_.emplace_back(i_maker.clone());}

    inline AbsParticleCollectionMaker<Particle,Rng>& maker(
        const unsigned i)
        {return *makers_.at(i);}

    inline const AbsParticleCollectionMaker<Particle,Rng>& maker(
        const unsigned i) const
        {return *makers_.at(i);}

    // The first element of the returned pair is the combined sequence
    // of particles. The second element specifies how many particles
    // were produced by each maker.
    std::pair<std::vector<Particle>, std::vector<unsigned> >
    make(Rng& eng) const
    {
        std::pair<std::vector<Particle>, std::vector<unsigned> > result;
        const unsigned nMk = makers_.size();
        result.second.reserve(nMk);
        for (unsigned i=0; i<nMk; ++i)
        {
            const std::vector<Particle>& made = makers_[i]->make(eng);
            result.first.insert(result.first.end(), made.begin(), made.end());
            result.second.push_back(made.size());
        }
        return result;
    }

private:
    std::vector<std::shared_ptr<AbsParticleCollectionMaker<Particle,Rng> > > makers_;
};

#endif // MULTIPARTICLECOLLECTIONMAKER_HH_
