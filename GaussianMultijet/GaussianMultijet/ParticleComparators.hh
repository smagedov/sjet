#ifndef PARTICLECOMPARATORS_HH_
#define PARTICLECOMPARATORS_HH_

template<class Particle>
struct LessByPt
{
    inline bool operator()(const Particle& l, const Particle& r) const
        {return l.pt() < r.pt();}  
};

template<class Particle>
struct LessByM
{
    inline bool operator()(const Particle& l, const Particle& r) const
        {return l.m() < r.m();}  
};

#endif // PARTICLECOMPARATORS_HH_
