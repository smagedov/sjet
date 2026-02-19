#ifndef SJET_OBJECTPT_
#define SJET_OBJECTPT_

namespace sjet {	
    template<class T>
    double objectPt(const T& particle)
    {
        return particle.pt();
    }
}

#endif // SJET_OBJECTPT_
