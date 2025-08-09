#ifndef SJET_DISTANCECALCULATOR_HH_
#define SJET_DISTANCECALCULATOR_HH_

#include <cmath>
#include <algorithm>

#include "sjet/diffusionFactor2.h"

namespace sjet {	
    struct DistanceCalculator {
        template<class Particle>
        inline double operator()(const Particle& p1, const Particle& p2) const
        {
            static constexpr double twopi = 2.0*M_PI;

            const double deta = p2.eta() - p1.eta();
            double dphi = p2.phi() - p1.phi();
            while (dphi > M_PI)
                dphi -= twopi;
            while (dphi < -M_PI)
                dphi += twopi;
            const double dR = sqrt(deta*deta + dphi*dphi);
            const double pt1 = p1.pt();
            const double pt2 = p2.pt();
            return dR*diffusionFactor2(std::min(pt1, pt2)/std::max(pt1, pt2));
        }
    };
}

#endif // SJET_DISTANCECALCULATOR_HH_
