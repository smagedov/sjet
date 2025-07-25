#ifndef SM_DISTANCECALCULATOR_HH_
#define SM_DISTANCECALCULATOR_HH_

#include <algorithm>
#include "rk/rk.hh"
#include "diffusionFactor2.h"

namespace sjet {
	
	struct DistanceCalculator {
		inline double operator()(const rk::P4& p1, const rk::P4& p2) const {
			const double deltaR = geom3::deltaR(p1.momentum(), p2.momentum());
			const double pT1 = p1.pt();
			const double pT2 = p2.pt();
			return deltaR*diffusionFactor2(std::min(pT1,pT2)/std::max(pT1,pT2));
		}
	};
}

#endif // SM_DISTANCECALCULATOR_HH_
