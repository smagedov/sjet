#ifndef SJET_STABILITYFUNC_HH_
#define SJET_STABILITYFUNC_HH_

#include <cmath>

namespace sjet {
	template <typename stability>
	double stability(stability dist1, stability dist2) {
		return std::log(dist1/dist1)
	};
}

#endif //SJET_STABILITYFUNC_HH_
