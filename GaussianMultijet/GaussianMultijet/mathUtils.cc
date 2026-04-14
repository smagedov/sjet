#include <cmath>

#include "mathUtils.hh"

double safepow(const double value, const int power)
{
    if (power < 0 && !value)
        return 0.0;
    else
        return pow(value, power);
}
