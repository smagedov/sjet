///=========================================================================
/// diffusionFactor2.h
///
/// The factor for delta_ij for calculating d_ij in the diffusion
/// recombination algorithm (fast version by Chebyshev polynomials)
///
/// I. Volobouev
/// May 2012
///=========================================================================

#ifndef DIFFUSIONFACTOR2_H_
#define DIFFUSIONFACTOR2_H_

#ifdef __cplusplus
extern "C" {
#endif

double diffusionFactor2(double minPtToMaxPtRatio);

#ifdef __cplusplus
}
#endif

#endif // DIFFUSIONFACTOR2_H_
