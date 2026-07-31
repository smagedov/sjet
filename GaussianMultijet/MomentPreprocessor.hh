#ifndef MOMENTPREPROCESSOR_HH_
#define MOMENTPREPROCESSOR_HH_

#include <vector>
#include <cmath>

#include "InvariantMomentCalculator.hh"

struct RawMomentPoint
{
    double eta;
    double phi;
    double weight;
};

class MomentPreprocessor
{
public:

    static std::vector<MomentPoint>
    preprocess(const std::vector<RawMomentPoint>& particles)
    {
        std::vector<MomentPoint> output;

        if (particles.empty())
            return output;

        //----------------------------------------------------------
        // Compute weighted centroid
        //----------------------------------------------------------

        double totalWeight = 0.0;

        double etaCenter = 0.0;
        double phiX = 0.0;
        double phiY = 0.0;
	//double phiCenter = 0.0;

        for (const auto& p : particles)
        {
            totalWeight += p.weight;

            etaCenter += p.weight * p.eta;

            phiX += p.weight * std::cos(p.phi);
            phiY += p.weight * std::sin(p.phi);
	    //phiCenter += p.weight * p.phi;
        }

        if (totalWeight <= 0.0)
            return output;

        etaCenter /= totalWeight;
	//phiCenter /= totalWeight;

        double phiCenter = std::atan2(phiY, phiX);

        //----------------------------------------------------------
        // Compute covariance matrix
        //----------------------------------------------------------

        double mu20 = 0.0;
        double mu02 = 0.0;
        double mu11 = 0.0;

        for (const auto& p : particles)
        {
            double dx = p.eta - etaCenter;
            double dy = deltaPhi(p.phi, phiCenter);
	    //double dy = p.phi - phiCenter;

            mu20 += p.weight * dx * dx;
            mu02 += p.weight * dy * dy;
            mu11 += p.weight * dx * dy;
        }

        mu20 /= totalWeight;
        mu02 /= totalWeight;
        mu11 /= totalWeight;

        //----------------------------------------------------------
        // Principal axis rotation
        //----------------------------------------------------------

        double theta =
            0.5 *
            std::atan2(
                2.0 * mu11,
                mu20 - mu02);

        double ct = std::cos(theta);
        double st = std::sin(theta);

        //----------------------------------------------------------
        // Overall scale
        //----------------------------------------------------------

        double sigma =
            std::sqrt(mu20 + mu02);

        if (sigma < 1e-8)
            sigma = 1.0;

        //----------------------------------------------------------
        // Transform every point
        //----------------------------------------------------------

        output.reserve(particles.size());

        for (const auto& p : particles)
        {
            double x = p.eta - etaCenter;
            double y = deltaPhi(p.phi, phiCenter);
	    //double y = p.phi - phiCenter;

            // Rotate

            double xr =  ct * x + st * y;
            double yr = -st * x + ct * y;
	    //double xr = x;
	    //double yr = y;

            // Normalize

            //xr /= sigma;
            //yr /= sigma;

            output.push_back(
            {
                xr,
                yr,
                p.weight
            });
        }

        return output;
    }

private:

    static double deltaPhi(double phi1,
                           double phi2)
    {
        double dphi = phi1 - phi2;

        while (dphi > M_PI)
            dphi -= 2.0 * M_PI;

        while (dphi < -M_PI)
            dphi += 2.0 * M_PI;

        return dphi;
    }
};

#endif // MOMENTPREPROCESSOR_HH_
