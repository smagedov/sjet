#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <cmath>

#include "MomentPreprocessor.hh"
#include "InvariantMomentCalculator.hh"

using namespace std;


//------------------------------------------------------------
// Shift jet in detector coordinates
//------------------------------------------------------------

std::vector<RawMomentPoint>
transformJet(
    const std::vector<RawMomentPoint>& jet,
    double deta,
    double dphi)
{
    auto out = jet;

    for(auto& p : out)
    {
        p.eta += deta;

        p.phi += dphi;

        // Wrap phi into [-pi,pi]
        while(p.phi > M_PI)
            p.phi -= 2.0*M_PI;

        while(p.phi < -M_PI)
            p.phi += 2.0*M_PI;
    }

    return out;
}


//------------------------------------------------------------

Invariants calculate(
    const std::vector<RawMomentPoint>& jet)
{
    auto normalized =
        MomentPreprocessor::preprocess(jet);

    Moments moments =
        InvariantMomentCalculator::compute(
            normalized,
            5,
            5);

    return
        InvariantMomentCalculator::computeInvariants(
            moments);
}


//------------------------------------------------------------

double relativeError(
    double a,
    double b)
{
    double denom =
        std::max(std::fabs(a),1e-12);

    return std::fabs(a-b)/denom;
}


//------------------------------------------------------------

int main()
{

    //--------------------------------------------------------
    // Jet deliberately placed near phi boundary
    //--------------------------------------------------------

    vector<RawMomentPoint> jet =
    {
        { 0.22,  3.05, 12.0 },
        {-0.13,  3.10,  8.0 },
        { 0.44, -3.12,  5.0 },
        {-0.37, -3.05, 10.0 },
        { 0.06,  3.13,  6.0 },
        {-0.42, -3.10,  4.0 },
        { 0.28,  3.00,  3.0 }
    };


    //--------------------------------------------------------
    // Reference
    //--------------------------------------------------------

    Invariants reference =
        calculate(jet);


    //--------------------------------------------------------
    // Random detector transformations
    //--------------------------------------------------------

    mt19937 rng(12345);


    uniform_real_distribution<double>
        etaShift(-5.0,5.0);


    uniform_real_distribution<double>
        phiShift(-M_PI,M_PI);


    double maxError = 0.0;


    //--------------------------------------------------------
    // Monte Carlo
    //--------------------------------------------------------

    for(int trial=0; trial<1000; ++trial)
    {

        double deta =
            etaShift(rng);

        double dphi =
            phiShift(rng);


        auto transformed =
            transformJet(
                jet,
                deta,
                dphi);


        Invariants inv =
            calculate(transformed);


        for(size_t i=0;i<reference.value.size();++i)
        {

            double err =
                relativeError(
                    reference.value[i],
                    inv.value[i]);


            maxError =
                std::max(
                    maxError,
                    err);


            if(err > 1e-8)
            {
                cout
                    << "FAIL\n\n"
                    << "Trial      : "
                    << trial << "\n"
                    << "Invariant  : "
                    << i+1 << "\n"
                    << setprecision(16)
                    << "Reference  : "
                    << reference.value[i] << "\n"
                    << "Computed   : "
                    << inv.value[i] << "\n"
                    << "Rel Error  : "
                    << err
                    << endl;

                return 1;
            }
        }
    }


    //--------------------------------------------------------

    cout
        << "PASS\n\n";

    cout
        << "All 18 invariants remained invariant under\n"
        << "1000 random eta translations and phi rotations.\n\n";


    cout
        << scientific
        << "Maximum relative error = "
        << maxError
        << endl;


    return 0;
}
