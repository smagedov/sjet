#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "InvariantMomentCalculator.hh"

using namespace std;

std::vector<MomentPoint> rotate(
    const std::vector<MomentPoint>& points,
    double theta)
{
    std::vector<MomentPoint> out = points;

    double c = cos(theta);
    double s = sin(theta);

    for(auto &p : out)
    {
        double x = p.x;
        double y = p.y;

        p.x = c*x - s*y;
        p.y = s*x + c*y;
    }

    return out;
}

void printMoments(const Moments& M)
{
    cout << fixed << setprecision(10);

    cout << "M20 = " << M.value[2][0] << endl;
    cout << "M02 = " << M.value[0][2] << endl;
    cout << "M11 = " << M.value[1][1] << endl;

    cout << endl;
}

void printInvariants(const Invariants& I)
{
    cout << fixed << setprecision(12);

    for(size_t i=0;i<I.value.size();i++)
    {
        cout
            << "I"
            << i+1
            << " = "
            << I.value[i]
            << endl;
    }
}

int main()
{
    //-------------------------------------------------------
    // A completely arbitrary point cloud
    //-------------------------------------------------------

    vector<MomentPoint> cloud =
    {
        {-0.40,-0.20,2.0},
        { 0.35,-0.10,5.0},
        {-0.10, 0.45,4.0},
        { 0.55, 0.20,3.0},
        {-0.25, 0.30,7.0},
        { 0.15,-0.45,6.0}
    };

    //-------------------------------------------------------
    // Original
    //-------------------------------------------------------

    Moments M0 =
        InvariantMomentCalculator::compute(
            cloud,
            5,
            5);

    Invariants I0 =
        InvariantMomentCalculator::computeInvariants(
            M0);

    //-------------------------------------------------------
    // Rotate 90 degrees
    //-------------------------------------------------------

    auto rotated =
        rotate(cloud,M_PI/2.0);

    Moments M1 =
        InvariantMomentCalculator::compute(
            rotated,
            5,
            5);

    Invariants I1 =
        InvariantMomentCalculator::computeInvariants(
            M1);

    //-------------------------------------------------------
    // Print low-order moments
    //-------------------------------------------------------

    cout << "Original moments\n";
    printMoments(M0);

    cout << "Rotated moments\n";
    printMoments(M1);

    //-------------------------------------------------------
    // Compare invariants
    //-------------------------------------------------------

    cout << "\nInvariant comparison\n\n";

    double maxError = 0.0;

    for(size_t i=0;i<I0.value.size();i++)
    {
        double err =
            fabs(I0.value[i]-I1.value[i]);

        maxError =
            max(maxError,err);

        cout
            << setw(2)
            << i+1
            << "  "
            << setw(18)
            << I0.value[i]
            << "  "
            << setw(18)
            << I1.value[i]
            << "  "
            << err
            << endl;
    }

    cout << "\nMaximum absolute error = "
         << maxError
         << endl;

    return 0;
}
