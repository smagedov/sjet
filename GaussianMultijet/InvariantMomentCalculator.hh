#ifndef INVARIANTMOMENTCALCULATOR_HH_
#define INVARIANTMOMENTCALCULATOR_HH_

#include <vector>
#include <cmath>

struct MomentPoint
{
    double x;       // Normalized x coordinate
    double y;       // Normalized y coordinate
    double weight;  // Particle weight (typically pT)
};

struct Moments
{
    std::vector<std::vector<double>> value;
};

struct Invariants
{
    std::vector<double> value;
};


class InvariantMomentCalculator
{
public:

    static Moments compute(const std::vector<MomentPoint>& points,
                           int Nmax,
                           int Mmax)
    {
        Moments M;
        M.value.assign(Nmax + 1,
                       std::vector<double>(Mmax + 1, 0.0));

        for (const auto& p : points)
        {
            std::vector<double> psiX(Nmax + 1);
            std::vector<double> psiY(Mmax + 1);

            for (int n = 0; n <= Nmax; ++n)
                psiX[n] = hermitePsi(n, p.x);

            for (int m = 0; m <= Mmax; ++m)
                psiY[m] = hermitePsi(m, p.y);

            for (int n = 0; n <= Nmax; ++n)
            {
                for (int m = 0; m <= Mmax; ++m)
                {
                    M.value[n][m] +=
                        p.weight * psiX[n] * psiY[m];
                }
            }
        }

        return M;
    }

    static Invariants computeInvariants(const Moments& M)
    {
	    Invariants I;
	    I.value.reserve(18);
	    
	    const double M02 = M.value[0][2];
	    const double M03 = M.value[0][3];
	    const double M04 = M.value[0][4];
	    const double M05 = M.value[0][5];
	    
	    const double M11 = M.value[1][1];
	    const double M12 = M.value[1][2];
	    const double M13 = M.value[1][3];
	    const double M14 = M.value[1][4];
	    
	    const double M20 = M.value[2][0];
	    const double M21 = M.value[2][1];
	    const double M22 = M.value[2][2];
	    const double M23 = M.value[2][3];
	    
	    const double M30 = M.value[3][0];
	    const double M31 = M.value[3][1];
	    const double M32 = M.value[3][2];
	    
	    const double M40 = M.value[4][0];
	    const double M41 = M.value[4][1];
	    
	    const double M50 = M.value[5][0];
	    
	    const double A = M30 + M12;
	    const double B = M21 + M03;
	    const double C = M31 + M13;
	    const double D = M31 - M13;
	    const double E = M40 - M04;
	    const double F = M40 - 6.0*M22 + M04;
	    const double G = M50 + 2.0*M32 + M14;
	    const double H = M41 + 2.0*M23 + M05;
	    const double J = M50 - 2.0*M32 - 3.0*M14;
	    const double K = 3.0*M41 + 2.0*M23 - M05;
	    const double L = M50 - 10.0*M32 + 5.0*M14;
	    const double R = 5.0*M41 - 10.0*M23 + M05;

	    //Second/Third-Order Invariants
	    I.value.push_back(M20 + M02);
	    I.value.push_back(A*A + B*B);
	    I.value.push_back((M20-M02)*(A*A-B*B)+4.0*M11*A*B);
	    I.value.push_back(M11*(A*A-B*B)-(M20-M02)*A*B);
	    I.value.push_back((M30-3*M12)*A*(A*A-3*B*B)+(M03-3*M21)*B*(B*B-3*A*A));
	    I.value.push_back((M30-3*M12)*B*(B*B-3*A*A)+(3*M21-M03)*A*(A*A-3*B*B));

	    //Fourth-Order Invariants
	    I.value.push_back(M40 + 2*M22 + M04);
	    I.value.push_back(E*(A*A-B*B)+4*C*A*B);
	    I.value.push_back(C*(A*A-B*B)-E*A*B);
	    I.value.push_back(F*(pow(A,4)-6*A*A*B*B+pow(B,4))+16*D*A*B*(A*A-B*B));
	    I.value.push_back(F*A*B*(B*B-A*A)+D*(pow(A,4)-6*A*A*B*B+pow(B,4)));

	    //Fifth-Order Invariants
	    I.value.push_back(G*G + H*H);
	    I.value.push_back(G*A + H*B);
	    I.value.push_back(H*A - G*B);
	    I.value.push_back(J*(A*A*A-3*A*B*B)-K*(B*B*B-3*A*A*B));
	    I.value.push_back(J*(B*B*B-3*A*A*B)+K*(A*A*A-3*A*B*B));
	    I.value.push_back(L*(pow(A,5)-10*pow(A,3)*B*B+5*A*pow(B,4))+R*(pow(B,5)-10*A*A*pow(B,3)+5*pow(A,4)*B));
	    I.value.push_back((M05-10*M23+5*M41)*(pow(A,5)-10*pow(A,3)*B*B+5*A*pow(B,4))-(5*M14-10*M32+M50)*(pow(B,5)-10*A*A*pow(B,3)+5*pow(A,4)*B));

	    return I;
    }

private:

    static double hermiteH(int n, double x)
    {
        if (n == 0)
            return 1.0;

        if (n == 1)
            return 2.0 * x;

        double Hnm1 = 1.0;
        double Hn   = 2.0 * x;

        for (int k = 1; k < n; ++k)
        {
            double Hnp1 =
                2.0 * x * Hn
                - 2.0 * k * Hnm1;

            Hnm1 = Hn;
            Hn   = Hnp1;
        }

        return Hn;
    }

    static double factorial(int n)
    {
        double f = 1.0;

        for (int i = 2; i <= n; ++i)
            f *= i;

        return f;
    }

    static double hermitePsi(int n, double x)
    {
        const double H = hermiteH(n, x);

        const double norm =
            std::sqrt(
                std::pow(2.0, n)
                * factorial(n)
                * std::sqrt(M_PI));

        return
            H
            * std::exp(-0.5 * x * x)
            / norm;
    }
};

#endif //INVARIENTMOMENTCALCULATOR_HH_
