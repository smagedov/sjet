#ifndef INVARIENTMOMENTS_HH_
#define INVARIENTMOMENTS_HH_

#include <vector>
#include <cmath>
#include <iostream>

#include "cnpy.h"
#include "ParticleDistances.hh"

struct Moments {
    std::vector<std::vector<double>> real;
    std::vector<std::vector<double>> imag;
};

template <class Event>
class InvarientMoments : public frw::AbsFrameworkAnalyzer<Event>
{
private:
        typedef typename Event::clust_seq_type::cluster_type cluster_type;
public:
        typedef Event event_type;
        typedef frw::AbsFrameworkAnalyzer<Event> Base;

        inline InvarientMoments(const std::string& i_label)
                : Base(i_label)
        {
        }

        inline virtual ~InvarientMoments() override {}

        inline virtual InvarientMoments* clone() const override
                {return new InvarientMoments(*this);}

	void getLeaves(int node,
                        const std::vector<cluster_type>& history,
                        std::vector<int>& leaves)
        {
                if (node < 0 || node >= static_cast<int>(history.size())) {
                        return;
                }

                const auto& c = history[node];

                int p1 = c.parent1();
                int p2 = c.parent2();

                if (p1 < 0 && p2 < 0) {
                        leaves.push_back(node);
                        return;
                }

                if (p1 >= 0) {
                        getLeaves(p1, history, leaves);
                }

                if (p2 >= 0) {
                        getLeaves(p2, history, leaves);
                }
        }

	double hermiteH(int n, double x) {
		if (n == 0) return 1.0;
		if (n == 1) return 2.0 * x;
		
		double Hnm1 = 1.0;
		double Hn = 2.0 * x;
		
		for (int k = 1; k < n; ++k) {
			double Hnp1 = 2.0 * x * Hn - 2.0 * k * Hnm1;
			Hnm1 = Hn;
			Hn = Hnp1;
		}
		
		return Hn;
	}

	double factorial(int n) {
		double f = 1.0;
		for (int i = 2; i <= n; ++i) f *= i;
		return f;
	}

	double hermitePsi(int n, double x) {
		double H = hermiteH(n, x);
		double norm = std::sqrt(std::pow(2.0, n) * factorial(n) * std::sqrt(M_PI));
		return (H * std::exp(-0.5 * x * x)) / norm;
	}

	double deltaPhi(double phi1, double phi2) {
		double dphi = phi1 - phi2;
		while (dphi > M_PI) {
			dphi -= 2.0 * M_PI;
		}
		while (dphi < -M_PI) {
			dphi += 2.0 * M_PI;
		}
		return dphi;
	}

	Moments computeMoments(const std::vector<int>& jet_particles, const std::vector<cluster_type>& history, int Nmax, int Mmax)
	{
		Moments M;
		M.real.assign(Nmax+1, std::vector<double>(Mmax+1, 0.0));
		M.imag.assign(Nmax+1, std::vector<double>(Mmax+1, 0.0));

		double totalPt = 0.0;
		double etaCenter = 0.0;
		double phiX = 0.0;
		double phiY = 0.0;
		for (int id : jet_particles) {
			double eta = history[id].p().eta();
			double phi = history[id].p().phi();
			double pt  = history[id].p().pt();
			
			totalPt += pt;
			etaCenter += pt * eta;
		       	phiX += pt * std::cos(phi);
			phiY += pt * std::sin(phi);
		}
		etaCenter /= totalPt;
		double phiCenter = std::atan2(phiY, phiX);

		double etaWidth2 = 0.0;
		for (int id : jet_particles) {
			const auto& p = history[id].p();
			double deta = p.eta() - etaCenter;
			etaWidth2 += p.pt() * deta * deta;
		}
		etaWidth2 /= totalPt;
		double etaWidth = std::sqrt(etaWidth2);
		
		// Prevent division by zero
		if (etaWidth < 1e-8) {
			etaWidth = 1.0;
		}

		for (int id : jet_particles) {
			const auto& p = history[id].p();

			double eta = (p.eta() - etaCenter) / etaWidth;
			double phi = deltaPhi(p.phi(), phiCenter);;
			double pt  = p.pt();
			
			double exp_eta[Nmax+1];
			
			for (int n = 0; n <= Nmax; ++n) {
				exp_eta[n] = hermitePsi(n, eta);
			}
			
			for (int m = 0; m <= Mmax; ++m) {
				double c = std::cos(m * phi);
				double s = std::sin(m * phi);
				for (int n = 0; n <= Nmax; ++n) {
					M.real[n][m] += pt * exp_eta[n] * c;
					M.imag[n][m] += pt * exp_eta[n] * s;
				}
			}
		}
		
		return M;
	}

        inline virtual bool analyze(const Event& evt) override
        {
                assert(evt.diffusionSequenceReady);
                assert(evt.simpleDiffusionJetsReady);

                const unsigned nJetNodes = evt.simpleDiffusionNodes.size();

                std::string evtnum = std::to_string(evt.number());

                assert(nJetNodes == evt.simpleDiffusionJets.size());
                const std::vector<cluster_type>& history = evt.diffusionSequence.clustHist();
                const unsigned size = history.size();

		std::vector<double> im;
		std::vector<double> imdist;
                for (unsigned i=0; i<size; ++i) {
                        std::vector<int> jet_particles;
                        getLeaves(i, history, jet_particles);
                        if (jet_particles.size() < 2) {
                                continue;
                        }

                        Moments moments = computeMoments(jet_particles, history, 5, 5);
			
			for (int n = 0; n <= 5; ++n) {
				for (int k = 0; k <= 5; ++k) {
					double invariant = std::hypot(moments.real[n][k], moments.imag[n][k]);
					im.push_back(invariant);
					imdist.push_back(history[i].dist());
				}
			}
                }

                cnpy::npy_save("npyarrays/im/" + evtnum + "_im.npy", &im[0], {im.size()}, "w");
                cnpy::npy_save("npyarrays/imdist/" + evtnum + "_imdist.npy", &imdist[0], {imdist.size()}, "w");

                im.clear();
                imdist.clear();

                return true;
        }
};

#endif // INVARIENTMOMENTS_HH_

