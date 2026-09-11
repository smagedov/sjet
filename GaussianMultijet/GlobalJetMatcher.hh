#ifndef GLOBALJETMATCHER_HH_
#define GLOBALJETMATCHER_HH_

#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <stdexcept>

namespace GlobalJetMatching
{

struct MatchCandidate
{
    unsigned cluster;
    double score;
    double deltaR;
    double ptRatio;
    double oracleMass;
    double clusterMass;
    double massPenalty;
};

struct MatchResult
{
    bool found = false;

    double totalScore = std::numeric_limits<double>::infinity();

    double maxRegret = std::numeric_limits<double>::infinity();

    double totalRegret = std::numeric_limits<double>::infinity();

    std::vector<unsigned> clusterForOracle;

    // Diagnostics
    std::vector<bool> protectedOracle;
    std::vector<double> bestLocalScore;
};


namespace detail
{
inline bool goodLocalMatch(
		const MatchCandidate& c,
		double maxDeltaR,
		double maxAbsLogPtRatio)
{
	return c.deltaR < maxDeltaR && std::abs(std::log(c.ptRatio)) < maxAbsLogPtRatio;
}


inline bool betterSolution(
    double maxRegret,
    double totalRegret,
    double totalScore,
    const MatchResult& bestResult)
{
    constexpr double eps = 1.0e-12;

    if (!bestResult.found)
        return true;

    // Primary:
    // minimize worst sacrifice of any oracle jet.
    if (maxRegret < bestResult.maxRegret - eps)
        return true;

    if (maxRegret > bestResult.maxRegret + eps)
        return false;

    // Secondary:
    // minimize total sacrifice.
    if (totalRegret < bestResult.totalRegret - eps)
        return true;

    if (totalRegret > bestResult.totalRegret + eps)
        return false;

    // Final tie-breaker:
    // minimize original matching score.
    return totalScore < bestResult.totalScore - eps;
}


inline void searchGlobalMatches(
    unsigned depth,
    const std::vector<unsigned>& oracleOrder,
    const std::vector<std::vector<MatchCandidate>>& candidates,
    const std::vector<std::vector<bool>>& incompatible,
    const std::vector<double>& bestLocalScore,
    const std::vector<bool>& protectedOracle,
    double protectedMaxRegret,
    std::vector<unsigned>& currentAssignment,
    std::vector<unsigned>& chosenClusters,
    double currentMaxRegret,
    double currentTotalRegret,
    double currentTotalScore,
    MatchResult& bestResult,
    unsigned invalidCluster)
{
    // ============================================================
    // Complete assignment.
    // ============================================================

    if (depth == oracleOrder.size())    {
        if (betterSolution(currentMaxRegret, currentTotalRegret, currentTotalScore, bestResult)) {
            bestResult.found = true;

            bestResult.maxRegret = currentMaxRegret;

            bestResult.totalRegret = currentTotalRegret;

            bestResult.totalScore = currentTotalScore;

            bestResult.clusterForOracle = currentAssignment;
        }

        return;
    }


    // ============================================================
    // Branch-and-bound.
    // ============================================================

    if (bestResult.found) {
        constexpr double eps = 1.0e-12;

        if (currentMaxRegret > bestResult.maxRegret + eps) {
            return;
        }

        if (std::abs(currentMaxRegret - bestResult.maxRegret) <= eps && currentTotalRegret > bestResult.totalRegret + eps) {
            return;
        }
    }


    const unsigned oracle = oracleOrder[depth];


    // ============================================================
    // Try candidates for this oracle jet.
    // ============================================================

    for (const MatchCandidate& candidate : candidates[oracle]) {
        const unsigned cluster = candidate.cluster;


        // --------------------------------------------------------
        // Calculate local regret.
        // --------------------------------------------------------

        double regret = candidate.score - bestLocalScore[oracle];

        if (regret < 0.0 && regret > -1.0e-12) {
            regret = 0.0;
        }


        // --------------------------------------------------------
        // Protect jets that already have a genuinely good local
        // reconstruction.
        //
        // These jets are not allowed to degrade arbitrarily just
        // to rescue a bad oracle jet elsewhere in the event.
        // --------------------------------------------------------

        if (protectedOracle[oracle] && regret > protectedMaxRegret) {
            continue;
        }


        // --------------------------------------------------------
        // Check compatibility with previously selected nodes.
        // --------------------------------------------------------

        bool compatible = true;

        for (unsigned chosen : chosenClusters) {
            if (incompatible[cluster][chosen]) {
                compatible = false;
                break;
            }
        }

        if (!compatible)
            continue;


        const double newMaxRegret = std::max(currentMaxRegret, regret);

        const double newTotalRegret = currentTotalRegret + regret;

        const double newTotalScore = currentTotalScore + candidate.score;


        // --------------------------------------------------------
        // Additional branch-and-bound.
        // --------------------------------------------------------

        if (bestResult.found)
        {
            constexpr double eps = 1.0e-12;

            if ( newMaxRegret > bestResult.maxRegret + eps) {
                continue;
            }

            if (std::abs(newMaxRegret - bestResult.maxRegret) <= eps && newTotalRegret > bestResult.totalRegret + eps) {
                continue;
            }
        }


        currentAssignment[oracle] = cluster;

        chosenClusters.push_back(cluster);


        searchGlobalMatches(
            depth + 1,
            oracleOrder,
            candidates,
            incompatible,
            bestLocalScore,
            protectedOracle,
            protectedMaxRegret,
            currentAssignment,
            chosenClusters,
            newMaxRegret,
            newTotalRegret,
            newTotalScore,
            bestResult,
            invalidCluster);


        chosenClusters.pop_back();

        currentAssignment[oracle] = invalidCluster;
    }
}

} // namespace detail


template <
    typename OracleJetContainer,
    typename HistoryContainer,
    typename DistanceFunctor,
    typename IncompatibilityMarker> MatchResult match(
		    const OracleJetContainer& oracleJets,
		    const HistoryContainer& history,
		    unsigned nJets,
		    double alpha,
		    double beta,
		    DistanceFunctor distance,
		    IncompatibilityMarker markIncompatible,
		    
		    // Number of candidates retained for protected jets.
		    unsigned maxCandidates = 100,
		    
		    bool printCandidates = false,
		    
		    // Definition of a genuinely good local match.
		    double goodMatchMaxDeltaR = 0.10,
		    double goodMatchMaxAbsLogPtRatio = 0.10,
		    // Maximum allowed degradation for protected jets.
		    double protectedMaxRegret = std::numeric_limits<double>::infinity()) 
{
    const unsigned historySize = static_cast<unsigned>(history.size());


    if (nJets > oracleJets.size()) {
        throw std::runtime_error("GlobalJetMatcher: nJets exceeds oracle jet collection size.");
    }


    MatchResult result;

    result.clusterForOracle.assign(nJets, historySize);

    result.protectedOracle.assign(nJets, false);

    result.bestLocalScore.assign(nJets, std::numeric_limits<double>::infinity());


    if (nJets == 0) {
        result.found = true;

        result.totalScore = 0.0;
        result.maxRegret = 0.0;
        result.totalRegret = 0.0;

        return result;
    }


    if (historySize == 0)
        return result;



    // ============================================================
    // Build incompatibility matrix.
    // ============================================================

    std::vector<std::vector<bool>> incompatible(historySize, std::vector<bool>(historySize, false));


    for (unsigned j = 0; j < historySize; ++j) {
        std::vector<bool> available(historySize, true);


        markIncompatible(j, available);


        for (unsigned k = 0; k < historySize; ++k) {
            if (!available[k]) {
                incompatible[j][k] = true;
                incompatible[k][j] = true;
            }
        }


        // Do not allow exact same node twice.
        incompatible[j][j] = true;
    }



    // ============================================================
    // Generate FULL candidate lists first.
    //
    // Important:
    // We need to determine whether each jet has a good local match
    // before deciding whether its candidate list should be trimmed.
    // ============================================================

    std::vector<std::vector<MatchCandidate>> candidates(nJets);


    for (unsigned i = 0; i < nJets; ++i) {
        const double oraclePt = oracleJets[i].pt();


        if (oraclePt <= 0.0)
            continue;


        candidates[i].reserve(historySize);


        for (unsigned j = 0; j < historySize; ++j) {
            const double clusterPt = history[j].p().pt();


            if (clusterPt <= 0.0)
                continue;


            const double deltaR = distance(oracleJets[i], history[j].p());
            const double ptRatio = oraclePt / clusterPt;
	    const double oracleMass = std::max(0.0, oracleJets[i].m());
	    const double clusterMass = std::max(0.0, history[j].p().m());
	    
	    // Prevent problems for nearly massless jets.
	    constexpr double massEpsilon = 1.0;
	    // Only penalize candidates that are MORE massive than the corresponding oracle jet.
	    const double massPenalty = std::abs(std::log((clusterMass + massEpsilon) / (oracleMass + massEpsilon)));


            const double score = deltaR + alpha*std::abs(std::log(ptRatio)) + beta*massPenalty;


            candidates[i].push_back({j, score, deltaR, ptRatio, oracleMass, clusterMass, massPenalty});
        }


        std::sort(candidates[i].begin(), candidates[i].end(),
            [](const MatchCandidate& a, const MatchCandidate& b) {
                return a.score < b.score;
            });
    }



    // ============================================================
    // Every oracle jet needs at least one candidate.
    // ============================================================

    for (unsigned i = 0; i < nJets; ++i) {
        if (candidates[i].empty())
            return result;
    }



    // ============================================================
    // Determine best local score and protection status.
    // ============================================================

    for (unsigned i = 0; i < nJets; ++i) {
        result.bestLocalScore[i] = candidates[i][0].score;

        result.protectedOracle[i] = detail::goodLocalMatch(candidates[i][0], goodMatchMaxDeltaR, goodMatchMaxAbsLogPtRatio);
    }



    // ============================================================
    // Candidate trimming.
    //
    // Protected jet:
    //
    //   its local reconstruction is already good.
    //   Keep maxCandidates for efficiency.
    //
    // Unprotected jet:
    //
    //   its local reconstruction is already poor.
    //   KEEP THE ENTIRE HISTORY.
    //
    // This allows the bad jet to move very far down its candidate
    // list rather than destroying a good jet.
    // ============================================================

    for (unsigned i = 0; i < nJets;++i) {
        if (result.protectedOracle[i] && candidates[i].size() > maxCandidates) {
            candidates[i].resize(maxCandidates);
        }
    }



    // ============================================================
    // Candidate diagnostics.
    // ============================================================

    if (printCandidates) {
        for (unsigned i = 0; i < nJets; ++i) {
            std::cout << "\nOracle jet " << i << " candidate matches:";

            if (result.protectedOracle[i]) {
                std::cout << " [PROTECTED]";
            } else {
                std::cout << " [UNPROTECTED]";
            }

            std::cout << std::endl;


            const unsigned nPrint = std::min<unsigned>(5, candidates[i].size());


            for (unsigned k = 0; k < nPrint; ++k) {
                const MatchCandidate& c = candidates[i][k];


                const double regret = c.score - result.bestLocalScore[i];


                std::cout
		<< "  #" << k
    		<< " cluster " << c.cluster
    		<< " score " << c.score
    		<< " regret " << regret
    		<< " dR " << c.deltaR
    		<< " pT ratio " << c.ptRatio
    		<< " oracle mass " << c.oracleMass
    		<< " cluster mass " << c.clusterMass
    		<< " mass penalty " << c.massPenalty
    		<< std::endl;
            }
        }


        std::cout << std::endl;
    }




    // ============================================================
    // Determine recursive search order.
    //
    // Protected jets first.
    //
    // Among protected jets, use best local score first.
    //
    // Unprotected jets are handled last so they adapt around the
    // already-good reconstruction.
    // ============================================================

    std::vector<unsigned> oracleOrder( nJets);

    std::iota(oracleOrder.begin(), oracleOrder.end(), 0);

    std::sort(
        oracleOrder.begin(),
        oracleOrder.end(),
        [&](unsigned a, unsigned b) {
            if (result.protectedOracle[a] != result.protectedOracle[b]) {
                return result.protectedOracle[a] > result.protectedOracle[b];
            }
	    return result.bestLocalScore[a] < result.bestLocalScore[b];
        });



    // ============================================================
    // Run protected global search.
    // ============================================================

    std::vector<unsigned> currentAssignment(nJets, historySize);
    std::vector<unsigned> chosenClusters;
    chosenClusters.reserve(nJets);


    detail::searchGlobalMatches(
        0,
        oracleOrder,
        candidates,
        incompatible,
        result.bestLocalScore,
        result.protectedOracle,
        protectedMaxRegret,
        currentAssignment,
        chosenClusters,
        0.0,// currentMaxRegret
        0.0,// currentTotalRegret
        0.0,// currentTotalScore
        result,
        historySize);



    // ============================================================
    // Diagnostics.
    // ============================================================

    if (printCandidates) {
        if (result.found) {
            std::cout << "GLOBAL MATCH:" << std::endl;

            std::cout << "  max regret: "  << result.maxRegret << std::endl;

            std::cout << "  total regret: " << result.totalRegret << std::endl;

            std::cout << "  total score: " << result.totalScore << std::endl;

            for (unsigned i = 0; i < nJets; ++i) {
                const unsigned selectedCluster = result.clusterForOracle[i];


                bool printed = false;


                for (const MatchCandidate& c : candidates[i]) {
                    if (c.cluster != selectedCluster) {
                        continue;
                    }


                    const double regret = c.score - result.bestLocalScore[i];


                    std::cout << "  oracle " << i;

                    if (result.protectedOracle[i]) {
                        std::cout << " [PROTECTED]";
                    } else {
                        std::cout << " [UNPROTECTED]";
                    }


                    std::cout
			<< " -> cluster " << selectedCluster
    			<< " score " << c.score
    			<< " regret " << regret
    			<< " dR " << c.deltaR
    			<< " pT ratio " << c.ptRatio
    			<< std::endl;

                    printed = true;
                    break;
                }


                if (!printed) {
                    std::cout << "  oracle " << i  << " -> cluster " << selectedCluster << " [candidate diagnostics unavailable]" << std::endl;
                }
            }


            std::cout << std::endl;
        } else {
            std::cout << "GLOBAL MATCH: " << "no complete assignment exists while " << "respecting protected-match constraints." << std::endl;

            std::cout << "  protectedMaxRegret = " << protectedMaxRegret << std::endl;

            std::cout << std::endl;
        }
    }


    return result;
}


} // namespace GlobalJetMatching

#endif
