#ifndef DiffusionDistFastJet_hh
#define DiffusionDistFastJet_hh

//////////////////////////////////////////
/// Wrapper to interface the diffusion distance
/// jet algorithm and fast jet
//////////////////////////////////////////

#include <string>
#include <sstream>
#include <cfloat>
#include <iostream>

#include "diffusionFactor2.h"

#include "fastjet/NNH.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;

namespace stab {

class DiffBriefJet {
public:
  void init(const fastjet::PseudoJet &jet ) {
    std::cout << "We are in init." << std::endl;
    pt = jet.pt();
    eta = jet.eta();
    phi = jet.phi();
  }

  double distance( const DiffBriefJet *jt ) const {
    std::cout << "We are in distance." << std::endl;
    double dEta = eta - jt->eta;
    double dPhi = deltaPhi(phi,jt->phi);
    double ratio = ( pt < jt->pt ? pt : jt->pt ) / ( pt < jt->pt ? jt->pt : pt );
    return sqrt(dEta*dEta+dPhi*dPhi) * diffusionFactor2(ratio);
  }

  double beam_distance() const {
    return DBL_MAX;
  }

private:
  double pt,eta,phi;

  inline double deltaPhi(const double phi1, const double phi2) const {
    double dPhi = phi1 - phi2;
    while ( dPhi < -M_PI )
      dPhi += 2.*M_PI;
    while ( dPhi > M_PI )
      dPhi -= 2.*M_PI;

    return dPhi;
  }

};

void DiffusionDistFastJet::run_clustering(ClusterSequence &seq) const
{
  // This function should fill the rest of the ClusterSequence using it member function
  // plugin_do_ij_recombination()
  // plugin_do_iB_recombination()

  // follow the example of the EECambridgePlugin included in the FastJet package
  std::cout << "maybe?" << std::endl;
  int njets = seq.jets().size();
  std::cout << njets << std::endl;
  fastjet::NNH<DiffBriefJet> nnh(seq.jets());

  // the ClusterSequence has already been filled with the initial event
  const std::vector<PseudoJet> & initialEvent = seq.jets();

  double dij=-1;
  while ( njets > 1 ) {
    int i,j,k;
    dij = nnh.dij_min(i,j); // i,j are return values

    if ( j >=0 ) {
      seq.plugin_record_ij_recombination(i,j,dij,k);
      nnh.merge_jets(i,j,seq.jets()[k],k);
      njets--;
    } else {
      seq.plugin_record_iB_recombination(i,dij);
      nnh.remove_jet(i);
    }

  }

  const std::vector<fastjet::ClusterSequence::history_element> & history = seq.history();
  const unsigned histSize = history.size();
  for ( unsigned h=0; h != histSize; h++ ) {
    const fastjet::ClusterSequence::history_element & el = history[h];

    // if child node is < 0, marke as final state jet by doing a beam merge
    if ( el.child < 0 ) {
      seq.plugin_record_iB_recombination(h,DBL_MAX);
      nnh.remove_jet(h);
    }

  }


};

class DiffusionDistFastJet: public JetDefinition::Plugin {

public:
  inline DiffusionDistFastJet(const double R): JetDefinition::Plugin()
    ,m_R(R)
    { }

  inline virtual ~DiffusionDistFastJet() { }

  // returns a textual description of the jet algorithm
  inline virtual std::string description() const {

    std::stringstream desc;
    desc << "Diffusion Distance (Fast Jet Plugin) Description:\n";
    desc << "The diffusion distance is defined as the minimum width Gaussian filter\n";
    desc << "for which only one energy deposit is observed from the original two.\n";
    desc << "It is analogous to the sqrt(t) at which two drops of honey\n";
    desc << "diffuse into a one."; 
    
    return desc.str();
  }

  inline virtual double R() const { return m_R; }

  // run the clustering
  virtual void run_clustering(ClusterSequence &) const;


private:

    // how the DiffusionDist class calculates the distance
    //inline double operator()(const jet &a, const jet &b ) const
    //{
    //   double pT1 = a.pt();
    //   double pT2 = b.pt();
    //   double ratio =  (pT1 < pT2 ? pT1 : pT2 )/(pT1 < pT2 ? pT2 : pT1 );
    //   return geom3::deltaR(a, b) * diffusionFactor2(ratio);
    //}

  double m_R;

};



} // end stab namespace


#endif
