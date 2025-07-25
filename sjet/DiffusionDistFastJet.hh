#ifndef DiffusionDistFastJet_hh
#define DiffusionDistFastJet_hh

//////////////////////////////////////////
/// Wrapper to interface the diffusion distance
/// jet algorithm and fast jet
//////////////////////////////////////////

#include <string>
#include <sstream>

#include "diffusionFactor2.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;

namespace stab {


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
