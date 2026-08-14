#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h" //after v09_44 release
//#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NumuCuts.h"
#include "TVector3.h"
#include "TRandom3.h"

#include <algorithm>
#include <vector>

#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "stdio.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TF1.h"

// ROOT
#include "TGraph2D.h"
#include "TMath.h"

using namespace ana;


double MU_LENGHT = 50.;
double PRO_KINETIC_ENERGY = 50.;
double CONTAINMENT_CUT = 10.;
double TRACKSCORE_HIGH = 0.5;
double TRACKSCORE_LOW = 0.4;

// OLD VERSION
const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* spill){
  //for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
  //  double min_time =-1; double max_time =-1;
  //  if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
  //  if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
  //  if(match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return true;} 
  //}
  return true;
});

namespace wiremod 
{
  static constexpr float kWireSpacing = 0.3f;

  static constexpr float kIntegral_A =  1.35696f;
  static constexpr float kIntegral_B =  0.0786976f;
  static constexpr float kIntegral_C = -24.1874f;

  inline bool WireModHitCut(float phi, float pitch, float integral, int plane)
  {
    if (pitch <= 0.f || plane < 0 || plane > 2) return true;

    const float cosG      = kWireSpacing / pitch;
    const float thetaXW   = std::abs(std::atan(std::cos(phi) / cosG));
    const float thetaDeg  = thetaXW * (180.f / M_PI);

    const float thr_integral = kIntegral_A * std::exp(kIntegral_B * thetaDeg) + kIntegral_C;
    if (integral <= thr_integral) return true;

    return false;
  }
}

bool isInContained (double x, double y, double z, double dist)
{
	bool containment = false;

    if (std::isnan(x) || std::isnan(y) || std::isnan(z)) return false;

	// Check if point is in the triangles.
	// The limits are defined as 5 cm away from the hypotenuse
	// and they are calculated from geometric information of
	// coordinates of the 26th wire from the first and last
	// one in the planes
	if (	y <  1.732007 * z - 1687.5114	||
			y > -1.732007 * z + 1640.6114	||
			y >  1.732007 * z + 1640.6114	||
			y < -1.732007 * z - 1687.5114	)
	{
		return containment;
	}

	// Check distance from edges
	if (    (   ( x < -61.94  - dist && x > -358.49 + dist ) ||
				( x >  61.94  + dist && x <  358.49 - dist )     ) &&
			(   ( y > -181.86 + dist && y <  134.96 - dist ) &&
				( z > -894.95 + dist && z <  894.95 - dist )     )		)
	{
		containment = true;
	}

	return containment;
}

bool isInFV (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
  //need to add a check to avoid having the vtx slice in one cryo and the end track in the other cryo

    //fiducial volume for the dangling cable 
  if(x>210.0 && y > 60.0 && z> 290.0 && z< 390.0 ) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
			( x >  61.94 + 25 && x <  358.49 - 25 )) &&
		  ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
		  ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
}

bool isInActive (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 && x > -358.49 ) ||
			( x >  61.94 && x <  358.49)) &&
		  ( ( y > -181.86 && y < 134.96)  &&
		    ( z > -894.95 && z < 894.95) ));
}

bool all_contained ( const caf::Proxy<caf::SRSlice>& islc ) { 

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
    {
        if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
        if((islc.reco.pfp[ipfp].trk.start.x*islc.vertex.x)<0){return false;} //not contained if they cross cryostats
        //if not contained return false
        if(!isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,CONTAINMENT_CUT)){return false;}
    }   
    
return true;
}


int find_truth_muon ( const caf::Proxy<caf::SRSlice>& islc, int dist_mucut) { 

    //Select muon as longest track
    double max_length=-1.0;
    int ipfp_mu = -1;
    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 RecoStart;

    if(islc.is_clear_cosmic)return -1;

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){

        if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;

        RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);

        if(islc.reco.pfp[ipfp].trackScore<TRACKSCORE_HIGH)continue;

        if(islc.reco.pfp[ipfp].trk.len>max_length && ((RecoVtx-RecoStart).Mag()<dist_mucut) && islc.reco.pfp[ipfp].trk.len>MU_LENGHT && 
        std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==13 && 
        isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,CONTAINMENT_CUT) && 
        isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,CONTAINMENT_CUT) &&
        (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && islc.reco.pfp[ipfp].parent_is_primary){
        max_length=islc.reco.pfp[ipfp].trk.len;
        ipfp_mu=ipfp;
            }
        }//loop of pfp to find muon
return ipfp_mu;
}


bool is_truth_muon ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_mucut) 
{ 

    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 RecoStart;

    if(islc.is_clear_cosmic)return -1;



        if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return false;

        RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);

        if(islc.reco.pfp[ipfp].trackScore<TRACKSCORE_HIGH)return false;

        if( ((RecoVtx-RecoStart).Mag()<dist_mucut) && 
            islc.reco.pfp[ipfp].trk.len>MU_LENGHT && 
            std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==13 && 
            isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,CONTAINMENT_CUT) && 
            isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,CONTAINMENT_CUT) &&
            (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && 
            islc.reco.pfp[ipfp].parent_is_primary
          )return true;

    return false;
}

int id_pfp_truth ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_cut, int shower_plane = 2 ) { 
    //return 1 PROTONS
    //return 2 PIONS
    //return 3 SHOWER
    //return 9 other -> nan, not primary, too far, below energy threshold... 

    if(
        !isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,CONTAINMENT_CUT) //||
        //!isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,CONTAINMENT_CUT)
      )return 9;

    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);

    if (!(islc.reco.pfp[ipfp].parent_is_primary ))return 9;
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
 
    TVector3 start(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
    TVector3 end(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);

    double min_dist = ((start-RecoVtx).Mag()< (end-RecoVtx).Mag() ? (start-RecoVtx).Mag() : (end-RecoVtx).Mag());

    if(min_dist>50.0) return 9;
    
    if(islc.reco.pfp[ipfp].trackScore>=TRACKSCORE_HIGH)
    {
      if (std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
      if (std::isnan(islc.reco.pfp[ipfp].trk.start.y) || std::isnan(islc.reco.pfp[ipfp].trk.start.z)|| std::isnan(islc.reco.pfp[ipfp].trk.end.y) || std::isnan(islc.reco.pfp[ipfp].trk.end.z) )return 9;
    
      //skip low energy tagged pions
      TVector3 Start_mom_v;
      if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);}
      if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(139.570,2)+pow(Start_mom_v.Mag()*1000,2))-139.570)>=25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 2;}}
      
      //skip low energy protons
      if(islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 ){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);}
      if(islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(938.3,2)+pow(Start_mom_v.Mag()*1000,2))-938.3)>=PRO_KINETIC_ENERGY){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 1;}}
            
    }
    if(islc.reco.pfp[ipfp].trackScore<TRACKSCORE_HIGH)
    {
      if(islc.reco.pfp[ipfp].trackScore>=TRACKSCORE_LOW && islc.reco.pfp[ipfp].trackScore<TRACKSCORE_HIGH && islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 )
      {
        TVector3 Start_mom_v2;
        Start_mom_v2.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
        if((sqrt(pow(938.3,2)+pow(Start_mom_v2.Mag()*1000,2))-938.3)>=PRO_KINETIC_ENERGY && ((RecoVtx-start).Mag()<dist_cut) && islc.reco.pfp[ipfp].parent_is_primary){return 1;}
      }
      if(!(islc.reco.pfp[ipfp].trackScore>=TRACKSCORE_LOW && islc.reco.pfp[ipfp].trackScore<TRACKSCORE_HIGH && islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 )){

      if(std::isnan(islc.reco.pfp[ipfp].shw.plane[shower_plane].energy))return 9;
      if(islc.reco.pfp[ipfp].shw.plane[shower_plane].energy*1000<25.0)return 9;
      if(islc.reco.pfp[ipfp].shw.plane[shower_plane].energy*1000>25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 3;}}} 
    }
    
    return 9;
}


bool all_contained_truth ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ) 
{ 
    //Check only those pfp that are visible 
                
                for ( auto const& ipart : islc.truth.prim ){
                    if ( ipart.G4ID < 0 )  continue;
                    if ( ipart.cryostat < 0 )  continue;
                    int iG4ID_parent;
                    double dep_E=0;   
                    //check if charged primaries are contained: 
                    if(abs(ipart.pdg)==13 || abs(ipart.pdg)==2212 || abs(ipart.pdg)==211 || abs(ipart.pdg)==11){
                        if(isInContained(ipart.end.x,ipart.end.y,ipart.end.z,CONTAINMENT_CUT)==false){return false;}
                    }                   
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
			  if ( itrue.cryostat < 0 )  continue;
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){
                        if(abs(itrue.pdg)==13 || abs(itrue.pdg)==2212 || abs(itrue.pdg)==211 || abs(itrue.pdg)==11){
			  if(itrue.end.x!=-9999 && itrue.end.y!=-9999 && itrue.end.z!=-9999){if(isInContained(itrue.end.x,itrue.end.y,itrue.end.z,CONTAINMENT_CUT)==false){return false;}}
                    }                             
                        }
                        }
                    
                        } 
 
                }  
return true;
}

inline bool same_g4id(unsigned int parent, int g4id)
{
    if (g4id < 0) return false;
    return parent == static_cast<unsigned int>(g4id);
}

bool is_charged_lepton_or_hadron(int pdg) 
{
  if(
    std::abs(pdg)==13 ||
    std::abs(pdg)==2212 || 
    std::abs(pdg)==211 ||
    std::abs(pdg)==11
  )
  return true;

  return false;
}

bool all_contained_mc(const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRTrueInteraction>& nu)
{
    for (const auto& ipart : nu.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        int pid = std::abs(ipart.pdg);

        // Only visible charged primaries
        if (is_charged_lepton_or_hadron(pid)) {

            if (!isInContained(ipart.end.x, ipart.end.y, ipart.end.z, CONTAINMENT_CUT))
            {
                //MyFile_truth << "this is not contained " << endl;
                return false;}
        }

        // Daughters
        for (const auto& d : sr->true_particles) {

            if (!same_g4id(d.parent, ipart.G4ID) || d.cryostat < 0) continue;

            int dpid = std::abs(d.pdg);

            if (is_charged_lepton_or_hadron(dpid)) {

                if (d.end.x == -9999 || d.end.y == -9999 || d.end.z == -9999) continue;

                if (!isInContained(d.end.x, d.end.y, d.end.z, CONTAINMENT_CUT))
                    {
                        //MyFile_truth << "Daughter is not contained " <<  TruthPIDToString(dpid) << endl; 
                    return false;}
            }
        }
    }

    return true;
}

bool classification_type_MC(const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRTrueInteraction>& nu)
{
    // ----------------- Sanity -----------------
    if (std::isnan(nu.position.x) ||
        std::isnan(nu.position.y) ||
        std::isnan(nu.position.z)) {
        return false;
    }

    // ----------------- Neutrino + volume -----------------
    if (std::abs(nu.pdg) != 14 || !nu.iscc) {
        return false;
    }

    if (!isInActive(nu.position.x,
                    nu.position.y,
                    nu.position.z)) {
        return false;
    }

    if (!isInFV(nu.position.x,
                nu.position.y,
                nu.position.z)) {
        return false;
    }

    // ----------------- Counters -----------------
    int    n_muons = 0;
    int    n_protons_above = 0;
    double muon_length = 0.0;

    const int use_plane = 2;

    // ----------------- Primary loop -----------------
    for (const auto& ipart : nu.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        int pid = std::abs(ipart.pdg); 
        double depE = ipart.plane[ipart.cryostat][use_plane].visE * 1000.0;

        // -------- Muons --------
        if (pid == 13) {
            ++n_muons;
            muon_length = ipart.length;
        }

        // -------- Charged pions --------
        if (pid == 211 && depE > 25) {
            return false;
        }

        // -------- Neutral pions --------
        if (pid == 111) {
            for (const auto& d : sr->true_particles) {
                if (same_g4id(d.parent, ipart.G4ID) &&
                    std::abs(d.pdg) == 22) {

                    double eg = d.plane[ipart.cryostat][use_plane].visE * 1000.0;

                    if (eg > 25) {
                        return false;
                    }
                }
            }
        }

        // -------- Photons --------
        if (pid == 22) {
            for (const auto& d : sr->true_particles){
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;
            }

                if (depE > 25) {
                    return false;
                }
            
        }

        // -------- Protons --------
        if (pid == 2212) {
            for (const auto& d : sr->true_particles){
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;
            }
                if (depE > PRO_KINETIC_ENERGY) {
                    ++n_protons_above;
                }
            
        }
    }

    // ----------------- Containment -----------------
    if (!all_contained_mc(sr, nu)) {
        return false;
    }

    // ----------------- Final classification -----------------

    if (n_muons == 1 && muon_length > 50) {

        if (n_protons_above == 1) {
                //MyFile_truth2 << "\n================ NEW EVENT ================\n";
                //MyFile_truth2 << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " 1" << endl; 

            return true;}

        if (n_protons_above > 1) {
                //MyFile_truth2 << "\n================ NEW EVENT ================\n";
                //MyFile_truth2 << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " " << n_protons_above << endl; 

            return true;}

        return false;
    }

    return false;
}


std::string classification_type_generic ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ) 
{

    // the output is a string with format abcde_process
    // a : # muons
    // b : # protons above 50 MeV
    // c : # pions above 25 MeV 
    // d : # gammas above 25 MeV
    // e : # pi0 with bot gammas above 25 MeV
    // bad_slice : slice does not satisfy : truth_matching_eff > 0.05 && is_clear_cosmic && true-reco vertex distance < 100 or Nan values

    std::string return_string="bad_slice";
    
    if(islc.truth.index<0)return "cosmic";

    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);

    //if( !(islc.tmatch.eff > 0.05 && (!islc.is_clear_cosmic) && (vertex_true-vertex_reco).Mag()<100.)) return "bad_slice";

    if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)) return "reco_vertex_isNAN";
    
    if(std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)) return "true_vertex_isNAN";

    //if (  abs(islc.truth.pdg) == 14 /*&& islc.truth.iscc && isInActive(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z) && islc.truth.index<0*/ )
    if(abs(islc.truth.pdg) != 14)return "not_numu"; 
          //{          
                if(!isInActive(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z))return "not_in_Active";

                if(!isInFV(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z))return "not_in_FV";
            //{
                if(!all_contained_truth(sr, islc))return "not_contained";
                
                int num_protons_above50 = 0;
                int num_pions_above25 = 0;
                int num_gammas_above25 = 0;
                int num_neutral_pions_both_above25 = 0; 
                int num_muons = 0; 
                //int interaction_type = (int)islc.truth.genie_inttype;
                             
                double dep_E=0;
                for ( auto const& ipart : islc.truth.prim )
                {
                  if ( ipart.G4ID < 0 )  continue;

                  //MUONS
                  if(abs(ipart.pdg) == 13 && ipart.length > 50){num_muons+=1;}  // muons

                  int iG4ID_parent;  
                  int use_plane = 2;
                    
                  //PRIMARY NEUTRAL PIONS
                  if(abs(ipart.pdg)==111) //Neutral pions - reject if any of its gamma > 25 MeV
                  {            
                    if(ipart.daughters.size()>0)
                    {
                      for ( auto const& itrue2 : sr->true_particles )
                      {
                        iG4ID_parent=itrue2.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID && abs(itrue2.pdg) == 22)
                        {
                          if(itrue2.plane[ipart.cryostat][use_plane].visE*1000>25)
                          {
                            num_neutral_pions_both_above25++;
                          }
                        }
                      }
                    }                        
                  }

                  //PRIMARY PHOTONS
                  if(abs(ipart.pdg) == 22)
                  {                    
                    if(ipart.daughters.size()>0)
                    {
                      for ( auto const& itrue : sr->true_particles )
                      {
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID )
                        {
                          dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }
                      }
                    }
                    dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                  } 
                  if(abs(ipart.pdg)== 22 && dep_E>25.0){num_gammas_above25++;}   
                  dep_E=0;

                  //PRIMARY PROTONS
                  if(abs(ipart.pdg)== 2212)
                  {                    
                    if(ipart.daughters.size()>0)
                    {
                      for ( auto const& itrue : sr->true_particles )
                      {
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID )
                        {
                          dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }
                      }
                    }
                    dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                  }
                  if(abs(ipart.pdg)== 2212 && dep_E>PRO_KINETIC_ENERGY){num_protons_above50+=1;} //protons
                  dep_E=0;  

                  //PRIMARY CHARGED PIONS
                  /*
                  if(abs(ipart.pdg)== 211)
                  {                    
                    if(ipart.daughters.size()>0)
                    {
                      for ( auto const& itrue : sr->true_particles )
                      {
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID )
                        {
                          dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }
                      }
                    }
                    dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                  }
                  if(abs(ipart.pdg)== 211 && dep_E>25.0){num_pions_above25+=1;} //pions
                  dep_E=0;
                  */
                  if(abs(ipart.pdg)== 211)
                  {                    
                    dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                  }
                  if(abs(ipart.pdg)== 211 && dep_E>25.0){num_pions_above25+=1;} //pions
                  dep_E=0;

                }//all true particles in the slice

                //slices that satisfy: non Nan true and reco vertex, in Active, In Fiducial, All contained, nu mu interaction
                return_string = Form("%d%d%d%d%d",num_muons,num_protons_above50,num_pions_above25,num_gammas_above25,num_neutral_pions_both_above25);

              //}//fiducial
            //}//muon neutrino
         
    return return_string; 
}


double bar_flash (const caf::SRSpillProxy* spill,const caf::Proxy<caf::SRSlice>& islc)
{
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if((islc.vertex.x*match.flashPosition.x>0) && match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return match.flashPosition.z;} 

  }
  return -10000;
};

double bar_flash_x (const caf::SRSpillProxy* spill,const caf::Proxy<caf::SRSlice>& islc)
{
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if((islc.vertex.x*match.flashPosition.x>0) && match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return match.flashPosition.x;} 

  }
  return 0;
};
