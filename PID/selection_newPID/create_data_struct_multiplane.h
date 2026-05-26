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
#include "TTree.h"
#include "TBranch.h"

#include "data_struct.h" //--> definition of data struct
#include "/exp/icarus/app/users/nsommagg/NicolaICARUS/PID/selection1muNp_PIDfinale/helper_newPID.h" //--> PID helper functions

using namespace ana;

//TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
//auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
//auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
//auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
//auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

//FNAL
TFile* file = TFile::Open("/exp/icarus/data/users/nsommagg/RefCurvesChi2.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

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

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
    //if (!(islc.reco.pfp[ipfp].parent_is_primary ))continue; //skip secondaries
    //if(islc.reco.pfp[ipfp].trackScore<0.4)continue; //Want to check only tracks??
    if((islc.reco.pfp[ipfp].trk.start.x*islc.vertex.x)<0){return false;} //not contained if they cross cryostats
    //if not contained return false
    if(!isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0)){return false;}
      
    }   
    
return true;
}

bool all_contained_truth ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ) { 
    //Check only those pfp that are visible 
                
                for ( auto const& ipart : islc.truth.prim ){
                    if ( ipart.G4ID < 0 )  continue;
                    if ( ipart.cryostat < 0 )  continue;
                    int iG4ID_parent;
                    double dep_E=0;   
                    //check if charged primaries are contained: 
                    if(abs(ipart.pdg)==13 || abs(ipart.pdg)==2212 || abs(ipart.pdg)==211 || abs(ipart.pdg)==11){
                        if(isInContained(ipart.end.x,ipart.end.y,ipart.end.z,5.)==false){return false;}
                    }                   
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
			  if ( itrue.cryostat < 0 )  continue;
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){
                        if(abs(itrue.pdg)==13 || abs(itrue.pdg)==2212 || abs(itrue.pdg)==211 || abs(itrue.pdg)==11){
			  if(itrue.end.x!=-9999 && itrue.end.y!=-9999 && itrue.end.z!=-9999){if(isInContained(itrue.end.x,itrue.end.y,itrue.end.z,5.)==false){return false;}}
                    }                             
                        }
                        }
                    
                        } 
 
                }  
return true;
}


std::vector<double> chi2_ALG(std::vector<double> &dEdx,std::vector<double> &RR, double rr_min, double rr_max)
{
    //The output is chi2s

    ////////////////// correzzione dEdx mc ///////////////////////////////////////////
    /*
  correction_function->SetParameters(-1.08107e+00, -1.38302e+00, -2.84635e-02, 5.87297e-02, 7.17814e+00, -3.51461e+00);
    for(int i=0; i<int(dEdx.size()); i++)
    {
        double factor = (2+correction_function->Eval(dEdx[i]))/(2-correction_function->Eval(dEdx[i]));
        dEdx[i]=dEdx[i]*factor;
    }
    */
    //////////////////////////////////////////////////////////////////////////////////////

    double threshold = 0.5;
    double max_rr = rr_max;
    double min_rr = rr_min;

    std::vector<float> trkdedx;
    std::vector<float> trkres;
    std::vector<double> vpida;
    

    for(std::size_t i(0); i<dEdx.size(); ++i){
      if(i==0 || i==dEdx.size()-1)continue;
        if(RR[i]<max_rr && RR[i]>rr_min ){trkdedx.push_back(dEdx[i]);trkres.push_back(RR[i]);}       
    }


    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double avgdedx = 0;
    double PIDA = 0;

    int used_trkres = 0;
    for (unsigned i = 0; i < trkdedx.size(); ++i) { 
      avgdedx += trkdedx[i];
      if (trkres[i] < 26) {
        PIDA += trkdedx[i] * pow(trkres[i], 0.42);
        vpida.push_back(trkdedx[i] * pow(trkres[i], 0.42));
        used_trkres++;
      }
      if (trkdedx[i] > 100 || trkdedx[i]<threshold) continue; //protect against large pulse height
      int bin = dedx_range_pro->FindBin(trkres[i]);
      if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
        double bincpro = dedx_range_pro->GetBinContent(bin);
        if (bincpro < 1e-6) { //for 0 bin content, using neighboring bins
          bincpro =
            (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
        }
        double bincka = dedx_range_ka->GetBinContent(bin);
        if (bincka < 1e-6) {
          bincka =
            (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;
        }
        double bincpi = dedx_range_pi->GetBinContent(bin);
        if (bincpi < 1e-6) {
          bincpi =
            (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
        }
        double bincmu = dedx_range_mu->GetBinContent(bin);
        if (bincmu < 1e-6) {
          bincmu =
            (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
        }
        double binepro = dedx_range_pro->GetBinError(bin);
        if (binepro < 1e-6) {
          binepro =
            (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;
        }
        double bineka = dedx_range_ka->GetBinError(bin);
        if (bineka < 1e-6) {
          bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;
        }
        double binepi = dedx_range_pi->GetBinError(bin);
        if (binepi < 1e-6) {
          binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;
        }
        double binemu = dedx_range_mu->GetBinError(bin);
        if (binemu < 1e-6) {
          binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;
        }
        //double errke = 0.05*trkdedx[i];   //5% KE resolution
        double errdedx = 0.04231 + 0.0001783 * trkdedx[i] * trkdedx[i]; //resolution on dE/dx
        errdedx *= trkdedx[i];
        chi2pro += pow((trkdedx[i] - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
        chi2ka += pow((trkdedx[i] - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
        chi2pi += pow((trkdedx[i] - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
        chi2mu += pow((trkdedx[i] - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
        //std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
        ++npt;
      }
    } //hits
        std::vector<double> chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt};

    return chi2s;
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

bool all_contained_mc(const caf::SRSpillProxy* sr,
                         const caf::Proxy<caf::SRTrueInteraction>& nu)
{
    for (const auto& ipart : nu.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        int pid = std::abs(ipart.pdg);

        // Only visible charged primaries
        if (is_charged_lepton_or_hadron(pid)) {

            if (!isInContained(ipart.end.x, ipart.end.y, ipart.end.z, 5))
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

                if (!isInContained(d.end.x, d.end.y, d.end.z, 5))
                    {
                        //MyFile_truth << "Daughter is not contained " <<  TruthPIDToString(dpid) << endl; 
                    return false;}
            }
        }
    }

    return true;
}

bool classification_type_MC(const caf::SRSpillProxy* sr,
                                      const caf::Proxy<caf::SRTrueInteraction>& nu)
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
                if (depE > 50) {
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

std::string classification_type_generic ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ) {

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
                  if(abs(ipart.pdg)== 2212 && dep_E>50.0){num_protons_above50+=1;} //protons
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

double bar_flash (const caf::SRSpillProxy* spill,const caf::Proxy<caf::SRSlice>& islc){
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if((islc.vertex.x*match.flashPosition.x>0) && match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return match.flashPosition.z;} 

  }
  return -10000;
};

double bar_flash_x (const caf::SRSpillProxy* spill,const caf::Proxy<caf::SRSlice>& islc){
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if((islc.vertex.x*match.flashPosition.x>0) && match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return match.flashPosition.x;} 

  }
  return 0;
};

////////////////////////////////////////
////// HELPER FUNCTIONS NEW PID   /////
///////////////////////////////////////

int find_muon_newPID ( const caf::Proxy<caf::SRSlice>& islc, int dist_mucut, int mu1_pfp = -1, int mu2_pfp = -1, double stitched_len = -1) { 

    int ipfp_mu = -1;

    double max_p_as_mu = -1;

    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 RecoStart;

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
    {
      if((int)ipfp == mu2_pfp)continue; //--> in this first part we care only about the start of the track

      if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
        
      RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
        
      bool is_primary = islc.reco.pfp[ipfp].parent_is_primary;

      if(islc.reco.pfp[ipfp].trackScore<0.5)continue;
   
      double track_length = -1;
      int current_pfp = ipfp;
      if(current_pfp == mu1_pfp && mu2_pfp>=0) //--> now we care about the last part of the track
      {
        track_length = stitched_len; 
        current_pfp = mu2_pfp;
      }
      else{track_length = islc.reco.pfp[ipfp].trk.len;}

      //PID --> look for the most mu-like pfp in the slice that it is not a proton

      int pid_pred = PIDclass(islc,current_pfp);
      if(pid_pred == 2 || pid_pred == 3)continue;

      std::vector<double> prediction_proba;
      prediction_proba = PIDproba(islc,current_pfp); //probabilities for the track to be each of the 6 classes

      double imu_prob = std::max(prediction_proba[0], prediction_proba[1]);
 
      if(imu_prob>max_p_as_mu && 
        ((RecoVtx-RecoStart).Mag()<dist_mucut) && 
        track_length > 50. && 
        isInContained(islc.reco.pfp[current_pfp].trk.end.x,islc.reco.pfp[current_pfp].trk.end.y,islc.reco.pfp[current_pfp].trk.end.z,5.0) && 
        (islc.reco.pfp[current_pfp].trk.end.x*islc.vertex.x)>0 &&
        is_primary)
        {
          max_p_as_mu=imu_prob;
          ipfp_mu=ipfp;
        }
      }//loop of pfp to find muon

      return ipfp_mu;
}



int id_pfp_newPID ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_cut) { 
    //return 1 PROTONS
    //return 2 PIONS
    //return 3 SHOWER
    //return 9 other -> nan, not primary, too far, below energy threshold... 


    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);

    //skip secondaries
    if (!(islc.reco.pfp[ipfp].parent_is_primary ))return 9;

    //check for non nan values
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.start.y) || std::isnan(islc.reco.pfp[ipfp].trk.start.z)|| std::isnan(islc.reco.pfp[ipfp].trk.end.y) || std::isnan(islc.reco.pfp[ipfp].trk.end.z) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;

    TVector3 start(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
    TVector3 end(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);

    //start - vertex distance
    double min_dist = ((start-RecoVtx).Mag()< (end-RecoVtx).Mag() ? (start-RecoVtx).Mag() : (end-RecoVtx).Mag());

    if(min_dist>50.0) return 9;
    
    //PID
    int PIDprediction = -1;
    PIDprediction = PIDclass(islc,ipfp); //it gives you the predicted label for the track

    TVector3 Start_mom_v;
    //PIONS
    if(islc.reco.pfp[ipfp].trackScore>=0.5)
    {
      TVector3 Start_mom_v;
      if(PIDprediction == 0 || PIDprediction == 1 || PIDprediction == 4 || PIDprediction == 5 )
      { 
        Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);

        if( ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(139.570,2)+pow(Start_mom_v.Mag()*1000,2))-139.570)>=25.0 )return 2;
      }
    }
    //PROTONS
    if(islc.reco.pfp[ipfp].trackScore>=0.4)
    {
      if(PIDprediction == 2 || PIDprediction == 3 )
      { 
        Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);

        if( ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(938.3,2)+pow(Start_mom_v.Mag()*1000,2))-938.3)>=50.0 )return 1;
      }       
    }

    //SHOWERS
    if(islc.reco.pfp[ipfp].trackScore<0.5)
    {
      if(std::isnan(islc.reco.pfp[ipfp].shw.plane[2].energy))return 9;
      if(islc.reco.pfp[ipfp].shw.plane[2].energy*1000<25.0)return 9;
      if(islc.reco.pfp[ipfp].shw.plane[2].energy*1000>25.0)return 3;
    }

    return 9;
}

std::string selection_1muNp_newPID ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_cut, int cut_baryc, int mu1_pfp = -1, int mu2_pfp = -1, double stitched_len = -1){

  std::string returned_string;

  if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z))return "bad_slice";

  if(!isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z))return "bad_slice";

  if(!(all_contained(islc)))return "bad_slice";

  double bar_flash_z = bar_flash(sr,islc); 
	double bar_falsh_x = bar_flash_x(sr,islc);
  double delta=fabs(bar_flash(sr,islc)-islc.charge_center.z);

  if(bar_falsh_x*islc.vertex.x <= 0)return "bad_slice";
  if(bar_flash_z <= -10000)return "bad_slice";
  if(delta <= 0 || delta >= cut_baryc) return "bad_slice";

  int ipfp_mu = -1;
  ipfp_mu = find_muon_newPID(islc,dist_cut,mu1_pfp,mu2_pfp,stitched_len);

  //if(ipfp_mu!=-1)cout << "muon found" << endl;

  if(ipfp_mu!=-1)returned_string = "1";
  else{returned_string = "0";}

  int num_protons =0;
  int num_pions =0;
  int num_showers =0; 

  for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
  {
    if(int(ipfp)==ipfp_mu || int(ipfp)==mu1_pfp || int(ipfp)==mu2_pfp)continue;
    if(id_pfp_newPID(islc, ipfp,dist_cut)==1){num_protons+=1;}
    if(id_pfp_newPID(islc, ipfp,dist_cut)==2){num_pions+=1;}
    if(id_pfp_newPID(islc, ipfp,dist_cut)==3){num_showers+=1;}
  }

  //cout << num_protons << " protons found " << num_pions << " pions found " << num_showers << " showers found " << endl; 

  returned_string = returned_string + Form("%d%d%d",num_protons,num_pions,num_showers);

  return returned_string;
}


bool is1muNp(std::string reco_class)
{
  if(reco_class.compare(0,9,"bad_slice")==0)return false;

  std::vector<int> reco;
  for(char c : reco_class){reco.push_back(c - '0');}

  if(reco[0]==1 && reco[1]>0 && reco[2]==0 and reco[3]==0)return true;

  return false;
}

bool is1muNp_true(std::string true_class)
{
  if(
    true_class == "bad_slice" || 
    true_class == "not_in_FV" ||
    true_class == "not_in_Active" ||
    true_class == "not_contained" ||
    true_class == "cosmic" ||
    true_class == "reco_vertex_isNAN" ||
    true_class == "true_vertex_isNAN" ||
    true_class == "not_numu"
  )return false;

  std::vector<int> true_c;
  for(char c : true_class){true_c.push_back(c - '0');}

  if(true_c[0]==1 && true_c[1]>0 && true_c[2]==0 and true_c[3]==0 and true_c[4]==0)return true;

  return false;
}

bool is1mu_true(std::string true_class)
{
  if(
    true_class == "bad_slice" || 
    true_class == "not_in_FV" ||
    true_class == "not_in_Active" ||
    true_class == "not_contained" ||
    true_class == "cosmic" ||
    true_class == "reco_vertex_isNAN" ||
    true_class == "true_vertex_isNAN" ||
    true_class == "not_numu"
  )return false;

  std::vector<int> true_c;
  for(char c : true_class){true_c.push_back(c - '0');}

  if(true_c[0]==1 && true_c[1]==0 && true_c[2]==0 and true_c[3]==0 and true_c[4]==0)return true;

  return false;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////                                                 Data Struct                                                  //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//std::vector<RecoSlice> slices; 
//int tot_1muNp=0;

//std::string txt_name = "/exp/icarus/data/users/nsommagg/data_struct_multiplane/dumpNuE_provaMaria.txt";
//ofstream dumpNuE(txt_name.c_str());

TTree * tree;
RecoSlice thisslice;

const SpillMultiVar DataLoader([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  std::vector<double> vector_active;
  //if(sr->hdr.run==9435){cout << "trovato" << endl;}
  if(sr->hdr.run>0){ 

  //count how many true 1muNp neutrino interactions
  //int nu_multiplicity=0;
  //for ( auto const& nu : sr->mc.nu )
  //{
    //bool is_1muNp = classification_type_MC(sr,nu);
    //if(is_1muNp)
    //{
      //tot_1muNp++;
      //dumpNuE << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " " << nu.E << " " << nu.index << " " << nu_multiplicity << endl;
      //nu_multiplicity = nu_multiplicity + 1;
      //if(nu.index != 0) cout << "FOUND" << endl; 
    //}
  //}

  int slice_number=-1;
  for (auto const& islc : sr->slc)
  {     
    RecoSlice islice;
    //asking fot non NAN reco and true vertex
    //if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z))continue;
    //if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z))continue;
    
    slice_number++;
    //saving run, subrun, evt
    islice.slice_number = slice_number;
    islice.run = sr->hdr.run;
    islice.subrun = sr->hdr.subrun;
    islice.evt = sr->hdr.evt;

    //slice true classification
    islice.true_slice_classifications = classification_type_generic(sr,islc);

    //slice reco class
    //islice.reco_slice_classification = selection_1muNp_newPID(sr,islc,10,100,-1,-1,-1);

    //islice.is_reco_1muNp = is1muNp(islice.reco_slice_classification);
    islice.is_true_1muNp = is1muNp_true(islice.true_slice_classifications);
    islice.is_true_1mu0p0pi = is1muNp_true(islice.true_slice_classifications);

    //slice interaction type
    islice.is_charged_current = islc.truth.iscc;
    islice.is_neutral_current = islc.truth.isnc;
    islice.genie_mode = islc.truth.genie_mode;
    islice.neutrino_pdg = abs(islc.truth.pdg);
    islice.is_in_FV = isInFV(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z);

    //true neutrino energy
    islice.true_neutrino_energy = islc.truth.E;

    //truth matching completeness and purity
    islice.truth_matching_efficiency = islc.tmatch.eff;
    islice.truth_matching_purity = islc.tmatch.pur;

    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 TrueVtx;
    TrueVtx.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);

    //saving vertex informations
    islice.reco_vertex = {islc.vertex.x, islc.vertex.y, islc.vertex.z};
    islice.true_vertex = {islc.truth.position.x, islc.truth.position.y, islc.truth.position.z};

    //distance between true and reco vertex
    islice.reco_true_vertex_distance = (RecoVtx-TrueVtx).Mag();

    //is in FV ? 
    islice.isinFV = isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z);

    //cahrge-light barycenter
    islice.barycenter=bar_flash(sr,islc);  //This is z!!!
	  islice.barycenter_x=bar_flash_x(sr,islc);
	  islice.bar_charge=islc.charge_center.z;
    islice.delta=fabs(bar_flash(sr,islc)-islc.charge_center.z);

	  //double delta=fabs(barycenter-bar_charge);
    //  if(barycenter_x*islc.vertex.x<=0) return false;
    //  if(barycenter<=-10000)return false;
    //  if(delta <=0 || delta>=100)return false;

    //all contained ? 
    islice.all_tracks_contained = all_contained(islc);

    //is clear cosmic
    islice.is_clear_cosmic = islc.is_clear_cosmic;

    std::vector<track> itracks;
    //cout << islc.reco.npfp << endl;

    int ipfp_mu = find_muon_newPID(islc,10,-1,-1,-1);

    for(std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp)
    {
      track itrack;
      itrack.reco.ipfp = ipfp;
      itrack.truth.ipfp = ipfp;
      itrack.reco.trackscore = islc.reco.pfp[ipfp].trackScore;
      itrack.reco.length = islc.reco.pfp[ipfp].trk.len;
      itrack.reco.is_start_contained = isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,5.0);
      itrack.reco.is_end_contained = isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0);
      itrack.reco.all_in_1_tpc = ((islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0);
      itrack.reco.is_primary = islc.reco.pfp[ipfp].parent_is_primary;

      TVector3 RecoStart;
      TVector3 RecoEnd;
      RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
      RecoEnd.SetXYZ(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);
      itrack.reco.start_distance_from_reco_vertex = (RecoStart-RecoVtx).Mag();
      itrack.reco.end_distance_from_reco_vertex = (RecoEnd-RecoVtx).Mag();

      itrack.reco.start = {islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z};
      itrack.reco.end = {islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z}; 

      TVector3 Start_mom_v_pion;
      Start_mom_v_pion.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x, (islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y, (islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);
      itrack.reco.energy_deposited_as_pi = sqrt(pow(139.570,2)+pow(Start_mom_v_pion.Mag()*1000,2))-139.570;
      TVector3 Start_mom_v_proton;
      Start_mom_v_proton.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
      itrack.reco.energy_deposited_as_pro = sqrt(pow(938.3,2)+pow(Start_mom_v_proton.Mag()*1000,2))-938.3;
      TVector3 Start_mom_v_muon;
      Start_mom_v_muon.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.z);
      itrack.reco.energy_deposited_as_mu = sqrt(pow(105.658,2)+pow(Start_mom_v_muon.Mag()*1000,2))-105.658;

      itrack.reco.shower_energy = islc.reco.pfp[ipfp].shw.plane[2].energy;

      std::vector<double> output;
      std::vector<double> dedx;
      std::vector<double> rr;
      for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
      {
        //if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr>30)continue;
        dedx.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
        rr.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
      } // calo points
      output = chi2_ALG(dedx,rr,0.0,25.0);
      itrack.reco.chi2_as_mu = output[0];
      itrack.reco.chi2_as_pro = output[1];

      itrack.reco.dedx = dedx;
      itrack.reco.rr = rr;
      //itrack.reco.nhits_coll = islc.reco.pfp[ipfp].trk.calo[2].points.size();
      
      std::vector<double> dedx_ind1;
      std::vector<double> rr_ind1;
      std::vector<double> pitch_ind1;
      for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[0].points.size(); ++ihit )
      {
        //if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr>30)continue;
        dedx_ind1.push_back(islc.reco.pfp[ipfp].trk.calo[0].points[ihit].dedx);
        rr_ind1.push_back(islc.reco.pfp[ipfp].trk.calo[0].points[ihit].rr);
        pitch_ind1.push_back(islc.reco.pfp[ipfp].trk.calo[0].points[ihit].pitch);
      } // calo points

      itrack.reco.dedx_ind1 = dedx_ind1;
      itrack.reco.rr_ind1 = rr_ind1;

      std::vector<double> dedx_ind2;
      std::vector<double> rr_ind2;
      std::vector<double> pitch_ind2;
      for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[1].points.size(); ++ihit )
      {
        //if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr>30)continue;
        dedx_ind2.push_back(islc.reco.pfp[ipfp].trk.calo[1].points[ihit].dedx);
        rr_ind2.push_back(islc.reco.pfp[ipfp].trk.calo[1].points[ihit].rr);
        pitch_ind2.push_back(islc.reco.pfp[ipfp].trk.calo[1].points[ihit].pitch);
      } // calo points

      itrack.reco.dedx_ind2 = dedx_ind2;
      itrack.reco.rr_ind2 = rr_ind2;



      itrack.reco.likelihood_ratios = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);


      int bestplane = find_best_plane(islc,ipfp);
      itrack.reco.bestplane = bestplane;

      if(bestplane!=-1)
      {
        itrack.reco.deposited_energy = compute_depE(islc,ipfp,bestplane);

        std::vector<double> temp_daughter_vars = compute_daughter_vars(islc,ipfp);
        itrack.reco.daughter_depE = compute_daughter_vars(islc,ipfp)[0];
        itrack.reco.daughter_angle_end = compute_daughter_vars(islc,ipfp)[1];
        itrack.reco.ek_bestplane = islc.reco.pfp[ipfp].trk.calo[bestplane].ke;
      }
      
      itrack.truth.pdg = islc.reco.pfp[ipfp].trk.truth.p.pdg;
      itrack.truth.true_class = true_selection(sr,islc,ipfp);

      //itrack.reco.predicted_class = PIDclass(islc,ipfp);
      //itrack.reco.predicted_probabilities = PIDproba(islc,ipfp);


      if((int)ipfp == ipfp_mu)itrack.reco.id_ipfp_reco = 0;
      else
      {
        if(id_pfp_newPID(islc,ipfp,10)==1)itrack.reco.id_ipfp_reco = 1;
        if(id_pfp_newPID(islc,ipfp,10)==2)itrack.reco.id_ipfp_reco = 2;
        if(id_pfp_newPID(islc,ipfp,10)==3)itrack.reco.id_ipfp_reco = 3;
      }

      itrack.reco.dir = {islc.reco.pfp[ipfp].trk.dir.x, islc.reco.pfp[ipfp].trk.dir.y, islc.reco.pfp[ipfp].trk.dir.z};
      
      itrack.reco.ek_coll = islc.reco.pfp[ipfp].trk.calo[2].ke;
      itrack.reco.ek_ind1 = islc.reco.pfp[ipfp].trk.calo[0].ke;
      itrack.reco.ek_ind2 = islc.reco.pfp[ipfp].trk.calo[1].ke;

      

      itracks.push_back(itrack);
    }
    
    islice.tracks = itracks;
    //slices.push_back(islice);

    thisslice = islice;
    tree->Fill();
    thisslice = RecoSlice();

  }//loop over all slices
} //only run 9435

  return vector_active;
});
