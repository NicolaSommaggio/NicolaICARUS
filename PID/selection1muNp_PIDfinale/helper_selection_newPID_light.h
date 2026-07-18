#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h" 

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

#include "sbnanaobj/StandardRecord/SRVector3D.h"
#include "TrackMomentumCalculator.h"
#include "helper_newPID.h"
#include "slice_struct.h"

double PROTON_KINETIC_E = 50.;
double CONTAINMENT_CUT = 10.;
double DEP_E_P_RISING_CUT = 0.;

using namespace ana;

//CNAF
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

const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* spill){
  return true;
});

constexpr double PION_MASS            = 139.570; // MeV
constexpr double PROTON_MASS          = 938.3;   // MeV
constexpr double MUON_MASS            = 105.658; // MeV

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

std::vector<double> compute_chi2(const caf::Proxy<caf::SRSlice>& islc,
                                  std::size_t ipfp,
                                  int plane)
{
    std::vector<double> dedx, rr;

    const auto& calo = islc.reco.pfp[ipfp].trk.calo[plane].points;
    if (calo.empty()) return {};

    dedx.reserve(calo.size());
    rr.reserve(calo.size());

    for (const auto& pt : calo) {
        dedx.push_back(pt.dedx);
        rr.push_back(pt.rr);
    }

    return chi2_ALG(dedx, rr, 0, 25);
}

inline double kinetic_energy(double mass, double momentum_GeV)
{
    double p = momentum_GeV * 1000.0; // → MeV
    return std::sqrt(mass*mass + p*p) - mass;
}

int id_pfp_chi2(const caf::Proxy<caf::SRSlice>& islc,
             int ipfp,
             double dist_cut = 10.)
{
    // =========================================================
    // Particle ID strategy:
    //   → only primary PFPs
    //   → vertex associated
    //   → track-like → p/π via χ² + KE
    //   → shower-like → EM shower
    // =========================================================

    const auto& pfp = islc.reco.pfp[ipfp];

    // =========================================================
    // CUT 0 — must be primary and valid
    // =========================================================
    if (!pfp.parent_is_primary){ return 9;}
    if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x) ||
        std::isnan(pfp.trk.len)){ return 9;}

    const TVector3 reco_vtx(islc.vertex.x,
                            islc.vertex.y,
                            islc.vertex.z);

    const TVector3 start(pfp.trk.start.x,
                         pfp.trk.start.y,
                         pfp.trk.start.z);

    const TVector3 end(pfp.trk.end.x,
                       pfp.trk.end.y,
                       pfp.trk.end.z);

    const double min_dist =
        std::min((start - reco_vtx).Mag(), (end - reco_vtx).Mag());

    // =========================================================
    // CUT 1 — vertex association
    // =========================================================
    if (min_dist > 50){return 9;}

    // =========================================================
    // CUT 2 — calorimetry availability
    // =========================================================
    /*
    if (pfp.trk.calo[DEFAULT_PLANE].nhit <= 2 ||
        pfp.trk.calo[DEFAULT_PLANE].points.empty()) {
        {MyFile << "No calo check - Unknown" << endl; return PFPID::Unknown;}
    }*/

    std::vector<double> chi2 = compute_chi2(islc, ipfp, 2);
    if (chi2.size() < 2) {return 9;}

    const double chi2_proton = chi2[1];

    // =========================================================
    // CUT 3 — good track-like objects
    // =========================================================
    if (pfp.trackScore >= 0.5) {

        // ---------- pion hypothesis ----------
        if (chi2_proton >= 100.0) {

            TVector3 p_pi(pfp.trk.dir.x * pfp.trk.rangeP.p_pion,
              pfp.trk.dir.y * pfp.trk.rangeP.p_pion,
              pfp.trk.dir.z * pfp.trk.rangeP.p_pion);

            double ke = kinetic_energy(PION_MASS, p_pi.Mag());

            if ((reco_vtx - start).Mag() < dist_cut && ke >= 25.)
                {return 2;}
        }

        // ---------- proton hypothesis ----------
        if (chi2_proton < 100.0) {

            TVector3 p_pr(pfp.trk.dir.x * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.y * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.z * pfp.trk.rangeP.p_proton);

            double ke = kinetic_energy(PROTON_MASS, p_pr.Mag());

            if ((reco_vtx - start).Mag() < dist_cut && ke >= 50.)
                {return 1;}
        }
    }

    // =========================================================
    // CUT 4 — marginal tracks (proton recovery)
    // =========================================================
    if (pfp.trackScore >= 0.4 &&
        pfp.trackScore <  0.5 &&
        chi2_proton < 100.0) {

        TVector3 p_pr(pfp.trk.dir.x * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.y * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.z * pfp.trk.rangeP.p_proton);

        double ke = kinetic_energy(PROTON_MASS, p_pr.Mag());

        if ((reco_vtx - start).Mag() < dist_cut && ke >= 50.)
            {return 1;}
    }

    // =========================================================
    // CUT 5 — shower-like objects
    // =========================================================
    if (!(pfp.trackScore >= 0.4 &&
        pfp.trackScore <  0.5 &&
        chi2_proton < 100.0) && pfp.trackScore < 0.5){
    //if (pfp.trackScore < PRIMARY_TRACK_SCORE) {

        //if (pfp.shw.plane[DEFAULT_PLANE].nHits <= 2){MyFile << "Shw no calo - Unknown" << endl; return PFPID::Unknown;}
        if (std::isnan(pfp.shw.plane[2].energy)){return 9;}

        double E = pfp.shw.plane[2].energy * 1000.0; // MeV

        if (E > 25.)
            {return 3;}
    }

    return 9;
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

bool all_contained ( const caf::Proxy<caf::SRSlice>& islc ) 
{ 

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
    //if (!(islc.reco.pfp[ipfp].parent_is_primary ))continue; //skip secondaries
    //if(islc.reco.pfp[ipfp].trackScore<0.4)continue; //Want to check only tracks??
    if((islc.reco.pfp[ipfp].trk.start.x*islc.vertex.x)<0){return false;} //not contained if they cross cryostats
    //if not contained return false
    if(!isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,CONTAINMENT_CUT)){return false;}
      
    }   
    
    return true;
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
                if (depE > PROTON_KINETIC_E) {
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
                  if(abs(ipart.pdg)== 2212 && dep_E>PROTON_KINETIC_E){num_protons_above50+=1;} //protons
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

bool is_good_truth_matching(const caf::Proxy<caf::SRSlice>& islc)
{
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);

    if((islc.tmatch.eff > 0.05 && (!islc.is_clear_cosmic) && (vertex_true-vertex_reco).Mag()<100.)) return true;

    return false;
}

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
        isInContained(islc.reco.pfp[current_pfp].trk.end.x,islc.reco.pfp[current_pfp].trk.end.y,islc.reco.pfp[current_pfp].trk.end.z,CONTAINMENT_CUT) && 
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

        if( ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(938.3,2)+pow(Start_mom_v.Mag()*1000,2))-938.3)>=PROTON_KINETIC_E )return 1;
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

  if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z))return "bad_slice_VERTEX_NAN";

  if(!isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z))return "bad_slice_FV";

  if(!(all_contained(islc)))return "bad_slice_CONTAINMENT";

  // test without light -->
  //  double bar_flash_z = bar_flash(sr,islc); 
  //	double bar_falsh_x = bar_flash_x(sr,islc);
  //  double delta=fabs(bar_flash(sr,islc)-islc.charge_center.z);

  //  if(bar_falsh_x*islc.vertex.x <= 0)return "bad_slice_CRT_PMT_and_SAME_CRYO";
  //  if(bar_flash_z <= -10000)return "bad_slice_BARYCENTER";
  //  if(delta <= 0 || delta >= cut_baryc) return "bad_slice_BARYCENTER";
  // <--

  int ipfp_mu = -1;
  ipfp_mu = find_muon_newPID(islc,dist_cut,mu1_pfp,mu2_pfp,stitched_len);

  if(ipfp_mu!=-1)returned_string = "1";
  else{returned_string = "0";}

  int num_protons =0;
  int num_pions =0;
  int num_showers =0; 

  for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
  {

    if(int(ipfp)==ipfp_mu || int(ipfp)==mu1_pfp || int(ipfp)==mu2_pfp){; continue; }

    int pfp_id = id_pfp_newPID(islc, ipfp,dist_cut);

    if(pfp_id==1)
    {
      //modification to suppress cosmic background -->
      int bestplane_temp = find_best_plane(islc,ipfp);
      double depE_temp = compute_depE(islc,ipfp,bestplane_temp);
      int class_temp = PIDclass(islc,ipfp);  
      if(class_temp == 2)
      { 
        if(depE_temp > DEP_E_P_RISING_CUT){num_protons+=1;}
      }
      else{num_protons+=1;}
      //<--
    }
    if(pfp_id==2){num_pions+=1; }
    if(pfp_id==3){num_showers+=1; }

  }

  returned_string = returned_string + Form("%d%d%d",num_protons,num_pions,num_showers);

  return returned_string;
}

std::vector<int> get_protons_pfp ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu)
{
  std::vector<int> v_ipfp_pro;

  //int ipfp_mu = -1;
  //ipfp_mu = find_muon_newPID(islc,10,-1,-1,-1);

  for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
  {
    if(int(ipfp)==ipfp_mu)continue;
    if(id_pfp_newPID(islc,ipfp,10)==1)
    { 
        int bestplane_temp = find_best_plane(islc,ipfp);
        double depE_temp = compute_depE(islc,ipfp,bestplane_temp);
        int class_temp = PIDclass(islc,ipfp);  
        if(class_temp == 2)
        { 
            if(depE_temp > DEP_E_P_RISING_CUT){v_ipfp_pro.push_back(ipfp);}
        }
        else{v_ipfp_pro.push_back(ipfp);}
    }
  }

  return v_ipfp_pro;
}

int get_n_valid_hit(const caf::Proxy<caf::SRSlice>& islc, int ipfp)
{
  int nhit = 0;

  int bestplane = find_best_plane(islc,ipfp);

  if(bestplane == -1)return 0;

  for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
  {
    if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
    {
      nhit ++;
    }
  }  

  return nhit;
}

// input: slice, muon index e longest proton index
double Transverse_angle ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) 
{

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    double p_p_tot = sqrt(p_p_x*p_p_x+p_p_y*p_p_y+p_p_z*p_p_z);

    //Total momentum
    p_tot_x=p_p_x+p_mu_x;
    p_tot_y=p_p_y+p_mu_y;
    p_tot_z=p_p_z+p_mu_z;

    //TRANSVERSE PLANE - ANGLE!
    double norm_mu=sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y);
    double norm_pro=sqrt(p_p_x*p_p_x+p_p_y*p_p_y);


    return (p_mu_x*p_p_x+p_mu_y*p_p_y)/(norm_mu*norm_pro);
}

double T3D_angle_mup ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) 
{

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    TVector3 mu_vector(p_mu_x,p_mu_y,p_mu_z);

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    TVector3 pro_vector(p_p_x,p_p_y,p_p_z);


    return cos(mu_vector.Angle(pro_vector));
}

///////////////////////////////////////////////////////////////////////////////////
//////////                           SELECTION                           //////////
///////////////////////////////////////////////////////////////////////////////////



std::vector<std::string> slices_reco_class;
std::vector<std::string> slices_true_class;
int tot_1muNp = 0;

//ofstream dumpNuE("NEW_MC_dumpNuE_CONTAINMENT_10CM.txt");

bool ismc;

/*
struct _pfp
{
  std::vector<double> _lr;
  double _depE = -999;
  std::vector<double> _KE;
  double _length = -999;
  double _daughter_vars;
  std::vector<double> _dedx;
  std::vector<double> _rr;
  double _theta_xw = -999;
};

struct _slice
{
  int _run = -999;
  int _evt = -999;
  int _slice_counter = -999;
  std::vector<_pfp> _protons;
  _pfp _mu;
  std::string _reco_class = "N/D";
  std::string _true_class = "N/D";
  double _nuE = -999;
  double _mu_pro_angle = -999;
  double _nu_score = -999;

  double _crlongtrkdiry = -999;
  double _bar_flash = -999;
  double bar_flash_x = -999;
};
*/

std::vector<_slice> _slices;



//ofstream dump_reco_class("DUMP_selection_NO_LIGHT_NO_CRT_depE_CUT.txt");

const SpillMultiVar selection([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  ismc = sr->hdr.ismc;

  if(sr->hdr.ismc)
  {
    //count how many true 1muNp neutrino interactions
    for ( auto const& nu : sr->mc.nu )
    {
      bool is_1muNp = classification_type_MC(sr,nu);
      if(is_1muNp)
      {
        //dumpNuE << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " " << nu.E << " " << nu.index << endl;
        if(nu.index == 0)tot_1muNp++;
      }
    }
  }

    std::vector<double> vector_active;

    int slice_counter = -1;
    for(const auto &islc : sr->slc)
    { 
      slice_counter ++;

      std::string slice_class_true = "N/D";
      //double nuE = -999;
      if (sr->hdr.ismc)
      {
        //nuE = islc.truth.E;
        slice_class_true = classification_type_generic(sr,islc);
      }
      slices_true_class.push_back(slice_class_true);

      std::string slice_class = selection_1muNp_newPID(sr,islc,10,100,-1,-1,-1);
      slices_reco_class.push_back(slice_class);

      //double muon_length = -999;
      //double mu_pro_angle = -999;
      //double length_leading_proton = -999;

      _slice slice;

      slice._run = sr->hdr.run;
      slice._evt = sr->hdr.evt;
      slice._slice_counter = slice_counter;
      slice._true_class = slice_class_true;
      slice._reco_class = slice_class;
      slice._nuE = islc.truth.E;
      slice._nu_score = islc.nu_score;
      slice._bar_flash = bar_flash(sr,islc);
      slice._bar_flash_x = bar_flash_x(sr,islc);
      slice._crlongtrkdiry = islc.nuid.crlongtrkdiry;

      _pfp muone;
      std::vector<_pfp> protoni;

      if(is1muNp(slice_class))
      {
        int ipfp_mu = find_muon_newPID(islc,10,-1,-1,-1);

        //muon_length = islc.reco.pfp[ipfp_mu].trk.len;

        int mu_bestplane = find_best_plane(islc,ipfp_mu);
        muone._length = islc.reco.pfp[ipfp_mu].trk.len;
        muone._daughter_vars = compute_daughter_vars(islc,ipfp_mu);
        muone._KE = compute_ke(islc,ipfp_mu,mu_bestplane);
        muone._lr = likelihood(islc,ipfp_mu, prob_d_coll, prob_d_ind1, prob_d_ind2);
        muone._depE = compute_depE(islc, ipfp_mu, mu_bestplane);
        
        double pitch_temp_mu = 0;
        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[mu_bestplane].points.size(); ++ihit )
        {
          muone._rr.push_back(islc.reco.pfp[ipfp_mu].trk.calo[mu_bestplane].points[ihit].rr);
          muone._dedx.push_back(islc.reco.pfp[ipfp_mu].trk.calo[mu_bestplane].points[ihit].dedx);
          pitch_temp_mu += islc.reco.pfp[ipfp_mu].trk.calo[mu_bestplane].points[ihit].pitch;
        } 
        double avg_pitch_mu = pitch_temp_mu / islc.reco.pfp[ipfp_mu].trk.calo[mu_bestplane].points.size();

        muone._theta_xw = std::atan( islc.reco.pfp[ipfp_mu].trk.dir.x * avg_pitch_mu / 0.3 );

        double maxL = 0;
        int ipfp_maxL = -1;
        std::vector<int> v_ipfp_pro = get_protons_pfp(sr,islc,ipfp_mu);
        for(const auto &ipfp_pro : v_ipfp_pro)
        {
          _pfp protone;

          int pro_bestplane = find_best_plane(islc,ipfp_pro);
          protone._length = islc.reco.pfp[ipfp_pro].trk.len;
          protone._daughter_vars = compute_daughter_vars(islc,ipfp_pro);
          protone._KE = compute_ke(islc,ipfp_pro,pro_bestplane);
          protone._lr = likelihood(islc,ipfp_pro, prob_d_coll, prob_d_ind1, prob_d_ind2);
          protone._depE = compute_depE(islc, ipfp_pro, pro_bestplane);
        
          double pitch_temp_pro = 0;
          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_pro].trk.calo[pro_bestplane].points.size(); ++ihit )
          {
            protone._rr.push_back(islc.reco.pfp[ipfp_pro].trk.calo[pro_bestplane].points[ihit].rr);
            protone._dedx.push_back(islc.reco.pfp[ipfp_pro].trk.calo[pro_bestplane].points[ihit].dedx);
            pitch_temp_pro += islc.reco.pfp[ipfp_pro].trk.calo[pro_bestplane].points[ihit].pitch;
          } 
          double avg_pitch_pro = pitch_temp_pro / islc.reco.pfp[ipfp_pro].trk.calo[pro_bestplane].points.size();

          protone._theta_xw = std::atan( islc.reco.pfp[ipfp_pro].trk.dir.x * avg_pitch_pro / 0.3 );

          protoni.push_back(protone);

          if(islc.reco.pfp[ipfp_pro].trk.len > maxL)
          {
            maxL = islc.reco.pfp[ipfp_pro].trk.len;
            ipfp_maxL = ipfp_pro;
          }
        }


        //mu_pro_angle = T3D_angle_mup(islc,ipfp_mu,ipfp_maxL);
        slice._mu_pro_angle = T3D_angle_mup(islc,ipfp_mu,ipfp_maxL);
        //length_leading_proton = maxL;

        slice._mu = muone;
        slice._protons = protoni;

        _slices.push_back(slice);
      }

      //if(sr->hdr.ismc)
      //{
        //dump_reco_class << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " " << nuE << " " << slice_class << " " << slice_class_true << " " << mu_pro_angle << " " << length_leading_proton << " " << muon_length << " " << islc.nu_score << endl;
      //}
      //else
      //{
        //if(is1muNp(slice_class))
        //{
          //dump_reco_class << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " " << nuE << " " << slice_class << " " << slice_class_true << " " << mu_pro_angle << " " << length_leading_proton << " " << muon_length << " " << islc.nu_score << endl;
        //}
      //}

    }//loop over all slices

  return vector_active;
});




/*

const SpillMultiVar files_check([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    return vector_active;
});
*/

/*
std::string file_sample = "DATA";

ofstream muon_candidates(Form("PIDSYST_CHI2_MUON_CANDIDATES_%s.txt",file_sample.c_str()));
ofstream selected_muons(Form("PIDSYST_CHI2_SELECTED_MUON_%s.txt",file_sample.c_str()));
ofstream selected_protons(Form("PIDSYST_CHI2_SELECTED_PROTONS_%s.txt",file_sample.c_str()));
ofstream selected_pions(Form("PIDSYST_CHI2_SELECTED_PIONS_%s.txt",file_sample.c_str()));

const SpillMultiVar dump_BDT_vars([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    int slice_counter = -1;
    for(const auto &islc : sr->slc)
    { 
      slice_counter ++;

      if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z))continue;

      if(!isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z))continue;

      if(!(all_contained(islc)))continue;

	    double bar_falsh_x = bar_flash_x(sr,islc);

      if(bar_falsh_x*islc.vertex.x <= 0)continue;

      // MUON SEARCH

      int ipfp_mu = -1;

      double max_p_as_mu = -1;

      TVector3 RecoVtx;
      RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
      TVector3 RecoStart;

      double max_track_length = -1;

      for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
      {

        if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
        
        RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
        
        bool is_primary = islc.reco.pfp[ipfp].parent_is_primary;

        if(islc.reco.pfp[ipfp].trackScore<0.5)continue;
        if((RecoVtx-RecoStart).Mag() >= 10)continue;
        if(islc.reco.pfp[ipfp].trk.len <= 50 )continue;
        if(!(isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,CONTAINMENT_CUT)))continue;
        if((islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x) <= 0)continue;
        if(!is_primary)continue;

        // --> at this point it is a muon candidate!

        // --> retrieve BDT variables

        int bestplane = 2;

        // --> if(bestplane == -1)continue; //no valid hits (RR < 25 cm and RR > 1 cm and dEdx < 30 MeV/cm and dEdx > 1 MeV/cm) in neither of the 3 planes

        std::vector<double> temp_dedx;
        std::vector<double> temp_rr;
        std::vector<double> temp_pitch;
        std::vector<double> temp_x;
        std::vector<double> temp_mult;
        // --> std::vector<double> temp_proba = PIDproba(islc,ipfp);

        double average_pitch = 0;
        double track_dir_x;

        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit )
        {
          temp_rr.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr);
          temp_dedx.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx);
          temp_pitch.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch);
          temp_x.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].x);
          temp_mult.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].mult);

          average_pitch = average_pitch + islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch;
        }

        average_pitch = average_pitch / islc.reco.pfp[ipfp].trk.calo[bestplane].points.size();
        track_dir_x = islc.reco.pfp[ipfp].trk.dir.x;

        double thetaXW = std::atan(track_dir_x * average_pitch/0.3);

        std::vector<double> temp_ke = compute_ke(islc,ipfp,bestplane);

        // --> std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
        double depE = compute_depE(islc,ipfp,bestplane);
        // --> std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

        muon_candidates << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
        //for(const auto &l : lr){muon_candidates << l << " ";}
        //muon_candidates << depE << " " << dvars[0] << " " << dvars[1] << " ";
        muon_candidates << depE << " ";
        //muon_candidates << bestplane << " ";
        //for(const auto &p : temp_proba){muon_candidates << p << " ";}
        for(const auto &dedx : temp_dedx){muon_candidates << dedx << "_";}
        muon_candidates << " ";
        for(const auto &rr : temp_rr){muon_candidates << rr << "_";}
        muon_candidates << " ";
        for(const auto &pitch : temp_pitch){muon_candidates << pitch << "_";}
        muon_candidates << " ";
        for(const auto &x : temp_x){muon_candidates << x << "_";}
        muon_candidates << " ";
        for(const auto &m : temp_mult){muon_candidates << m << "_";}
        muon_candidates << " ";
        for(const auto &ke : temp_ke){muon_candidates << ke << " ";}
        muon_candidates << thetaXW << " " << track_dir_x;
        muon_candidates << endl;
        


        // --> query the PID

        // --> int pid_pred = PIDclass(islc,ipfp);
        // --> if(pid_pred == 2 || pid_pred == 3)continue;

        // --> std::vector<double> prediction_proba;
        // --> prediction_proba = PIDproba(islc,ipfp); //probabilities for the track to be each of the 6 classes

        //double imu_prob = std::max(prediction_proba[0], prediction_proba[1]);
 
        //if(imu_prob > max_p_as_mu)
        //{
          //max_p_as_mu=imu_prob;
          //ipfp_mu=ipfp;
        //}

        std::vector<double> output_chi2 = compute_chi2(islc,ipfp,2);
        //std::vector<double> output_chi2 = chi2_ALG(temp_dedx,temp_rr,0,25);
        if(output_chi2.size() > 1 && islc.reco.pfp[ipfp].trk.len > max_track_length && output_chi2[0] < 30. && output_chi2[1] > 60.)
        {
          ipfp_mu = ipfp;
          max_track_length = islc.reco.pfp[ipfp].trk.len;
        }



      }//loop of pfp to find muon

      if(ipfp_mu != -1)
      {
        // --> int bestplane = find_best_plane(islc,ipfp_mu);

        int bestplane = 2;

        std::vector<double> temp_dedx;
        std::vector<double> temp_rr;
        std::vector<double> temp_pitch;
        std::vector<double> temp_x;
        std::vector<double> temp_mult;
        // --> std::vector<double> temp_proba = PIDproba(islc,ipfp_mu);

        double average_pitch = 0;
        double track_dir_x;

        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points.size(); ++ihit )
        {
          temp_rr.push_back(islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points[ihit].rr);
          temp_dedx.push_back(islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points[ihit].dedx);
          temp_pitch.push_back(islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points[ihit].pitch);
          temp_x.push_back(islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points[ihit].x);
          temp_mult.push_back(islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points[ihit].mult);

          average_pitch = average_pitch + islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points[ihit].pitch;
        }

        average_pitch = average_pitch / islc.reco.pfp[ipfp_mu].trk.calo[bestplane].points.size();
        track_dir_x = islc.reco.pfp[ipfp_mu].trk.dir.x;

        double thetaXW = std::atan(track_dir_x * average_pitch/0.3);

        std::vector<double> temp_ke = compute_ke(islc,ipfp_mu,bestplane);

        // --> std::vector<double> lr = likelihood(islc,ipfp_mu, prob_d_coll, prob_d_ind1, prob_d_ind2);
        double depE = compute_depE(islc,ipfp_mu,bestplane);
        // --> std::vector<double> dvars = compute_daughter_vars(islc,ipfp_mu);

        selected_muons << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
        // --> for(const auto &l : lr){selected_muons << l << " ";}
        // --> selected_muons << depE << " " << dvars[0] << " " << dvars[1] << " ";

        selected_muons << depE << " ";

        // --> selected_muons << bestplane << " ";
        // --> for(const auto &p : temp_proba){selected_muons << p << " ";}
        for(const auto &dedx : temp_dedx){selected_muons << dedx << "_";}
        selected_muons << " ";
        for(const auto &rr : temp_rr){selected_muons << rr << "_";}
        selected_muons << " ";
        for(const auto &pitch : temp_pitch){selected_muons << pitch << "_";}
        selected_muons << " ";
        for(const auto &x : temp_x){selected_muons << x << "_";}
        selected_muons << " ";
        for(const auto &m : temp_mult){selected_muons << m << "_";}
        selected_muons << " ";
        for(const auto &ke : temp_ke){selected_muons << ke << " ";}
        selected_muons << thetaXW << " " << track_dir_x;
        selected_muons << endl;
      }

      // PION AND PROTON FOUND

      for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
      {
        if((int)ipfp == ipfp_mu)continue;

        // --> int id_pfp = id_pfp_newPID(islc, ipfp, 10);

        int id_pfp = id_pfp_chi2(islc,ipfp);


        if(id_pfp == 1)
        {
          // --> int bestplane = find_best_plane(islc,ipfp);

          int bestplane = 2;

          std::vector<double> temp_dedx;
          std::vector<double> temp_rr;
          std::vector<double> temp_pitch;
          std::vector<double> temp_x;
          std::vector<double> temp_mult;
          // --> std::vector<double> temp_proba = PIDproba(islc,ipfp);

          double average_pitch = 0;
          double track_dir_x;

          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit )
          {
            temp_rr.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr);
            temp_dedx.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx);
            temp_pitch.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch);
            temp_x.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].x);
            temp_mult.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].mult);

            average_pitch = average_pitch + islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch;
          }

          average_pitch = average_pitch / islc.reco.pfp[ipfp].trk.calo[bestplane].points.size();
          track_dir_x = islc.reco.pfp[ipfp].trk.dir.x;

          double thetaXW = std::atan(track_dir_x * average_pitch/0.3);

          std::vector<double> temp_ke = compute_ke(islc,ipfp,bestplane);

          // --> std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
          double depE = compute_depE(islc,ipfp,bestplane);
          // --> std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

          selected_protons << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
          // --> for(const auto &l : lr){selected_protons << l << " ";}
          // --> selected_protons << depE << " " << dvars[0] << " " << dvars[1] << " ";
          selected_protons << depE << " ";
          // --> selected_protons << bestplane << " ";
          // --> for(const auto &p : temp_proba){selected_protons << p << " ";}
          for(const auto &dedx : temp_dedx){selected_protons << dedx << "_";}
          selected_protons << " ";
          for(const auto &rr : temp_rr){selected_protons << rr << "_";}
          selected_protons << " ";
          for(const auto &pitch : temp_pitch){selected_protons << pitch << "_";}
          selected_protons << " ";
          for(const auto &x : temp_x){selected_protons << x << "_";}
          selected_protons << " ";
          for(const auto &m : temp_mult){selected_protons << m << "_";}
          selected_protons << " ";
          for(const auto &ke : temp_ke){selected_protons << ke << " ";}
          selected_protons << thetaXW << " " << track_dir_x;
          selected_protons << endl;
        }
        if(id_pfp == 2)
        {
          // --> int bestplane = find_best_plane(islc,ipfp);

          int bestplane = 2;

          std::vector<double> temp_dedx;
          std::vector<double> temp_rr;
          std::vector<double> temp_pitch;
          std::vector<double> temp_x;
          std::vector<double> temp_mult;
          // --> std::vector<double> temp_proba = PIDproba(islc,ipfp);

          double average_pitch = 0;
          double track_dir_x;

          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit )
          {
            temp_rr.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr);
            temp_dedx.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx);
            temp_pitch.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch);
            temp_x.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].x);
            temp_mult.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].mult);

            average_pitch = average_pitch + islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch;
          }

          average_pitch = average_pitch / islc.reco.pfp[ipfp].trk.calo[bestplane].points.size();
          track_dir_x = islc.reco.pfp[ipfp].trk.dir.x;

          double thetaXW = std::atan(track_dir_x * average_pitch/0.3);

          std::vector<double> temp_ke = compute_ke(islc,ipfp,bestplane);

          // --> std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
          double depE = compute_depE(islc,ipfp,bestplane);
          // --> std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

          selected_pions << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
          // --> for(const auto &l : lr){selected_pions << l << " ";}
          // --> selected_pions << depE << " " << dvars[0] << " " << dvars[1] << " ";

          selected_pions << depE << " ";

          // --> selected_pions << bestplane << " ";
          // --> for(const auto &p : temp_proba){selected_pions << p << " ";}
          for(const auto &dedx : temp_dedx){selected_pions << dedx << "_";}
          selected_pions << " ";
          for(const auto &rr : temp_rr){selected_pions << rr << "_";}
          selected_pions << " ";
          for(const auto &pitch : temp_pitch){selected_pions << pitch << "_";}
          selected_pions << " ";
          for(const auto &x : temp_x){selected_pions << x << "_";}
          selected_pions << " ";
          for(const auto &m : temp_mult){selected_pions << m << "_";}
          selected_pions << " ";
          for(const auto &ke : temp_ke){selected_pions << ke << " ";}
          selected_pions << thetaXW << " " << track_dir_x;
          selected_pions << endl;
        }

      }


    }//loop on slices

    return vector_active;
});

*/
