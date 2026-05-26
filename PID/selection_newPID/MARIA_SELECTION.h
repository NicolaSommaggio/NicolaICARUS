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
#include "TProfile.h"
using namespace ana;



//ofstream MyFile_truth("Debug_endy.txt");
//ofstream MyFile_truth2("Selected_clean_truenu.txt");

//ofstream MyFile2("Selected_clean_true_evt.txt");

TFile* file = TFile::Open("/exp/icarus/app/users/marterop/dev_areas/dEdxrestemplates.root");

//TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");


std::vector<double> chi2_ALG(std::vector<double> &dEdx,std::vector<double> &RR, double rr_min, double rr_max)
{
    //The output is chi2s



    double threshold = 0.5;
    double max_rr = rr_max;
    double min_rr = rr_min;

    std::vector<float> trkdedx;
    std::vector<float> trkres;
    std::vector<double> vpida;
    

    for(std::size_t i(0); i<dEdx.size(); ++i){
      if(i==0 || i==dEdx.size()-1)continue;
        if(RR[i]<max_rr && RR[i]>rr_min ){
            if(!std::isnan(dEdx[i])&&!std::isnan(RR[i])){trkdedx.push_back(dEdx[i]);trkres.push_back(RR[i]);}}       
    }


    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double avgdedx = 0;
    double PIDA = 0;

    int used_trkres = 0;
    for (unsigned i = 0; i < trkdedx.size(); ++i) { //hits
      //ignore the first and the last point
      //if (i == 0 || i == trkdedx.size() - 1) continue;
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



bool isInFV (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
  //need to add a check to avoid having the vtx slice in one cryo and the end track in the other cryo

    //fiducial volume for the dangling cable 
  if(x>210.0 && y > 60.0 && z> 290.0 && z< 390.0 ) return false;
  //Cathode padding
  //if(std::abs(x)>208.79 && std::abs(x)<211.79) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
			( x >  61.94 + 25 && x <  358.49 - 25 )) &&
		  ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
		  ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
}

bool isInContained (double x, double y, double z, double dist)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

//remove triangle regions -- corners 
  if( y<(1.732007*z-1687.5114))return false;
  if( y>(-1.732007*z+1640.6114))return false;
  if( y>(1.732007*z+1640.6114))return false;
  if( y<(-1.732007*z-1687.5114))return false;

//5 cm containment for Y in both cryo
  return (( ( x < -61.94 - dist && x > -358.49 + dist ) ||
			( x >  61.94 + dist && x <  358.49 - dist )) &&
		  //( ( y > -181.86 + (dist+1) && y < 134.96 - (dist+1) ) &&
		  ( ( y > -181.86 + (dist) && y < 134.96 - (dist) ) &&

		  ( z > -894.95 + dist && z < 894.95 - dist ) ));

}


bool isInActive (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 && x > -358.49 ) ||
			( x >  61.94 && x <  358.49)) &&
		  ( ( y > -181.86 && y < 134.96)  &&
		    ( z > -894.95 && z < 894.95) ));
}

bool isInContained (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

//remove triangle regions -- corners 
    if( y<(1.732007*z-1687.5114))return false;
    if( y>(-1.732007*z+1640.6114))return false;
    if( y>(1.732007*z+1640.6114))return false;
    if( y<(-1.732007*z-1687.5114))return false;

  return (( ( x < -61.94 - 5 && x > -358.49 + 5 ) ||
			( x >  61.94 + 5 && x <  358.49 - 5 )) &&
		  ( ( y > -181.86 + 5 && y < 134.96 - 5 ) &&
		  ( z > -894.95 + 5 && z < 894.95 - 5 ) ));

}


// ----------------- Tunable cuts -----------------
constexpr int    DEFAULT_PLANE        = 2;
constexpr double PRIMARY_TRACK_SCORE  = 0.5;
constexpr double SECONDARY_TRACK_LOW  = 0.4;
constexpr double VTX_MAX_DIST         = 50.0; // cm

constexpr double MIN_MUON_LENGTH   = 50.0; // cm
constexpr double MAX_CHI2_MUON     = 30.0;
constexpr double MIN_CHI2_PROTON   = 60.0;
constexpr double CONTAINMENT_CUT   = 10.0;  // cm
//constexpr double CONTAINMENT_CUT   = 5.0;  // cm


constexpr double CALO_RR_MIN          = 0.0;
constexpr double CALO_RR_MAX          = 25.0;

constexpr double PION_KE_MIN          = 25.0; // MeV
constexpr double PROTON_KE_MIN        = 50.0; // MeV

constexpr double PION_MASS            = 139.570; // MeV
constexpr double PROTON_MASS          = 938.3;   // MeV
constexpr double MUON_MASS            = 105.658; // MeV
constexpr double PROTON_BINDING_ENERGY = 30.9; // MeV (argon effective)


enum class Chi2PID {
    Muon   = 0,
    Proton = 1,
    Kaon   = 2,
    Pion   = 3
};

struct ParticleCounts {
    int protons  = 0;
    int pions    = 0;
    int showers  = 0;
};

enum class PFPID {
    Unknown = 0,
    Proton  = 1,
    Pion    = 2,
    Shower  = 3, 
};


enum class TruthClass {
    kOneMuOneP,   
    kOneMuNp,     
    kOther,      
    kCosmic,     
    kInvalid     
};


enum class TruthPID {
    Muon,
    Proton,
    ChargedPion,
    NeutralPion,
    Photon,
    Electron,
    Other
};

inline std::string TruthClassToString(TruthClass cls)
{
    switch (cls) {
        case TruthClass::kOneMuOneP: return "1mu1p";
        case TruthClass::kOneMuNp:   return "1muNp";
        case TruthClass::kOther:     return "Other";
        case TruthClass::kCosmic:     return "Cosmic";
        case TruthClass::kInvalid:   return "Invalid";
        default:                     return "Unknown";
    }
}

inline TruthPID classify_truth_pid(int pdg)
{
    switch (std::abs(pdg)) {
        case 13:   return TruthPID::Muon;
        case 2212: return TruthPID::Proton;
        case 211:  return TruthPID::ChargedPion;
        case 111:  return TruthPID::NeutralPion;
        case 22:   return TruthPID::Photon;
        case 11:   return TruthPID::Electron;
        default:   return TruthPID::Other;
    }
}

inline std::string TruthPIDToString(TruthPID pid)
{
    switch (pid) {
        case TruthPID::Muon:         return "Muon";
        case TruthPID::Proton:       return "Proton";
        case TruthPID::ChargedPion:  return "ChargedPion";
        case TruthPID::NeutralPion:  return "NeutralPion";
        case TruthPID::Photon:       return "Photon";
        case TruthPID::Electron:     return "Electron";
        default:                     return "Other";
    }
}



inline bool is_charged_lepton_or_hadron(TruthPID pid)
{
    switch (pid) {
        case TruthPID::Muon:
        case TruthPID::Proton:
        case TruthPID::ChargedPion:
        case TruthPID::Electron:
            return true;
        default:
            return false;
    }
}

inline bool same_g4id(unsigned int parent, int g4id)
{
    if (g4id < 0) return false;
    return parent == static_cast<unsigned int>(g4id);
}


bool all_contained ( const caf::Proxy<caf::SRSlice>& islc ) { 

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
    if(
        std::isnan(islc.reco.pfp[ipfp].trk.start.x) || 
        std::isnan(islc.reco.pfp[ipfp].trk.end.x) || 
        std::isnan(islc.reco.pfp[ipfp].trk.len)
        ) continue;
    //not contained if they cross cryostats
    if((islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)<0){return false;} 
    //if not contained return false
    if(!isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,CONTAINMENT_CUT)){return false;}
      
    }   
    
return true;
}

bool all_contained_truth(const caf::SRSpillProxy* sr,
                         const caf::Proxy<caf::SRSlice>& islc)
{
    for (const auto& ipart : islc.truth.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        TruthPID pid = classify_truth_pid(ipart.pdg);

        // Only visible charged primaries
        if (is_charged_lepton_or_hadron(pid)) {

            if (!isInContained(ipart.end.x, ipart.end.y, ipart.end.z, CONTAINMENT_CUT))
            {
                //MyFile_truth << "Primary is not contained " <<  TruthPIDToString(pid) << endl; 
                return false;}
        }

        // Daughters
        for (const auto& d : sr->true_particles) {

            if (!same_g4id(d.parent, ipart.G4ID) || d.cryostat < 0) continue;

            TruthPID dpid = classify_truth_pid(d.pdg);

            if (is_charged_lepton_or_hadron(dpid)) {

                if (d.end.x == -9999 || d.end.y == -9999 || d.end.z == -9999) continue;

                if (!isInContained(d.end.x, d.end.y, d.end.z, CONTAINMENT_CUT))
                    {
                       // MyFile_truth << "Daughter is not contained " <<  TruthPIDToString(dpid) << endl; 
                    return false;}
            }
        }
    }

    return true;
}

bool all_contained_mc(const caf::SRSpillProxy* sr,
                         const caf::Proxy<caf::SRTrueInteraction>& nu)
{
    for (const auto& ipart : nu.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        TruthPID pid = classify_truth_pid(ipart.pdg);

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

            TruthPID dpid = classify_truth_pid(d.pdg);

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


TruthClass classification_type_debug(const caf::SRSpillProxy* sr,
                                      const caf::Proxy<caf::SRSlice>& islc)
{
   // MyFile_truth << "\n================ NEW EVENT ================\n";
   // MyFile_truth << "RUN EVT " << sr->hdr.run << " " << sr->hdr.evt << endl; 


    // ----------------- Sanity -----------------
    if(islc.truth.index<0){        
     //   MyFile_truth << "FAIL: Cosmic\n";
        return TruthClass::kCosmic;}

    if (std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)) {
     //   MyFile_truth << "FAIL: reco vertex is NaN\n";
        return TruthClass::kInvalid;
    }

    if (std::isnan(islc.truth.position.x) ||
        std::isnan(islc.truth.position.y) ||
        std::isnan(islc.truth.position.z)) {
      //  MyFile_truth << "FAIL: truth vertex is NaN\n";
        return TruthClass::kInvalid;
    }

    // ----------------- Neutrino + volume -----------------
    if (std::abs(islc.truth.pdg) != 14 || !islc.truth.iscc) {
     //   MyFile_truth << "FAIL: not numu CC (pdg=" << islc.truth.pdg
     //       << ", iscc=" << islc.truth.iscc << ")\n";
        return TruthClass::kOther;
    }

    if (!isInActive(islc.truth.position.x,
                    islc.truth.position.y,
                    islc.truth.position.z)) {
      //  MyFile_truth << "FAIL: interaction not in active volume\n";
        return TruthClass::kOther;
    }

    if (!isInFV(islc.truth.position.x,
                islc.truth.position.y,
                islc.truth.position.z)) {
      //  MyFile_truth << "FAIL: interaction not in fiducial volume\n";
        return TruthClass::kOther;
    }

    // ----------------- Counters -----------------
    int    n_muons = 0;
    int    n_protons_above = 0;
    double muon_length = 0.0;

    const int use_plane = DEFAULT_PLANE;

    // ----------------- Primary loop -----------------
    for (const auto& ipart : islc.truth.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        TruthPID pid = classify_truth_pid(ipart.pdg); 
        double depE = ipart.plane[ipart.cryostat][use_plane].visE * 1000.0;
/*
        MyFile_truth << "Primary G4ID=" << ipart.G4ID
            << " PDG=" << ipart.pdg
            << " depE=" << depE << " MeV\n";
*/
        // -------- Muons --------
        if (pid == TruthPID::Muon) {
            ++n_muons;
            muon_length = ipart.length;
           // MyFile_truth << "  -> muon, length = " << muon_length << " cm\n";
        }

        // -------- Charged pions --------
        if (pid == TruthPID::ChargedPion && depE > PION_KE_MIN) {
          //  MyFile_truth << "FAIL: charged pion above " << PION_KE_MIN << " MeV\n";
            return TruthClass::kOther;
        }

        // -------- Neutral pions --------
        if (pid == TruthPID::NeutralPion) {
            for (const auto& d : sr->true_particles) {
                if (same_g4id(d.parent, ipart.G4ID) &&
                    classify_truth_pid(d.pdg) == TruthPID::Photon) {

                    double eg = d.plane[ipart.cryostat][use_plane].visE * 1000.0;
                  //  MyFile_truth << "  -> pi0 daughter gamma, E = " << eg << " MeV\n";

                    if (eg > PION_KE_MIN) {
                     //   MyFile_truth << "FAIL: pi0 gamma above " << PION_KE_MIN << " MeV\n";
                        return TruthClass::kOther;
                    }
                }
            }
        }

        // -------- Photons --------
        if (pid == TruthPID::Photon) {
            for (const auto& d : sr->true_particles){
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;
            }
           // MyFile_truth << "  -> total photon depE = " << depE << " MeV\n";

            if (depE > PION_KE_MIN) {
               // MyFile_truth << "FAIL: photon above " << PION_KE_MIN << " MeV\n";
                return TruthClass::kOther;
            }
        }

        // -------- Protons --------
        if (pid == TruthPID::Proton) {
            for (const auto& d : sr->true_particles){
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;
            }
          //  MyFile_truth << "  -> proton total depE = " << depE << " MeV\n";

            if (depE > PROTON_KE_MIN) {
                ++n_protons_above;
              //  MyFile_truth << "  -> proton above threshold\n";
            }
        }
    }

    // ----------------- Containment -----------------
    if (!all_contained_truth(sr, islc)) {
      //  MyFile_truth << "FAIL: not all contained\n";
        return TruthClass::kOther;
    }

    // ----------------- Final classification -----------------
    /*
    MyFile_truth << "SUMMARY: n_muons=" << n_muons
        << " n_protons_above=" << n_protons_above
        << " muon_length=" << muon_length << " cm\n";
*/
    if (n_muons == 1 && muon_length > MIN_MUON_LENGTH) {

        if (n_protons_above == 1) {
          //  MyFile_truth << "PASS: classified as 1mu1p\n";
            return TruthClass::kOneMuOneP;
        }

        if (n_protons_above > 1) {
           // MyFile_truth << "PASS: classified as 1muNp\n";
            return TruthClass::kOneMuNp;
        }

       // MyFile_truth << "FAIL: no proton above threshold\n";
        return TruthClass::kOther;
    }

   // MyFile_truth << "FAIL: muon multiplicity or muon length cut\n";
    return TruthClass::kOther;
}

TruthClass classification_type_MC(const caf::SRSpillProxy* sr,
                                      const caf::Proxy<caf::SRTrueInteraction>& nu)
{
    // ----------------- Sanity -----------------
    if (std::isnan(nu.position.x) ||
        std::isnan(nu.position.y) ||
        std::isnan(nu.position.z)) {
        return TruthClass::kInvalid;
    }

    // ----------------- Neutrino + volume -----------------
    if (std::abs(nu.pdg) != 14 || !nu.iscc) {
        return TruthClass::kOther;
    }

    if (!isInActive(nu.position.x,
                    nu.position.y,
                    nu.position.z)) {
        return TruthClass::kOther;
    }

    if (!isInFV(nu.position.x,
                nu.position.y,
                nu.position.z)) {
        return TruthClass::kOther;
    }

    // ----------------- Counters -----------------
    int    n_muons = 0;
    int    n_protons_above = 0;
    double muon_length = 0.0;

    const int use_plane = DEFAULT_PLANE;

    // ----------------- Primary loop -----------------
    for (const auto& ipart : nu.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        TruthPID pid = classify_truth_pid(ipart.pdg); 
        double depE = ipart.plane[ipart.cryostat][use_plane].visE * 1000.0;

        // -------- Muons --------
        if (pid == TruthPID::Muon) {
            ++n_muons;
            muon_length = ipart.length;
        }

        // -------- Charged pions --------
        if (pid == TruthPID::ChargedPion && depE > PION_KE_MIN) {
            return TruthClass::kOther;
        }

        // -------- Neutral pions --------
        if (pid == TruthPID::NeutralPion) {
            for (const auto& d : sr->true_particles) {
                if (same_g4id(d.parent, ipart.G4ID) &&
                    classify_truth_pid(d.pdg) == TruthPID::Photon) {

                    double eg = d.plane[ipart.cryostat][use_plane].visE * 1000.0;

                    if (eg > PION_KE_MIN) {
                        return TruthClass::kOther;
                    }
                }
            }
        }

        // -------- Photons --------
        if (pid == TruthPID::Photon) {
            for (const auto& d : sr->true_particles){
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;
            }

                if (depE > PION_KE_MIN) {
                    return TruthClass::kOther;
                }
            
        }

        // -------- Protons --------
        if (pid == TruthPID::Proton) {
            for (const auto& d : sr->true_particles){
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;
            }
                if (depE > PROTON_KE_MIN) {
                    ++n_protons_above;
                }
            
        }
    }

    // ----------------- Containment -----------------
    if (!all_contained_mc(sr, nu)) {
        return TruthClass::kOther;
    }

    // ----------------- Final classification -----------------

    if (n_muons == 1 && muon_length > MIN_MUON_LENGTH) {

        if (n_protons_above == 1) {
                //MyFile_truth2 << "\n================ NEW EVENT ================\n";
                //MyFile_truth2 << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " 1" << endl; 

            return TruthClass::kOneMuOneP;}

        if (n_protons_above > 1) {
                //MyFile_truth2 << "\n================ NEW EVENT ================\n";
                //MyFile_truth2 << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " " << n_protons_above << endl; 

            return TruthClass::kOneMuNp;}

        return TruthClass::kOther;
    }

    return TruthClass::kOther;
}


inline double kinetic_energy(double mass, double momentum_GeV)
{
    double p = momentum_GeV * 1000.0; // → MeV
    return std::sqrt(mass*mass + p*p) - mass;
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
    //if(match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0)
    //{cout << "islc.vertex.x " << islc.vertex.x << " match.flashGateTime " << match.flashGateTime << " match.flashPosition.x " << match.flashPosition.x <<  endl;} 
    if((islc.vertex.x*match.flashPosition.x>0) && match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return match.flashPosition.x;} 

  }
  return 0;
};

double bar_flash_gatetime (const caf::SRSpillProxy* spill,const caf::Proxy<caf::SRSlice>& islc){
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if((islc.vertex.x*match.flashPosition.x>0) && match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return (match.flashTime_us + spill->hdr.triggerinfo.trigger_within_gate);} 

  }
  return 0;
};


const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* spill){
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if(match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){return true;} 




  }
  return false;
});

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

    return chi2_ALG(dedx, rr, CALO_RR_MIN, CALO_RR_MAX);
}

PFPID id_pfp(const caf::Proxy<caf::SRSlice>& islc,
             int ipfp,
             double dist_cut)
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
    if (!pfp.parent_is_primary){ return PFPID::Unknown;}
    if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x) ||
        std::isnan(pfp.trk.len)){ return PFPID::Unknown;}

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
    if (min_dist > VTX_MAX_DIST){return PFPID::Unknown;}

    // =========================================================
    // CUT 2 — calorimetry availability
    // =========================================================
    /*
    if (pfp.trk.calo[DEFAULT_PLANE].nhit <= 2 ||
        pfp.trk.calo[DEFAULT_PLANE].points.empty()) {
        {MyFile << "No calo check - Unknown" << endl; return PFPID::Unknown;}
    }*/

    std::vector<double> chi2 = compute_chi2(islc, ipfp, DEFAULT_PLANE);
    if (chi2.size() < 2) {return PFPID::Unknown;}

    const double chi2_proton = chi2[static_cast<int>(Chi2PID::Proton)];

    // =========================================================
    // CUT 3 — good track-like objects
    // =========================================================
    if (pfp.trackScore >= PRIMARY_TRACK_SCORE) {

        // ---------- pion hypothesis ----------
        if (chi2_proton >= 100.0) {

            TVector3 p_pi(pfp.trk.dir.x * pfp.trk.rangeP.p_pion,
              pfp.trk.dir.y * pfp.trk.rangeP.p_pion,
              pfp.trk.dir.z * pfp.trk.rangeP.p_pion);

            double ke = kinetic_energy(PION_MASS, p_pi.Mag());

            if ((reco_vtx - start).Mag() < dist_cut && ke >= PION_KE_MIN)
                {return PFPID::Pion;}
        }

        // ---------- proton hypothesis ----------
        if (chi2_proton < 100.0) {

            TVector3 p_pr(pfp.trk.dir.x * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.y * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.z * pfp.trk.rangeP.p_proton);

            double ke = kinetic_energy(PROTON_MASS, p_pr.Mag());

            if ((reco_vtx - start).Mag() < dist_cut && ke >= PROTON_KE_MIN)
                {return PFPID::Proton;}
        }
    }

    // =========================================================
    // CUT 4 — marginal tracks (proton recovery)
    // =========================================================
    if (pfp.trackScore >= SECONDARY_TRACK_LOW &&
        pfp.trackScore <  PRIMARY_TRACK_SCORE &&
        chi2_proton < 100.0) {

        TVector3 p_pr(pfp.trk.dir.x * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.y * pfp.trk.rangeP.p_proton,
              pfp.trk.dir.z * pfp.trk.rangeP.p_proton);

        double ke = kinetic_energy(PROTON_MASS, p_pr.Mag());

        if ((reco_vtx - start).Mag() < dist_cut && ke >= PROTON_KE_MIN)
            {return PFPID::Proton;}
    }

    // =========================================================
    // CUT 5 — shower-like objects
    // =========================================================
    if (!(pfp.trackScore >= SECONDARY_TRACK_LOW &&
        pfp.trackScore <  PRIMARY_TRACK_SCORE &&
        chi2_proton < 100.0) && pfp.trackScore < PRIMARY_TRACK_SCORE){
    //if (pfp.trackScore < PRIMARY_TRACK_SCORE) {

        //if (pfp.shw.plane[DEFAULT_PLANE].nHits <= 2){MyFile << "Shw no calo - Unknown" << endl; return PFPID::Unknown;}
        if (std::isnan(pfp.shw.plane[DEFAULT_PLANE].energy)){return PFPID::Unknown;}

        double E = pfp.shw.plane[DEFAULT_PLANE].energy * 1000.0; // MeV

        if (E > PION_KE_MIN)
            {return PFPID::Shower;}
    }

    return PFPID::Unknown;
}


ParticleCounts count_particles(const caf::Proxy<caf::SRSlice>& islc,
                                int muon_index,
                                int dist_cut)
{
    ParticleCounts counts;

    for (std::size_t ipfp = 0; ipfp < islc.reco.npfp; ++ipfp) {

        if (static_cast<int>(ipfp) == muon_index){continue;}



        PFPID pid = static_cast<PFPID>(id_pfp(islc, ipfp, dist_cut));

        switch (pid) {
            case PFPID::Proton:  ++counts.protons;  break;
            case PFPID::Pion:    ++counts.pions;    break;
            case PFPID::Shower:  ++counts.showers;  break;
            default: break;
        }
    }

    return counts;
}


int find_muon(const caf::Proxy<caf::SRSlice>& islc, double dist_mucut)
{
    // =========================================================
    // Muon selection strategy:
    //   → longest primary contained track
    //   → muon-like dE/dx (χ²)
    //   → vertex-associated
    // =========================================================

    double max_length = -1.0;
    int muon_index = -1;

    const TVector3 reco_vtx(islc.vertex.x,
                            islc.vertex.y,
                            islc.vertex.z);

    for (std::size_t ipfp = 0; ipfp < islc.reco.npfp; ++ipfp) {

        const auto& pfp = islc.reco.pfp[ipfp];

        // =====================================================
        // CUT 0 — Basic track validity
        // =====================================================
        if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.len)) continue;
        if (pfp.trackScore < PRIMARY_TRACK_SCORE) continue;
        //if (pfp.trk.calo[DEFAULT_PLANE].nhit <= 2) continue;
        //if (pfp.trk.calo[DEFAULT_PLANE].points.empty()) continue;

        const TVector3 reco_start(pfp.trk.start.x,
                                  pfp.trk.start.y,
                                  pfp.trk.start.z);

        // =====================================================
        // CUT 1 — Vertex association
        // =====================================================
        if ((reco_vtx - reco_start).Mag() > dist_mucut) continue;

        // =====================================================
        // CUT 2 — Track quality
        // =====================================================
        if (pfp.trk.len < MIN_MUON_LENGTH) continue;
        if (!pfp.parent_is_primary) continue;

        // =====================================================
        // CUT 3 — Containment
        // =====================================================
        if (!isInContained(pfp.trk.end.x,
                           pfp.trk.end.y,
                           pfp.trk.end.z,
                           CONTAINMENT_CUT)) continue;

        // =====================================================
        // CUT 4 — Same cryostat
        // =====================================================
        if ((pfp.trk.end.x * islc.vertex.x) <= 0) continue;

        // =====================================================
        // CUT 5 — dE/dx PID (χ²)
        // =====================================================
        std::vector<double> chi2 = compute_chi2(islc, ipfp, DEFAULT_PLANE);

        if (chi2.size() < 2){ continue;}
        if(std::isnan(chi2[static_cast<int>(Chi2PID::Muon)]) || std::isnan(chi2[static_cast<int>(Chi2PID::Proton)]))continue;
        if (chi2[static_cast<int>(Chi2PID::Muon)]   > MAX_CHI2_MUON)   continue;
        if (chi2[static_cast<int>(Chi2PID::Proton)] < MIN_CHI2_PROTON) continue;

        // =====================================================
        // FINAL — choose longest muon-like track
        // =====================================================
        if (pfp.trk.len > max_length) {
            max_length = pfp.trk.len;
            muon_index = static_cast<int>(ipfp);
        }
    }

    return muon_index;
}



bool automatic_selection_1muNp_new(const caf::SRSpillProxy* sr,
                              const caf::Proxy<caf::SRSlice>& islc,
                              int dist_cut,
                              int cut_baryc)
{
    // =========================================================
    // CUT 0 — Sanity checks (reco & charge info must exist)
    // =========================================================
    if (std::isnan(islc.vertex.x) ||
        std::isnan(islc.vertex.y) ||
        std::isnan(islc.vertex.z) ||
        std::isnan(islc.charge_center.z)) {
        return false;
    }

    // =========================================================
    // CUT 1 — Fiducial volume
    // =========================================================
    if (!isInFV(islc.vertex.x, islc.vertex.y, islc.vertex.z)) {
        return false;
    }

    // =========================================================
    // CUT 2 — Barycenter–trigger consistency
    // =========================================================
    /*
    if (islc.barycenterFM.deltaZ_Trigger <= 0 ||
        islc.barycenterFM.deltaZ_Trigger >= cut_baryc) {
        return false;
    }
    */
   
    
      double barycenter=bar_flash(sr,islc);  //This is z!!!
	  double barycenter_x=bar_flash_x(sr,islc);

	  double bar_charge=islc.charge_center.z;
	  double delta=fabs(barycenter-bar_charge);
      // << "delta  " << delta << " barycenter Z " << barycenter << " X " << barycenter_x << " charge " << bar_charge << endl;
      if(barycenter_x*islc.vertex.x<=0){/*cout << "Other cryo " << endl;*/ return false;}
      //if(barycenter<=-10000){/*cout << "larger than 10k " << endl;*/return false;}
      //if(delta <=0 || delta>=100){/*cout << "Big delta or 0 " << endl;*/ return false;}

    

    // =========================================================
    // CUT 3 — Full containment
    // =========================================================
    if (!all_contained(islc)) {
        return false;
    }

    // =========================================================
    // CUT 3.1 — CRT PMT
    // =========================================================

    if (!kCRTPMTNeutrino(sr)) 
    {
        return false;
    }

    // =========================================================
    // CUT 4 — Exactly one reconstructed muon candidate
    // =========================================================
    const int muon_index = find_muon(islc, dist_cut);
    if (muon_index < 0) {
        return false;
    }

    // =========================================================
    // CUT 5 — Topology: 1μ + Np, no pions, no showers
    // =========================================================
    ParticleCounts counts = count_particles(islc, muon_index, dist_cut);
    
    //if(counts.protons > 0 && counts.pions   == 0 && counts.showers == 0){cout << "SELECTED " << endl;}
    // cout << " Protons " << counts.protons << " pions " << counts.pions << " showers " << counts.showers << endl;
    return (counts.protons > 0 &&
            counts.pions   == 0 &&
            counts.showers == 0);
}


double Neutrino_energy_reco_Np(const caf::Proxy<caf::SRSlice>& islc,
                               int ipfp_mu,
                               double dist_emucut)
{
    // =========================================================
    // CC Np energy reconstruction:
    //   Eν = Eμ + ΣTp + Eb
    // where Eb = 30.9 MeV effective argon binding
    // =========================================================

    if (ipfp_mu < 0 || ipfp_mu >= static_cast<int>(islc.reco.npfp))
        return -999.0;

    const auto& mu = islc.reco.pfp[ipfp_mu];

    // =========================================================
    // Muon energy
    // =========================================================
    TVector3 p_mu(mu.trk.dir.x * mu.trk.rangeP.p_muon,
                  mu.trk.dir.y * mu.trk.rangeP.p_muon,
                  mu.trk.dir.z * mu.trk.rangeP.p_muon);   // GeV

    const double p_mu_mag = p_mu.Mag();                    // GeV
    const double E_mu = std::sqrt((p_mu_mag * 1000.0)*(p_mu_mag * 1000.0)
                                   + MUON_MASS*MUON_MASS); // MeV

    // =========================================================
    // Proton kinetic energy sum
    // =========================================================
    double E_p_sum = 0.0; // MeV
    bool found_proton = false;

    for (std::size_t ipfp = 0; ipfp < islc.reco.npfp; ++ipfp) {

        if (static_cast<int>(ipfp) == ipfp_mu) continue;

        if (id_pfp(islc, ipfp, dist_emucut) != PFPID::Proton) continue;

        const auto& pfp = islc.reco.pfp[ipfp];

        TVector3 p_pr(pfp.trk.dir.x * pfp.trk.rangeP.p_proton,
                      pfp.trk.dir.y * pfp.trk.rangeP.p_proton,
                      pfp.trk.dir.z * pfp.trk.rangeP.p_proton); // GeV

        double T_p = kinetic_energy(PROTON_MASS, p_pr.Mag());

        E_p_sum += (T_p + PROTON_BINDING_ENERGY);
        found_proton = true;
    }

    if (!found_proton) return -999.0;

    // =========================================================
    // Neutrino energy (GeV)
    // =========================================================
    double E_nu = (E_mu + E_p_sum) / 1000.0;

    return E_nu;
}

double Neutrino_energy_visible_true_Np(const caf::SRSpillProxy* sr,
                                      const caf::Proxy<caf::SRTrueInteraction>& nu)
{

    // ----------------- Counters -----------------

    double E_p_sum = 0.0;
    double E_mu = 0.0;

    const int use_plane = DEFAULT_PLANE;

    // ----------------- Primary loop -----------------
    for (const auto& ipart : nu.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        TruthPID pid = classify_truth_pid(ipart.pdg); 
        double depE = ipart.plane[ipart.cryostat][use_plane].visE * 1000.0;

        // -------- Muons --------
        if (pid == TruthPID::Muon) {
            TVector3 p_mu(ipart.genp.x*1000.0,ipart.genp.y*1000.0,ipart.genp.z*1000.0); // MeV

            E_mu = std::sqrt(p_mu.Mag()*p_mu.Mag() + MUON_MASS*MUON_MASS); 
        }


        // -------- Protons --------
        if (pid == TruthPID::Proton) {
            for (const auto& d : sr->true_particles)
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;


            if (depE > PROTON_KE_MIN) {

            TVector3 p_pr(ipart.genp.x,ipart.genp.y,ipart.genp.z); // GeV

            double T_p = kinetic_energy(PROTON_MASS, p_pr.Mag());  // momentum in GeV!! 

        E_p_sum += (T_p + PROTON_BINDING_ENERGY);                
            }
        }
    }

        //std::cout << "Energy " << (E_p_sum+E_mu)/1000.0 << endl;
    return (E_p_sum+E_mu)/1000.0;
}


double Neutrino_energy_visible_true_Np_SLICE(const caf::SRSpillProxy* sr,
                                      const caf::Proxy<caf::SRSlice>& islc)
{

    // ----------------- Counters -----------------

    double E_p_sum = 0.0;
    double E_mu = 0.0;

    const int use_plane = DEFAULT_PLANE;

    // ----------------- Primary loop -----------------
    for (const auto& ipart : islc.truth.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        TruthPID pid = classify_truth_pid(ipart.pdg); 
        double depE = ipart.plane[ipart.cryostat][use_plane].visE * 1000.0;

        // -------- Muons --------
        if (pid == TruthPID::Muon) {
            TVector3 p_mu(ipart.genp.x*1000.0,ipart.genp.y*1000.0,ipart.genp.z*1000.0); // MeV

            E_mu = std::sqrt(p_mu.Mag()*p_mu.Mag() + MUON_MASS*MUON_MASS); 
        }


        // -------- Protons --------
        if (pid == TruthPID::Proton) {
            for (const auto& d : sr->true_particles)
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;


            if (depE > PROTON_KE_MIN) {

            TVector3 p_pr(ipart.genp.x,ipart.genp.y,ipart.genp.z); // GeV

            double T_p = kinetic_energy(PROTON_MASS, p_pr.Mag());  // momentum in GeV!! 

        E_p_sum += (T_p + PROTON_BINDING_ENERGY);                
            }
        }
    }

        //std::cout << "Energy " << (E_p_sum+E_mu)/1000.0 << endl;
    return (E_p_sum+E_mu)/1000.0;
}

int n_1mu1p = 0; int nu_1mu1p = 0;
int n_1muNp = 0; int nu_1muNp = 0;
int n_other = 0; int nu_other = 0;
int n_cosmic = 0; 
int n_invalid = 0; int nu_invalid = 0;



const SpillMultiVar kNeutrino_energy([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;
/*
    for ( auto const& nu : sr->mc.nu ){
        TruthClass cls_nu = classification_type_MC(sr,nu);

        switch (cls_nu) {
            case TruthClass::kOneMuOneP: ++nu_1mu1p; break;
            case TruthClass::kOneMuNp:   ++nu_1muNp; break;
            case TruthClass::kInvalid:   ++nu_invalid; break;
            default:                     ++nu_other; break;
        }

    }
    */

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(islc.reco.pfp[ipfp_mu].trk.end.y>125.0){
        TruthClass cls = classification_type_debug(sr, islc);
        
        if(ipfp_mu!=-1){
                /*
                MyFile_truth << " Muon end X " << islc.reco.pfp[ipfp_mu].trk.end.x << " true " << islc.reco.pfp[ipfp_mu].trk.truth.p.end.x<< endl;
                MyFile_truth << " Muon end Y " << islc.reco.pfp[ipfp_mu].trk.end.y << " true " << islc.reco.pfp[ipfp_mu].trk.truth.p.end.y<< endl;
                MyFile_truth << " Muon end Z " << islc.reco.pfp[ipfp_mu].trk.end.z << " true " << islc.reco.pfp[ipfp_mu].trk.truth.p.end.z<< endl;
*/
                  //MyFile2 << sr->hdr.run << " " << sr->hdr.evt << " " << ipfp_mu << " " << Neutrino_energy_reco_Np(islc,ipfp_mu,10) << " class " << TruthClassToString(cls) << endl; 
        switch (cls) {
            case TruthClass::kOneMuOneP: {++n_1mu1p; /*MyFile_truth << " 1mu1p interaction "<< endl;*/ break;}
            case TruthClass::kOneMuNp:   {++n_1muNp;/* MyFile_truth << " 1muNp interaction "<< endl;*/ break;}
            case TruthClass::kCosmic:    {++n_cosmic;/* MyFile_truth << " Cosmic "<< endl;*/ break;}
            //case TruthClass::kInvalid:   ++n_invalid; break;
            //default:                     ++n_other; break;
            default:
                // none of the above → check other conditions
                if (!isInFV(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z)) {
                   //  MyFile_truth << " OOFV interaction " << endl;
                }
                else if (std::abs(islc.truth.pdg) == 12) {  //nue
                   //  MyFile_truth << " nue interaction "<< endl;
                }   
                else if (islc.truth.isnc) {
                   //  MyFile_truth << " NC interaction "<< endl;
                }             
                else if (islc.truth.iscc && islc.truth.genie_mode == 0) {   //QE
                   //  MyFile_truth << " CCQE interaction "<< endl;
                }                 
                else if (islc.truth.iscc && islc.truth.genie_mode == 1) {   //RES
                   // MyFile_truth << " CCRES interaction "<< endl;
                }                     
                else if (islc.truth.iscc && islc.truth.genie_mode == 2) {   //DIS
                   //  MyFile_truth << " CCDIS interaction "<< endl;
                }        
                else if (islc.truth.iscc && (islc.truth.genie_mode == 3 || islc.truth.genie_mode == 4)) {   //COH
                   //  MyFile_truth << " CCCOH interaction "<< endl;
                } 
                else if (islc.truth.iscc && islc.truth.genie_mode == 10) {   //MEC
                   //  MyFile_truth << " CCMEC interaction "<< endl;
                }  
                else {
                    // optional fallback
                   //  MyFile_truth << " OTHER interaction "<< endl;
                }
                break;

        }

        vector_active.push_back(Neutrino_energy_reco_Np(islc,ipfp_mu,10)); }
    }
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});



std::vector<std::pair<int,int>> evt_pairs = {
    {10064, 55977},{9367,62724},{10064,88063}
};

bool contains_pair(const std::vector<std::pair<int,int>>& v, int a, int b)
{
    return std::any_of(v.begin(), v.end(),
                       [&](auto const& p){ return p.first == a && p.second == b; });
}


const SpillCut kSelectEvt([](const caf::SRSpillProxy* spill){
    if(contains_pair(evt_pairs, spill->hdr.run, spill->hdr.evt))
    {
        return true;
    }
    return false;
});


const SpillVar kEvt_M([](const caf::SRSpillProxy* sr) {return sr->hdr.evt;});
const SpillVar kSelectRun_M([](const caf::SRSpillProxy* sr) {return sr->hdr.run;});
