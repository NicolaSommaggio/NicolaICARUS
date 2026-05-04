
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
#include <numeric>

#include <fstream>
#include <iostream>
#include <sstream>

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

/* CNAF
TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");
*/

//FNAL
TFile* file = TFile::Open("/exp/icarus/data/users/nsommagg/RefCurvesChi2.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* spill){
    for(const auto& match: spill->crtpmt_matches) 
    {
        //Define the interval depending on Data or MC files
        double min_time =-1; double max_time =-1;
        if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
        if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
        //if(match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return true;} 
    }
    //return false;
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
			  if(itrue.end.x!=-9999 && itrue.end.y!=-9999 && itrue.end.z!=-9999){if(isInContained(ipart.end.x,ipart.end.y,ipart.end.z,5.)==false){return false;}}
                    }                             
                        }
                        }
                    
                        } 
 
                }  
return true;
}

std::vector<double> chi2_ALG(std::vector<double> &dEdx,std::vector<double> &RR, double rr_min, double rr_max)
{

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


struct SliceDrawInfo
{
    int run;
    int evt;
    std::vector<int> v_ipfp_mu;
    std::vector<int> v_ipfp_pro;
    std::vector<int> v_ipfp_pi;
    std::vector<std::vector<std::array<double,3>>> tracks_coordinate;
    std::vector<std::vector<double>> tracks_dedx;
    //std::vector<std::vector<std::array<double,3>>> showers_coordinate;
    //std::vector<std::vector<double>> showers_dedx;
    std::vector<std::pair<int,std::array<double,6>>> true_particles;
    std::vector<int> ipfp;
    std::array<double,3> vertex;
    std::vector<double> neutrino_interaction;
};

int find_idx(double rr)
{
    if(rr==25.)return 47;
    double a = std::round(rr);
    double max=-1;
    double min=-1;
    if(a>rr)
    {
        min=a-0.5;
        max=a;
    }
    else
    {
        min=a;
        max=a+0.5;
    }
    
    int n = (max-1.5)/0.5;

    return n;
}

int find_best_plane(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
    //there should be at least 1 hit with rr between 1 cm and 25 cm and with the dedx<30 MeV/cm
    //if the condition is satified in collection, collection is the best plane regardless
    int bestplane = 2;
    bool hasValidHits=false;

    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()>0)
    {
        //check if it has valid hits for for computing the likelihood
        for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
        {
            if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
            {
                hasValidHits=true;
            }
        }   
    }
    
    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()==0 || hasValidHits==false)
    {
        bestplane = islc.reco.pfp[ipfp].trk.bestplane;
        if(bestplane==-1)return -1;

        //check if it has valid hits for for computing the likelihood
        for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
        {
            if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
            {
              hasValidHits=true;
            }
        }

        if(hasValidHits==false)return -1;
    }

    return bestplane;
}


std::array<std::vector<TH1D*>,6> load_prob_densities(std::string plane, TFile *inputf)
{
  std::vector<TH1D*> probDensitiesMU_class0;
    std::vector<TH1D*> probDensitiesMU_class1;
    std::vector<TH1D*> probDensitiesPRO_class2;
    std::vector<TH1D*> probDensitiesPRO_class3;
    std::vector<TH1D*> probDensitiesPI_class4;
    std::vector<TH1D*> probDensitiesPI_class5;

    std::array<std::vector<TH1D*>,6> probDensities={
        probDensitiesMU_class0,
        probDensitiesMU_class1,
        probDensitiesPRO_class2,
        probDensitiesPRO_class3,
        probDensitiesPI_class4,
        probDensitiesPI_class5
    };

    std::vector<std::pair<std::string,int>> classes = {{"muon_class0",0}, {"muon_class1",1}, {"proton_class2",2}, {"proton_class3",3}, {"pion_class4",4}, {"pion_class5",5}};

    for(const auto &clas : classes)
    {
        for(double i=1.5; i<=25.0; i+=0.5)
        {
            TDirectory *dclass = (TDirectory*)inputf->Get(clas.first.c_str());
            TDirectory *d = (TDirectory*)dclass->Get(Form("%.1f",i));

            TH1D *dEdx_coll = (TH1D*)d->Get(Form("dEdx_%s_rr_%.1f", plane.c_str(), i));
            probDensities[clas.second].push_back(dEdx_coll);
        }
    }

    return probDensities;
}

std::vector<double> likelihood(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp, std::array<std::vector<TH1D*>,6> &coll_pd, std::array<std::vector<TH1D*>,6> &ind1_pd, std::array<std::vector<TH1D*>,6> &ind2_pd)
{    
  auto probDensities = coll_pd;

  std::vector<double> likelihood_return;
  for(int i=0; i<15; i++){likelihood_return.push_back(-1);}
  
  std::vector<double> likelihood_ratios;
            
  std::array<double,6> lkl; //it contains the likelihoods in the 6 hypotheses for the current track

  bool hasValidHits=false;

  //ciclo sulla particle hypothesis
  for(int j=0; j<6; j++)
  {
    auto density = probDensities[j];
    double log_lkh=0;
        
    for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit)
    {
      if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx<30.)
      {
        hasValidHits=true;
        int idx=find_idx(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
        int bin = density[idx]->FindBin(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
        log_lkh = log_lkh + (-1*std::log(density[idx]->GetBinContent(bin)/density[idx]->GetBinWidth(bin)));
        //if(std::isinf(log_lkh)){ cout << dat.track.rr->at(hit) << " " << dat.track.dE->at(hit) << " | " << density[idx]->GetBinContent(bin) << " " << std::log(density[idx]->GetBinContent(bin)) << " " << density[idx]->GetBinWidth(bin) << " | " << par << " " << subclass << " " << idx << "(" << dat.track.rr->at(hit) << ")" << " " << bin << endl; break;}
      }
    }
    lkl[j]=log_lkh;
  }//cycle over particle hypothesis

  if(!hasValidHits)
  {
    int use_plane = -1;
    if(islc.reco.pfp[ipfp].trk.bestplane==0 || islc.reco.pfp[ipfp].trk.bestplane==1 || islc.reco.pfp[ipfp].trk.bestplane==2){use_plane=islc.reco.pfp[ipfp].trk.bestplane;}
    else{use_plane=2;}

    if(use_plane == 0){probDensities = ind1_pd;}
    else if(use_plane == 1){probDensities = ind2_pd;}
    else if(use_plane == 2){probDensities = coll_pd;}

    for(int j=0; j<6; j++)
    { 
      auto density = probDensities[j];
      double log_lkh=0;

      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit)
      {
        if(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx<30.)
        {
          hasValidHits=true;
          int idx=find_idx(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
          int bin = density[idx]->FindBin(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
          log_lkh = log_lkh + (-1*std::log(density[idx]->GetBinContent(bin)/density[idx]->GetBinWidth(bin)));
          //if(std::isinf(log_lkh)){ cout << dat.track.rr->at(hit) << " " << dat.track.dE->at(hit) << " | " << density[idx]->GetBinContent(bin) << " " << std::log(density[idx]->GetBinContent(bin)) << " " << density[idx]->GetBinWidth(bin) << " | " << par << " " << subclass << " " << idx << "(" << dat.track.rr->at(hit) << ")" << " " << bin << endl; break;}
        }
      }
      lkl[j]=log_lkh;
    }//cycle over particle hypothesis
  }

  if(!hasValidHits){return likelihood_return;}

  //GETTING LIKELIHOOD RATIOS FOR EACH TRACK
  for(int k=0; k<6; k++)
  {
    for(int t=k+1; t<6; t++)
    {
      //if(isnan(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90))cout << "nan value ecountered " << lkl[k] << " " << lkl[t] << endl;
      likelihood_ratios.push_back(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90);
    }
  }  
  return likelihood_ratios;
}

double compute_depE(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp, int plane)
{
    double dep_E = 0;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[plane].points.size(); ++ihit )
    {
        if(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr <= 5.)
        {   
            dep_E = dep_E + islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].pitch;
        }
    }
    return dep_E;
}

std::vector<double> compute_daughter_vars(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
    int daughter_max_hits = 0;
    int daughter_id = -1;
    std::vector<double> daughter_rr;
    std::vector<double> daughter_dedx;
    double daughter_depE = 0;
    std::array<double,3> daughter_direction={-1,-1,-1};
    double angle_end = -1;

    for(int daughter=0; daughter<int(islc.reco.pfp[ipfp].ndaughters); daughter++)
    {  
        int d = islc.reco.pfp[ipfp].daughters[daughter];
    
        std::vector<double> temp_dedx;
        std::vector<double> temp_rr;
        double temp_depE = 0;
        std::array<double,3> temp_direction;
        for (std::size_t jpfp(0); jpfp < islc.reco.npfp; ++jpfp)
        {
            if(islc.reco.pfp[jpfp].id!=d)continue;

            int d_bestplane = islc.reco.pfp[jpfp].trk.bestplane;

            if(d_bestplane==-1)continue;

            for ( std::size_t ihit(0); ihit < islc.reco.pfp[jpfp].trk.calo[d_bestplane].points.size(); ++ihit )
            {
              if(islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].rr <= 5)
              {
                temp_rr.push_back(islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].rr);
                temp_dedx.push_back(islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].dedx);
                temp_depE = temp_depE + islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].dedx*islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].pitch;
              }
            }
            temp_direction = {islc.reco.pfp[jpfp].trk.dir.x, islc.reco.pfp[jpfp].trk.dir.y, islc.reco.pfp[jpfp].trk.dir.z};

            break;
        }

        if(temp_rr.empty())continue;
        if((int)temp_rr.size()>daughter_max_hits)
        {
          daughter_max_hits=(int)temp_rr.size();
          daughter_id = daughter;
          daughter_rr = temp_rr;
          daughter_dedx = temp_dedx;
          daughter_depE = temp_depE;
          daughter_direction = temp_direction;
        }
    }

    std::array<double,3> track_direction_end = {islc.reco.pfp[ipfp].trk.dir_end.x, islc.reco.pfp[ipfp].trk.dir_end.y, islc.reco.pfp[ipfp].trk.dir_end.z};
    angle_end = std::acos(track_direction_end[0]*daughter_direction[0] + track_direction_end[1]*daughter_direction[1] + track_direction_end[2]*daughter_direction[2]);

    if(daughter_id==-1)return {-1,-1};

    return {daughter_depE,angle_end};
}

TFile * f_prob_densities_coll = TFile::Open("/exp/icarus/data/users/nsommagg/PID_include_files/HISTO_prob_densities_6_classes_30percent_COLL.root", "READ");
TFile * f_prob_densities_ind1 = TFile::Open("/exp/icarus/data/users/nsommagg/PID_include_files/HISTO_prob_densities_6_classes_30percent_IND1.root", "READ");
TFile * f_prob_densities_ind2 = TFile::Open("/exp/icarus/data/users/nsommagg/PID_include_files/HISTO_prob_densities_6_classes_30percent_IND2.root", "READ");

std::array<std::vector<TH1D*>,6> prob_d_coll = load_prob_densities("coll",f_prob_densities_coll);
std::array<std::vector<TH1D*>,6> prob_d_ind1 = load_prob_densities("ind1",f_prob_densities_ind1);
std::array<std::vector<TH1D*>,6> prob_d_ind2 = load_prob_densities("ind2",f_prob_densities_ind2);


//PID model ------------------------------------------------------------------
std::string path_BDT_model = "/exp/icarus/data/users/nsommagg/PID_include_files/GBDT.txt";
#include "/exp/icarus/data/users/nsommagg/PID_include_files/GBDT_model_inport.h"
GBDTModel model = load_model(path_BDT_model.c_str());
//----------------------------------------------------------------------------

int PIDclass(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
  //-1 : unclassified
  //0 : muon rising
  //1 : muon mip
  //2 : proton rising
  //3 : proton interacting
  //4 : pion rising
  //5 : pion interacting 

  int bestplane = find_best_plane(islc,ipfp);

  if(bestplane == -1)return -1; //no valid hits (RR < 25 cm and RR > 1 cm and dEdx < 30 MeV/cm) in neither of the 3 planes

  std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
  double depE = compute_depE(islc,ipfp,bestplane);
  std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

  std::vector<double> track_features;

  track_features.insert(track_features.end(), lr.begin(), lr.end());
  track_features.push_back(depE);
  track_features.insert(track_features.end(), dvars.begin(), dvars.end());

  int prediction = model.predict_class(track_features); 
  
  return prediction;
}

std::vector<double> PIDproba(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
  //it return a 6D vector containing the probabilities for the track to be each of the 6 classes 

  int bestplane = find_best_plane(islc,ipfp);

  if(bestplane == -1)return {-1.,-1.,-1.,-1.,-1.,-1.}; //no valid hits (RR < 25 cm and RR > 1 cm and dEdx < 30 MeV/cm) in neither of the 3 planes

  std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
  double depE = compute_depE(islc,ipfp,bestplane);
  std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

  std::vector<double> track_features;

  track_features.insert(track_features.end(), lr.begin(), lr.end());
  track_features.push_back(depE);
  track_features.insert(track_features.end(), dvars.begin(), dvars.end());

  std::vector<double> prediction_proba = model.predict_proba(track_features); 
  
  return prediction_proba;
}

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


SliceDrawInfo GetSliceDrawInfo(const caf::Proxy<caf::SRSlice>& islc, std::vector<int> v_ipfp_mu, std::vector<int> v_ipfp_pro, std::vector<int> v_ipfp_pi, std::vector<std::pair<int,std::array<double,6>>> true_p)
{
    SliceDrawInfo thislice;

    thislice.true_particles=true_p;
    std::vector<std::vector<std::array<double,3>>> particles;
    std::vector<std::vector<double>> particles_dedx;
    thislice.vertex[0] = islc.vertex.x;
    thislice.vertex[1] = islc.vertex.y;
    thislice.vertex[2] = islc.vertex.z;

    thislice.v_ipfp_mu = v_ipfp_mu;
    thislice.v_ipfp_pro = v_ipfp_pro;
    thislice.v_ipfp_pi = v_ipfp_pi;

    for(std::size_t ipfp(0); ipfp<islc.reco.npfp; ipfp++)
    {
        int plane = find_best_plane(islc,ipfp);
        if(plane==-1)plane=2;

        thislice.ipfp.push_back(ipfp);

        std::vector<std::array<double,3>> track;
        std::vector<double> trackdedx;
        //std::vector<std::array<double,3>> shower;
        //std::vector<double> showerdedx;


        for(std::size_t ihit(0); ihit<islc.reco.pfp[ipfp].trk.calo[plane].points.size(); ihit++)
        {
            std::array<double,3> coordinate = {islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].x, islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].y, islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].z};
            track.push_back(coordinate);
            trackdedx.push_back(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx);
        }
        
        particles.push_back(track);
        particles_dedx.push_back(trackdedx);

        //int shower_bestplane = islc.reco.pfp[ipfp].shw.bestplane_for_dedx;
        //for(std::size_t ihit(0); ihit<islc.reco.pfp[ipfp].shw.plane[shower_bestplane].points.size(); ihit++)
   
    }

    thislice.tracks_coordinate = particles;
    thislice.tracks_dedx = particles_dedx;

    return thislice;
}


std::vector<std::pair<int,std::array<double,6>>> GetTrueParticlesInfo(const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc)
{
    //GET THE START AND END OF THE TRUE PRIMARIES AND SECONDARIES IN THE SLICE
    std::vector<std::pair<int,std::array<double,6>>> true_particles_inthislice;
    for(auto const & prim : islc.truth.prim)
    {

        if(!(std::abs(prim.pdg)==13 || std::abs(prim.pdg)==2212 || std::abs(prim.pdg)==211 || std::abs(prim.pdg)==11 || std::abs(prim.pdg)==111 || std::abs(prim.pdg)==321 || prim.pdg==22 ))continue;

         double dep_E=0;
        //SKIP LOW ENERGY PROTONS
        if(prim.cryostat!=-1){
        if(abs(prim.pdg)== 2212)
        {                    
            if(prim.daughters.size()>0)
            {
                for ( auto const& itrued : sr->true_particles )
                {
                    //sum depE daughters 
                    if(int(itrued.parent)==prim.G4ID ){dep_E+=itrued.plane[prim.cryostat][2].visE*1000;}
                }
            }
   
            dep_E += prim.plane[prim.cryostat][2].visE*1000;
        }
        //if(abs(prim.pdg)==2212 && dep_E<50.)continue;

        //SKIP LOW ENERGY PHOTONS
        if(abs(prim.pdg) == 22)
        {                    
            if(prim.daughters.size()>0)
            {
                for ( auto const& itrue : sr->true_particles )
                {
                    //sum depE daughters 
                    if(int(itrue.parent)==prim.G4ID ){dep_E+=itrue.plane[prim.cryostat][2].visE*1000;}
                }
            }
            dep_E += prim.plane[prim.cryostat][2].visE*1000;
        } 
        if(abs(prim.pdg)== 22 && dep_E<25.0)continue; 
        }

        std::pair<int,std::array<double,6>> thisparticle;
        thisparticle.first=prim.pdg;
        std::array<double,6> coordinate_temp={prim.start.x, prim.start.y, prim.start.z, prim.end.x, prim.end.y, prim.end.z};
        thisparticle.second=coordinate_temp;
        true_particles_inthislice.push_back(thisparticle);

        //GET SECONDARIES INFO
        for(auto const & itrue : sr->true_particles)
        {
            if(int(itrue.parent)==prim.G4ID && (std::abs(itrue.pdg)==13 || std::abs(itrue.pdg)==2212 || std::abs(itrue.pdg)==211 || std::abs(itrue.pdg)==11 || std::abs(itrue.pdg)==111 || std::abs(itrue.pdg)==321 /*|| itrue.pdg==22*/  ))
            {
                std::pair<int,std::array<double,6>> thisparticle_daughter;
                thisparticle_daughter.first=itrue.pdg;
                std::array<double,6> coordinate_temp_daughter={itrue.start.x, itrue.start.y, itrue.start.z, itrue.end.x, itrue.end.y, itrue.end.z};
                thisparticle_daughter.second=coordinate_temp_daughter;
                true_particles_inthislice.push_back(thisparticle_daughter);
            }
        }           
    }
    return true_particles_inthislice;
}

double mediana(std::vector<double> dummy)
{
    std::sort(dummy.begin(), dummy.end());
    double mediana=0;
    int idx=0;
    if(int(dummy.size())%2==1)
    {
        idx=(dummy.size()-1)/2;
        mediana=dummy.at(idx);
    }
    if(dummy.size()%2==0)
    {
        idx =dummy.size()/2-1;
        mediana=( dummy.at(dummy.size()/2) + dummy.at(idx) )/2.;
    }
    return mediana;
}

std::vector<std::pair<double,double>> rollingMedian(const caf::Proxy<caf::SRSlice>& islc, int ipfp)
{
    std::vector<double> dedx_temp;
    std::vector<double> rr_temp;
    std::vector<std::pair<double, double>> rm;
    for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit)
    {
        if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr<30)
        {
            dedx_temp.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
            rr_temp.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
        }
    }
    if(dedx_temp.size()<5)return {{-1,-1}};
    
    for(int i=3; i<int(dedx_temp.size())-3; i++)
    {
        std::vector<double> dummy = {dedx_temp[i-4], dedx_temp[i-3], dedx_temp[i-2], dedx_temp[i-1], dedx_temp[i] ,dedx_temp[i+1],dedx_temp[i+2], dedx_temp[i+3], dedx_temp[i+4]};
        double media_rr = (rr_temp[i-4] +  rr_temp[i-3] + rr_temp[i-2] + rr_temp[i-1] + rr_temp[i] + rr_temp[i+1] + rr_temp[i+2] + rr_temp[i+3] + rr_temp[i+4])/9.;
        rm.push_back({media_rr,mediana(dummy)});
    }
    return rm;
}

std::vector<SliceDrawInfo> slices3D;


ofstream dumpInfo("logSelection.txt");
int slice_counter=0;

std::vector<std::vector<std::pair<double,double>>> dumpRM;
std::vector<double> dumpMedian;
std::vector<std::vector<double>> dedx;
std::vector<std::vector<double>> rr;
std::vector<std::vector<double>> x;

//int RUN_DA_GUARDARE=-1;
//int EVT_DA_GUARDARE=-1;
//int PLANE_DA_GUARDARE=-1;
std::vector<int> runList;
std::vector<int> evtList;
std::vector<int> sliceList;
std::vector<std::string> filenames;


int found_events = 0;
const SpillMultiVar disegna([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    if(found_events >= 10)return vector_active;

    bool correct_file = false;
    for(const auto &f : filenames)
    {
        if(std::string(sr->hdr.sourceName) == f)correct_file = true;
    }
    
    if(correct_file)
    {    

    /*
    //bool is_true_signal=false;
    if(!((int)sr->hdr.run==RUN_DA_GUARDARE && (int)sr->hdr.evt==EVT_DA_GUARDARE)) return {}; 
    cout << sr->nslc << endl;
    */ 

    if(sr->hdr.run>0){

    int slice_number = -1;
    for (auto const& islc : sr->slc)
    {     
        slice_number ++ ;
        bool match_run_evt=false;
        for(int i=0; i<(int)runList.size(); i++)
        {
            if(runList[i]==(int)sr->hdr.run && evtList[i]==(int)sr->hdr.evt && slice_number==sliceList[i])match_run_evt=true;
        }
        if(!match_run_evt)continue;

        found_events++;

        dumpInfo << "*** RUN " << sr->hdr.run << " *** EVT " << sr->hdr.evt << " ***" << endl;

        //PROTON kinetic energy
        /* 
        std::vector<double> Ek_pro;
        for(auto const & prim : islc.truth.prim)
        {
            if(std::abs(prim.pdg)==2212)
            {
                double dep_E=0;
                if(prim.cryostat!=-1)
                {
                    if(abs(prim.pdg)== 2212)
                    {                    
                        if(prim.daughters.size()>0)
                        {
                            for ( auto const& itrued : sr->true_particles )
                            {
                                //sum depE daughters 
                                if(int(itrued.parent)==prim.G4ID ){dep_E+=itrued.plane[prim.cryostat][2].visE*1000;}
                            }
                        }
   
                        dep_E += prim.plane[prim.cryostat][2].visE*1000;
                    }
                }
                Ek_pro.push_back(dep_E);
                dumpInfo << endl << "MC proton: Ek= " << dep_E << endl;
                dumpInfo << "start: " << prim.start.x << " " << prim.start.y << " " << prim.start.z << endl;
                dumpInfo << "end: " << prim.end.x << " " << prim.end.y << " " << prim.end.z << endl;
            }
        }

        auto at_max_depE = std::max_element(Ek_pro.begin(), Ek_pro.end());
        double max_depE = -1;
        if (at_max_depE != Ek_pro.end()) {max_depE = *at_max_depE;}
        //cout << sr->hdr.run << " " << sr->hdr.evt << " " << slice_number << " " << Form("vertex=(%f,%f,%f)",(double)islc.vertex.x,(double)islc.vertex.y,(double)islc.vertex.z) << " " << max_depE << endl;
        cout << max_depE << endl;
        */

        cout << "evento trovato " << (islc.is_clear_cosmic ? "clear_cosmic" : Form("vertex=(%f,%f,%f)",(double)islc.vertex.x,(double)islc.vertex.y,(double)islc.vertex.z)) << endl;

        //if(sr->hdr.evt % 2 !=0)continue;

        bool is_true_signal=false;
        int ipfp_mu=-1;
        std::vector<int> v_ipfp_mu;
        std::vector<int> v_ipfp_pro;
        std::vector<int> v_ipfp_pi;

        slice_counter+=1;

        //ipfp_mu=find_truth_muon(islc,10);
        ipfp_mu = find_muon_newPID(islc,10);

        //if(ipfp_mu==-1)continue;

        v_ipfp_mu.push_back(ipfp_mu);
	    for(std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
		{
            if(int(ipfp)==ipfp_mu)continue;
		    //if(id_pfp_truth(islc, ipfp, 10,2)==1){v_ipfp_pro.push_back(int(ipfp));}
            //if(id_pfp_truth(islc, ipfp, 10,2)==2){v_ipfp_pi.push_back(int(ipfp));}
            if(id_pfp_newPID(islc,ipfp,10)==1){v_ipfp_pro.push_back(int(ipfp));}
            if(id_pfp_newPID(islc,ipfp,10)==2){v_ipfp_pi.push_back(int(ipfp));}
        }
        
        //cout << sr->hdr.run << " " << sr->hdr.evt << " " << v_ipfp_pro.size() << endl;

        //rolling median
        //std::vector<std::pair<double,double>> rm_thislice={{-1,-1}};
        //if(ipfp_mu!=-1)
        //{
        //    dumpRM.push_back(rollingMedian(islc,ipfp_mu));
        //    rm_thislice=rollingMedian(islc,ipfp_mu);
        //}
        ////median
        //double median=-1;
        //std::vector<double> dummy;
        //std::vector<double> dedx_temp;
        //std::vector<double> rr_temp;
        //std::vector<double> x_temp;
        //if(ipfp_mu!=-1)
        //{
        //    for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit) 
        //    {
        //        if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr<30)
        //        {
        //            dedx_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx);
        //            rr_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr);
        //            x_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x);
        //        }
        //        if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr<5.) 
        //        {
        //            dummy.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx);
        //        }
        //    }
        //    if(dummy.size()!=0) {dumpMedian.push_back(mediana(dummy)); median=mediana(dummy);}
        //    else {dumpMedian.push_back(-1);}
        //    dedx.push_back(dedx_temp);
        //    rr.push_back(rr_temp);
        //    x.push_back(x_temp);
        //}

        int max_length_pfp=-1;
        double max_length=-1;

        

        for(std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp){if(islc.reco.pfp[ipfp].trk.len>max_length){max_length_pfp=ipfp; max_length=islc.reco.pfp[ipfp].trk.len;}}

        for(std::size_t ipfp(1); ipfp < islc.reco.npfp; ++ipfp)
        {
            TVector3 RecoVtx;
            RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
            TVector3 RecoStart;
            RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
            TVector3 RecoEnd;
            RecoEnd.SetXYZ(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);
            TVector3 Start_mom_v_proton;
            TVector3 Start_mom_v_pion;

            dumpInfo << "new PFP new PFP new PFP new PFP new PFP new PFP new PFP new PFP new PFP new PFP new PFP new PFP" << endl << endl;

            dumpInfo << "start coordinates " << islc.reco.pfp[ipfp].trk.start.x << " " << islc.reco.pfp[ipfp].trk.start.y << " " << islc.reco.pfp[ipfp].trk.start.z << endl;
            dumpInfo << "end coordinates " << islc.reco.pfp[ipfp].trk.end.x << " " << islc.reco.pfp[ipfp].trk.end.y << " " << islc.reco.pfp[ipfp].trk.end.z << endl;
            dumpInfo << "collection ";
            for (std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit ){dumpInfo << islc.reco.pfp[ipfp].trk.calo[2].points[ihit].x << " ";}
            dumpInfo << endl;
            dumpInfo << "induction 1 ";
            for (std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[0].points.size(); ++ihit ){dumpInfo << islc.reco.pfp[ipfp].trk.calo[0].points[ihit].x << " ";}
            dumpInfo << endl;
            dumpInfo << "induction 2 ";
            for (std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[1].points.size(); ++ihit ){dumpInfo << islc.reco.pfp[ipfp].trk.calo[1].points[ihit].x << " ";}
            dumpInfo << endl;
            dumpInfo << "trackscore " << islc.reco.pfp[ipfp].trackScore << endl;
            if(int(ipfp)==max_length_pfp)dumpInfo << "is the longest, with length " << islc.reco.pfp[ipfp].trk.len << endl;
            else dumpInfo << "track length " << islc.reco.pfp[ipfp].trk.len << " track distance " << (RecoEnd - RecoStart).Mag() << endl;
            dumpInfo << "vertex - start distance " << (RecoVtx-RecoStart).Mag() << " vertex - end distance " << (RecoVtx-RecoEnd).Mag() << endl;
            dumpInfo << "pdg code " << islc.reco.pfp[ipfp].trk.truth.p.pdg << endl;
            dumpInfo << "is end contained " << isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) << endl;
            dumpInfo << "is start contained " << isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,5.0) << endl;
            dumpInfo << "is the same TPC as vertex " << ((islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0) << endl;
            dumpInfo << "parent is primary " << islc.reco.pfp[ipfp].parent_is_primary << endl;
            Start_mom_v_proton.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
            Start_mom_v_pion.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);
            dumpInfo << "energy lost as proton " << sqrt(pow(938.3,2)+pow(Start_mom_v_proton.Mag()*1000,2))-938.3 << " energy lost as pion " << sqrt(pow(139.570,2)+pow(Start_mom_v_pion.Mag()*1000,2))-139.570 << endl;
     
            //if(int(ipfp)==ipfp_mu)
            //{
            //    dumpInfo << "median " << median << endl;
            //    dumpInfo << "rolling median ";
            //    for(int i=0; i<int(rm_thislice.size()); i++) dumpInfo << rm_thislice[i].first << " ";
            //    dumpInfo << endl;
            //    for(int i=0; i<int(rm_thislice.size()); i++) dumpInfo << rm_thislice[i].second << " ";
            //    dumpInfo << endl;
            //}
            
            //if(int(ipfp)==ipfp_mu) dumpInfo << "MUON" << endl << endl;
            //else if(id_pfp_truth(islc, ipfp, 10,2)==1) dumpInfo << "PROTON" << endl << endl;
            //else if(id_pfp_truth(islc, ipfp, 10,2)==2) dumpInfo << "PION" << endl << endl;
            //else dumpInfo << "OTHER" << endl << endl;

        }

        dumpInfo << endl;
        dumpInfo << "*********************************** end slice ***********************************" << endl;
            
        //INFORMAZIONI PER DISEGNARE LA SLICE 3D
        std::vector<double> nu_int_temp;
        nu_int_temp.push_back(islc.truth.E);
        nu_int_temp.push_back(islc.truth.targetPDG);
        nu_int_temp.push_back(islc.truth.Q2);
        std::vector<std::pair<int,std::array<double,6>>> true_particles=GetTrueParticlesInfo(sr,islc);
        //SliceDrawInfo drawInfo = GetSliceDrawInfo(islc,v_ipfp_mu,v_ipfp_pro,v_ipfp_pi,true_particles,PLANE_DA_GUARDARE);
        SliceDrawInfo drawInfo = GetSliceDrawInfo(islc,v_ipfp_mu,v_ipfp_pro,v_ipfp_pi,true_particles);
        drawInfo.neutrino_interaction=nu_int_temp;
        drawInfo.run = sr->hdr.run;
        drawInfo.evt = sr->hdr.evt;
        slices3D.push_back(drawInfo);

    }

    }

    }
    return vector_active;
});
