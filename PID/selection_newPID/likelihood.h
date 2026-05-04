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

// OLD VERSION
const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* spill){
  //for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    //double min_time =-1; double max_time =-1;
    //if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    //if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    //if(match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return true;} 
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

        if(islc.reco.pfp[ipfp].trackScore<0.5)continue;

        if(islc.reco.pfp[ipfp].trk.len>max_length && ((RecoVtx-RecoStart).Mag()<dist_mucut) && islc.reco.pfp[ipfp].trk.len>50 && 
        std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==13 && 
        isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) && 
        isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,5.0) &&
        (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && islc.reco.pfp[ipfp].parent_is_primary){
        max_length=islc.reco.pfp[ipfp].trk.len;
        ipfp_mu=ipfp;
            }
        }//loop of pfp to find muon
return ipfp_mu;
}


int id_pfp_truth ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_cut) { 
    //return 1 PROTONS
    //return 2 PIONS
    //return 3 SHOWER
    //return 9 other -> nan, not primary, too far, below energy threshold... 

    if(
        !isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) //||
        //!isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,5.0)
      )return 9;

    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);

    if (!(islc.reco.pfp[ipfp].parent_is_primary ))return 9;
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
 
    TVector3 start(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
    TVector3 end(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);

    double min_dist = ((start-RecoVtx).Mag()< (end-RecoVtx).Mag() ? (start-RecoVtx).Mag() : (end-RecoVtx).Mag());

    if(min_dist>50.0) return 9;
    
    if(islc.reco.pfp[ipfp].trackScore>=0.5)
    {
      if (std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
      if (std::isnan(islc.reco.pfp[ipfp].trk.start.y) || std::isnan(islc.reco.pfp[ipfp].trk.start.z)|| std::isnan(islc.reco.pfp[ipfp].trk.end.y) || std::isnan(islc.reco.pfp[ipfp].trk.end.z) )return 9;
    
      //skip low energy tagged pions
      TVector3 Start_mom_v;
      if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);}
      if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(139.570,2)+pow(Start_mom_v.Mag()*1000,2))-139.570)>=25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 2;}}
      
      //skip low energy protons
      if(islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 ){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);}
      if(islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(938.3,2)+pow(Start_mom_v.Mag()*1000,2))-938.3)>=50.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 1;}}
            
    }
    if(islc.reco.pfp[ipfp].trackScore<0.5)
    {
      if(islc.reco.pfp[ipfp].trackScore>=0.3 && islc.reco.pfp[ipfp].trackScore<0.5 && islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 )
      {
        TVector3 Start_mom_v2;
        Start_mom_v2.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
        if((sqrt(pow(938.3,2)+pow(Start_mom_v2.Mag()*1000,2))-938.3)>=50.0 && ((RecoVtx-start).Mag()<dist_cut) && islc.reco.pfp[ipfp].parent_is_primary){return 1;}
      }
      if(!(islc.reco.pfp[ipfp].trackScore>=0.3 && islc.reco.pfp[ipfp].trackScore<0.5 && islc.reco.pfp[ipfp].trk.truth.p.pdg==2212 )){

      if(std::isnan(islc.reco.pfp[ipfp].shw.plane[2].energy))return 9;
      if(islc.reco.pfp[ipfp].shw.plane[2].energy*1000<25.0)return 9;
      if(islc.reco.pfp[ipfp].shw.plane[2].energy*1000>25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 3;}}} 
    }
    
    return 9;
}


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

int true_selection(const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
    //-1 : unclassified
    //0 : muon rising
    //1 : muon mip
    //2 : proton rising
    //3 : proton interacting
    //4 : pion rising
    //5 : pion interacting 

    int bestplane = 2;
    bool hasValidHits=false;

    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()>0)
    {
      //check if it has valid hits for for computing the likelihood
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
      {
        if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
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
        if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
        {
          hasValidHits=true;
        }
      }

      if(hasValidHits==false)return -1;
    }
    

    //controllo che vertice reco e true siano abbastanza vicini per eliminare la componente di cosmici
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    if((vertex_true-vertex_reco).Mag()>100.)return -1;

    //compute distance between the endpoints
    TVector3 end_true;
    end_true.SetXYZ(islc.reco.pfp[ipfp].trk.truth.p.end.x, islc.reco.pfp[ipfp].trk.truth.p.end.y, islc.reco.pfp[ipfp].trk.truth.p.end.z);
    TVector3 end_hit;
    int nhits = (int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size();
    double endx = islc.reco.pfp[ipfp].trk.calo[bestplane].points[nhits-1].x;
    double endy = islc.reco.pfp[ipfp].trk.calo[bestplane].points[nhits-1].y;
    double endz = islc.reco.pfp[ipfp].trk.calo[bestplane].points[nhits-1].z;
    end_hit.SetXYZ(endx,endy,endz);
    double end_distance = (end_hit-end_true).Mag();

    //looking at the energy match
    double daughter_electron_energy_match = -1;
    double daughter_proton_energy_match = -1;
    double daughter_pion_energy_match = -1;

    for (const auto& true_p : sr->true_particles)
    {
      for (auto const& match: islc.reco.pfp[ipfp].trk.truth.matches)
      {
        if(true_p.G4ID==match.G4ID)
        {
          if(std::abs(true_p.pdg)==11 && (int)true_p.parent == islc.reco.pfp[ipfp].trk.truth.p.G4ID)
          {
            if(match.energy/3. > daughter_electron_energy_match)daughter_electron_energy_match=match.energy/3.;
          }
          if(std::abs(true_p.pdg)==2212 && (int)true_p.parent == islc.reco.pfp[ipfp].trk.truth.p.G4ID)
          {
            if(match.energy/3. > daughter_proton_energy_match)daughter_proton_energy_match=match.energy/3.;
          }
          if(std::abs(true_p.pdg)==211 && (int)true_p.parent == islc.reco.pfp[ipfp].trk.truth.p.G4ID)
          {
            if(match.energy/3. > daughter_pion_energy_match)daughter_pion_energy_match=match.energy/3.;
          }
        }
      }
    } 

    int true_class=-1;
    //TRUE CLASSIFICATION
    if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==13)
    {
        if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness>=0.5 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity>=0.5)
        {
            if(daughter_electron_energy_match<=0.01)
            {
                if(end_distance<=3.)true_class = 0;       
                else true_class = 1;
            }
            else true_class=6;
        } 
        else true_class=-1;
    }
    if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==2212)
    {
        if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness>=0.3 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity>=0.3)
        {
            if
            ( 
                (islc.reco.pfp[ipfp].trk.truth.p.end_process==54 && end_distance <= 5.) ||
                (islc.reco.pfp[ipfp].trk.truth.p.end_process!=54 && (daughter_proton_energy_match>0.055 || daughter_pion_energy_match>0.055))
            ) true_class = 2;
            else true_class = 3;
        }
        else true_class=-1;
    }
    if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211)
    {
        if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness>=0.3 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity>=0.3)
        {
            if(( islc.reco.pfp[ipfp].trk.truth.p.end_process==3 || islc.reco.pfp[ipfp].trk.truth.p.end_process==45 ) && end_distance <= 1.5) true_class = 4;
            else true_class = 5;
        }
        else true_class=-1;
    }

    return true_class;
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
      if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx<30.)
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
        if(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx<30.)
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

void write(std::ofstream& file, int clas, std::vector<double> &vec, double val, double &val1, double &val2, int &val3)
{
    file << clas << " ";
    for(int i=0; i<(int)vec.size(); i++)
    {
        file << vec[i] << " ";
    }
    file << val << " " << val1 << " " << val2 << " " << val3 << endl;
}

void write_dedx(std::ofstream& file, int clas, std::vector<double> &rr, std::vector<double> &dedx, std::vector<double> &pitch, double &dirx, double &diry, double &dirz)
{
  std::string string_class;
  if(clas==0)string_class="muon_class0";
  if(clas==1)string_class="muon_class1";
  if(clas==2)string_class="proton_class2";
  if(clas==3)string_class="proton_class3";
  if(clas==4)string_class="pion_class4";
  if(clas==5)string_class="pion_class5";
  for(int i=0; i<(int)rr.size(); i++)
  {
    file << string_class << " " << rr[i] << " " << dedx[i] << " " << pitch[i] << " " << dirx << " " << diry << " " << dirz << endl;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////                                                 Likelihood                                                   //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//ofstream dump_lkl_ratios_gbdt_withdepE("/exp/icarus/data/users/nsommagg/dump_lkl_ratios_70percet_2di2_multiplane.txt");
ofstream dump_lkl_ratios_gbdt_withdepE("likelihoods_for_training.txt");

std::array<std::vector<TH1D*>,6> prob_d_coll;
std::array<std::vector<TH1D*>,6> prob_d_ind1;
std::array<std::vector<TH1D*>,6> prob_d_ind2;

//std::array<std::vector<TH1D*>,6> prob_d_coll_150;
//std::array<std::vector<TH1D*>,6> prob_d_ind1_150;
//std::array<std::vector<TH1D*>,6> prob_d_ind2_150;

int n_c0 = 0;
int n_c1 = 0;
int n_c2 = 0;
int n_c3 = 0;
int n_c4 = 0;
int n_c5 = 0;

int n_mu = 0;
int n_pro = 0;
int n_pi = 0;

std::vector<double> genp_mu;
std::vector<double> genp_pro;
std::vector<double> genp_pi;

const SpillMultiVar likelihood_dump([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  std::vector<double> vector_active;

  for(auto const &islc : sr->slc)
  {

    //cout << "RUN: " << sr->hdr.run << " EVT: " << sr->hdr.evt << endl; 
    //cout << "slice" << " vtxx: " << islc.truth.position.x << " vtxy: " << islc.truth.position.y << " vtxz: " << islc.truth.position.z << endl;
    //cout << "primary particles: ";
    //for(auto const &itrue : islc.truth.prim)
    //{
    //  cout << itrue.pdg << " ";
    //}
    //cout << endl;

    //cout << " in this slice: ";
    //for (std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp){cout << islc.reco.pfp[ipfp].trk.truth.p.pdg << " ";}
    //cout << endl;
 
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    
    if (islc.tmatch.eff <= 0.5 || (vertex_true-vertex_reco).Mag() >= 100 || islc.is_clear_cosmic)continue;

    int ipfp_mu=-1;
    ipfp_mu=find_truth_muon(islc,10);
    
    for (std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
    {   

        int bestplane = 2;
        bool hasValidHits=false;

        if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()>0)
        {
          //check if it has valid hits for for computing the likelihood
          for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
          {
            if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
            {
              hasValidHits=true;
            }
          }   
        }

        if(!hasValidHits)continue;
    
        //--> looking at hits from ind1 or ind2 if no valid hits in collection
        /*
        if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()==0 || hasValidHits==false)
        {
          bestplane = islc.reco.pfp[ipfp].trk.bestplane;
          if(bestplane==-1)continue;

          //check if it has valid hits for for computing the likelihood
          for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
          {
            if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
            {
              hasValidHits=true;
            }
          }

          if(hasValidHits==false)continue;
        }
        */

        // ---->
        
        std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
        double dep_E = 0;
        double dep_E_length_scaled = 0;
        std::vector<double> rrs_last5;
        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit )
        {
            if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr <= 5.)
            {   
                rrs_last5.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr);
                dep_E = dep_E + islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch;
            }
        }

        std::vector<double> rr;
        std::vector<double> dedx;
        std::vector<double> pitch;

        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit )
        {
          rr.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr);
          dedx.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx);
          pitch.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch);
        }


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

        if(daughter_id==-1)
        {
          daughter_depE = -1;
          angle_end = -1;
        }

        // <----

        int pfp_id = id_pfp_truth(islc, ipfp, 10);

        if(int(ipfp)==ipfp_mu)
        {
            n_mu ++;
            if(true_selection(sr,islc,ipfp)==0)
            { 
              //cout << "muon class 0 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
              n_c0 ++;
              write(dump_lkl_ratios_gbdt_withdepE,0,lr,dep_E,daughter_depE,angle_end,bestplane);
              //write_dedx(dump_dedx_range,0,rr,dedx,pitch,dirx,diry,dirz);
              //dump_depE << "muon_class0" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_pitch_scaled << " " << max_rr << " " << pitch_sum << endl; 
              //cout << "muon_class0" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_hit_scaled << endl; 
            }
            if(true_selection(sr,islc,ipfp)==1)
            {
              //cout << "muon class 1 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
              n_c1 ++;
              write(dump_lkl_ratios_gbdt_withdepE,1,lr,dep_E,daughter_depE,angle_end,bestplane);
              //write_dedx(dump_dedx_range,1,rr,dedx,pitch,dirx,diry,dirz);
              //dump_depE << "muon_class1" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_pitch_scaled << " " << max_rr << " " << pitch_sum << endl; 
              //cout << "muon_class1" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_hit_scaled << endl; 
            }
            //if(true_selection(sr,islc,ipfp)==-1)
            //{
            //  cout << "muon class -1 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
            //}
            continue;
        }
        else if(pfp_id==1)
        {
            n_pro ++;
            if(true_selection(sr,islc,ipfp)==2)
            {
              //cout << "proton class 2 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
              n_c2 ++;
              write(dump_lkl_ratios_gbdt_withdepE,2,lr,dep_E,daughter_depE,angle_end,bestplane);
              //write_dedx(dump_dedx_range,2,rr,dedx,pitch,dirx,diry,dirz);
              //dump_depE << "proton_class2" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_pitch_scaled << " " << max_rr << " " << pitch_sum << endl; 
              //cout << "proton_class2" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_hit_scaled << endl;
            }
            if(true_selection(sr,islc,ipfp)==3)
            {
              //cout << "proton class 3 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
              n_c3 ++;
              write(dump_lkl_ratios_gbdt_withdepE,3,lr,dep_E,daughter_depE,angle_end,bestplane);
              //write_dedx(dump_dedx_range,3,rr,dedx,pitch,dirx,diry,dirz);
              //dump_depE << "proton_class3" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_pitch_scaled << " " << max_rr << " " << pitch_sum << endl; 
              //cout << "proton_class3" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_hit_scaled << endl; 
            }
            //if(true_selection(sr,islc,ipfp)==-1)
            //{
            //  cout << "proton class -1 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
            //}
        }
        else if(pfp_id==2)
        {
            n_pi ++;
            if(true_selection(sr,islc,ipfp)==4)
            {
              //cout << "pion class 4 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
              n_c4 ++;
              write(dump_lkl_ratios_gbdt_withdepE,4,lr,dep_E,daughter_depE,angle_end,bestplane);
              //write_dedx(dump_dedx_range,4,rr,dedx,pitch,dirx,diry,dirz);
              //dump_depE << "pion_class4" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_pitch_scaled << " " << max_rr << " " << pitch_sum << endl; 
              //cout << "pion_class4" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_hit_scaled << endl; 
            }
            if(true_selection(sr,islc,ipfp)==5)
            {
              //cout << "pion class 5 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
              n_c5 ++;
              write(dump_lkl_ratios_gbdt_withdepE,5,lr,dep_E,daughter_depE,angle_end,bestplane);
              //write_dedx(dump_dedx_range,5,rr,dedx,pitch,dirx,diry,dirz);
              //dump_depE << "pion_class5" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_pitch_scaled << " " << max_rr << " " << pitch_sum << endl; 
              //cout << "pion_class5" << " " << dep_E << " " << dep_E_length_scaled << " " << dep_E_hit_scaled << endl;
            }
            //if(true_selection(sr,islc,ipfp)==-1)
            //{
            //  cout << "pion class -1 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
            //}
        }
        //else
        //{
        //  cout << "other pfp " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
        //}
    }
    //cout << endl;
  }

  return vector_active;
});
