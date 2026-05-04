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


int find_truth_muon ( const caf::Proxy<caf::SRSlice>& islc, int dist_mucut, int plane) { 

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


int id_pfp_truth ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_cut, int plane ) { 
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
     
    int use_plane = plane;
    
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

      int use_plane2=use_plane;
      if(std::isnan(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy))return 9;
      if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000<25.0)return 9;
      if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000>25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 3;}}} 
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

int true_selection(const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp, int plane)
{
    //-1 : unclassified
    //0 : muon rising
    //1 : muon mip
    //2 : proton rising
    //3 : proton interacting
    //4 : pion rising
    //5 : pion interacting 

    if((int)islc.reco.pfp[ipfp].trk.calo[plane].points.size()==0)return -1;

    //check if it has valid hits for for computing the likelihood
    bool hasValidHits=false;
    for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[plane].points.size(); ++ihit)
    {
      if(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx<30.)
      {
        hasValidHits=true;
      }
    }
    if(hasValidHits==false)return -1;

    //controllo che vertice reco e true siano abbastanza vicini per eliminare la componente di cosmici
    //TVector3 vertex_true;
    //vertex_true.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    //TVector3 vertex_reco;
    //vertex_reco.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    //if((vertex_true-vertex_reco).Mag()>100.)return -1;

    //compute distance between the endpoints
    TVector3 end_true;
    end_true.SetXYZ(islc.reco.pfp[ipfp].trk.truth.p.end.x, islc.reco.pfp[ipfp].trk.truth.p.end.y, islc.reco.pfp[ipfp].trk.truth.p.end.z);
    TVector3 end_hit;
    int nhits = (int)islc.reco.pfp[ipfp].trk.calo[plane].points.size();
    double endx = islc.reco.pfp[ipfp].trk.calo[plane].points[nhits-1].x;
    double endy = islc.reco.pfp[ipfp].trk.calo[plane].points[nhits-1].y;
    double endz = islc.reco.pfp[ipfp].trk.calo[plane].points[nhits-1].z;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////                              seleclet true class to make prob densities                                      //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> dedx_mu_class0;
std::vector<std::vector<double>> rr_mu_class0;
std::vector<std::vector<double>> dedx_mu_class1;
std::vector<std::vector<double>> rr_mu_class1;
std::vector<std::vector<double>> dedx_pro_class2;
std::vector<std::vector<double>> rr_pro_class2;
std::vector<std::vector<double>> dedx_pro_class3;
std::vector<std::vector<double>> rr_pro_class3;
std::vector<std::vector<double>> dedx_pi_class4;
std::vector<std::vector<double>> rr_pi_class4;
std::vector<std::vector<double>> dedx_pi_class5;
std::vector<std::vector<double>> rr_pi_class5;

std::vector<double> rr_temp;
std::vector<double> dedx_temp;

int plane;

const SpillMultiVar select_true_class([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  std::vector<double> vector_active;

  for(auto const &islc : sr->slc)
  {
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    
    if (islc.tmatch.eff <= 0.5 || (vertex_true-vertex_reco).Mag() >= 100 || islc.is_clear_cosmic)continue;

    int ipfp_mu=-1;
    ipfp_mu=find_truth_muon(islc,10,2);
    if(ipfp_mu!=-1)
    {
        
        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[plane].points.size(); ++ihit )
        {
            if(islc.reco.pfp[ipfp_mu].trk.calo[plane].points[ihit].rr>1 && islc.reco.pfp[ipfp_mu].trk.calo[plane].points[ihit].rr<25 && islc.reco.pfp[ipfp_mu].trk.calo[plane].points[ihit].dedx>1 && islc.reco.pfp[ipfp_mu].trk.calo[plane].points[ihit].dedx<30.)
            {    
                dedx_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[plane].points[ihit].dedx);                                  
                rr_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[plane].points[ihit].rr);
            }
        }
        if(true_selection(sr,islc,ipfp_mu,plane)==0){dedx_mu_class0.push_back(dedx_temp); rr_mu_class0.push_back(rr_temp);}
        if(true_selection(sr,islc,ipfp_mu,plane)==1){dedx_mu_class1.push_back(dedx_temp); rr_mu_class1.push_back(rr_temp);}
        dedx_temp.clear(); 
        rr_temp.clear();
    }
    
    for (std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
    {
        if(int(ipfp)==ipfp_mu)continue;

        if(id_pfp_truth(islc, ipfp, 10,2)==1)
        {
            for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[plane].points.size(); ++ihit )
            {
                if(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr>1 && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx>1 && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx<30.)
                {    
                    dedx_temp.push_back(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx);                                  
                    rr_temp.push_back(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr);
                }
            }
            if(true_selection(sr,islc,ipfp,plane)==2){dedx_pro_class2.push_back(dedx_temp); rr_pro_class2.push_back(rr_temp);}
            if(true_selection(sr,islc,ipfp,plane)==3){dedx_pro_class3.push_back(dedx_temp); rr_pro_class3.push_back(rr_temp);}
            dedx_temp.clear(); 
            rr_temp.clear();
        }

        if(id_pfp_truth(islc, ipfp, 10,2)==2)
        {
            for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[plane].points.size(); ++ihit )
            {
                if(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr>1 && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx>1 && islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx<30.)
                {    
                    dedx_temp.push_back(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx);                                  
                    rr_temp.push_back(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr);
                }
            }
            if(true_selection(sr,islc,ipfp,plane)==4){dedx_pi_class4.push_back(dedx_temp); rr_pi_class4.push_back(rr_temp);}
            if(true_selection(sr,islc,ipfp,plane)==5){dedx_pi_class5.push_back(dedx_temp); rr_pi_class5.push_back(rr_temp);}
            dedx_temp.clear(); 
            rr_temp.clear();
        }
    }
  }

  return vector_active;
});
