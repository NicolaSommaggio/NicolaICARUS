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

TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
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


bool isInContained (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
  return (( ( x < -61.94 - 5 && x > -358.49 + 5 ) ||
			( x >  61.94 + 5 && x <  358.49 - 5 )) &&
		  ( ( y > -181.86 + 5 && y < 134.96 - 5 ) &&
		  ( z > -894.95 + 5 && z < 894.95 - 5 ) ));

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



int find_muon ( const caf::Proxy<caf::SRSlice>& islc, int dist_mucut) { 

    //Select muon as longest track
    double max_length=-1.0;
    int ipfp_mu = -1;
    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 RecoStart;
    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
    {
        if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;
        RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
        if(islc.reco.pfp[ipfp].trackScore<0.5)continue;
        int use_plane = 2;

        //compute new chi2
        std::vector<double> output;
        std::vector<double> dedx;
        std::vector<double> rr;
        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit )
        {
            dedx.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
            rr.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);                    
        } // calo points

        //input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
        //output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt} 
        output = chi2_ALG(dedx,rr,0.0,25.0);

        if( islc.reco.pfp[ipfp].trk.len>max_length && 
            ((RecoVtx-RecoStart).Mag()<dist_mucut) && 
            islc.reco.pfp[ipfp].trk.len>50 && 
            output[0]<30 && output[1]>60 && 
            isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) && 
            (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && 
            islc.reco.pfp[ipfp].parent_is_primary
        )
        {
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

      int use_plane2=plane;
      if(std::isnan(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy))return 9;
      if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000<25.0)return 9;
      if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000>25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 3;}}} 
    }
    
    return 9;
}

int find_truth_muon ( const caf::Proxy<caf::SRSlice>& islc, int dist_mucut) { 

    //Select muon as longest track
    double max_length=-1.0;
    int ipfp_mu = -1;
    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 RecoStart;

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){

        if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) continue;

        RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);

        if(islc.reco.pfp[ipfp].trackScore<0.5)continue;

        if(islc.reco.pfp[ipfp].trk.len>max_length && ((RecoVtx-RecoStart).Mag()<dist_mucut) && islc.reco.pfp[ipfp].trk.len>50 && 
        std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==13 && 
        isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) && 
        (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && islc.reco.pfp[ipfp].parent_is_primary){
        max_length=islc.reco.pfp[ipfp].trk.len;
        ipfp_mu=ipfp;
            }
        }//loop of pfp to find muon
return ipfp_mu;
}


int id_pfp ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_cut ) { 
    //return 1 PROTONS
    //return 2 PIONS
    //return 3 SHOWER
    //return 9 other -> nan, not primary, too far, below energy threshold... 


    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);

    //skip secondaries
    if (!(islc.reco.pfp[ipfp].parent_is_primary ))return 9;
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;

    //consider only primary tracks which are 20cm close to the vertex, either vtx-start or vtx-end
    TVector3 start(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
    TVector3 end(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);
    //condition ? result_if_true : result_if_false
    double min_dist = ((start-RecoVtx).Mag()< (end-RecoVtx).Mag() ? (start-RecoVtx).Mag() : (end-RecoVtx).Mag());
    if(min_dist>50.0) return 9;
     
    int use_plane = 2;
    //compute new chi2
    std::vector<double> output;
    std::vector<double> dedx;
    std::vector<double> rr;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit )
    {
        dedx.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
        rr.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
    } // calo points
    //input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
    //output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt} 
    output = chi2_ALG(dedx,rr,0.0,25.0);
    if(islc.reco.pfp[ipfp].trackScore>=0.5)
    {
        if (std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
        if (std::isnan(islc.reco.pfp[ipfp].trk.start.y) || std::isnan(islc.reco.pfp[ipfp].trk.start.z)|| std::isnan(islc.reco.pfp[ipfp].trk.end.y) || std::isnan(islc.reco.pfp[ipfp].trk.end.z) )return 9;
    
        //skip low energy tagged pions
        TVector3 Start_mom_v;
        if(output[1]>=100 ){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);}
        if(output[1]>=100 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(139.570,2)+pow(Start_mom_v.Mag()*1000,2))-139.570)>=25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 2;}}
        //skip low energy protons
        if(output[1]<100 ){ Start_mom_v.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);}
        if(output[1]<100 && ((RecoVtx-start).Mag()<dist_cut) && (sqrt(pow(938.3,2)+pow(Start_mom_v.Mag()*1000,2))-938.3)>=50.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 1;}}       
    }
    if(islc.reco.pfp[ipfp].trackScore<0.5)
    {   
        if(islc.reco.pfp[ipfp].trackScore>=0.3 && islc.reco.pfp[ipfp].trackScore<0.5 && output[1]<100 )
        {
            TVector3 Start_mom_v2;
            Start_mom_v2.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
            if((sqrt(pow(938.3,2)+pow(Start_mom_v2.Mag()*1000,2))-938.3)>=50.0 && ((RecoVtx-start).Mag()<dist_cut) && islc.reco.pfp[ipfp].parent_is_primary){return 1;}
        }
        
        if(!(islc.reco.pfp[ipfp].trackScore>=0.3 && islc.reco.pfp[ipfp].trackScore<0.5 && output[1]<100 ))
        {
            //int use_plane2 = islc.reco.pfp[ipfp].trk.calo[2].nhit>islc.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;  
            int use_plane2 = 2;
            if(std::isnan(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy))return 9;
            if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000<25.0)return 9;
            if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000>25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 3;}}
        }
        
    }
    
    return 9;
}


int classification_type ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc ) { 

    // return 2 -> Visible and true 1mu1p
    // return 5 -> Visible and true 1muNp
    // return 3 -> Other 

            //CHECK WHICH TYPE OF INTERACTION 
    if(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)) return 99;
    if(std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)) return 99;
                if ( abs(islc.truth.pdg) == 14 && islc.truth.iscc && !std::isnan(islc.truth.position.x) && !std::isnan(islc.truth.position.y) && !std::isnan(islc.truth.position.z) &&
		        isInActive(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z) ){          
                if(isInFV(islc.truth.position.x,islc.truth.position.y,islc.truth.position.z)){

                 int num_protons_above50 = 0;
                 int num_muons = 0; 
                 double length_muon = 0;                 
                double dep_E=0;
                for ( auto const& ipart : islc.truth.prim ){
                    if ( ipart.G4ID < 0 )  continue;

                    if(abs(ipart.pdg) == 13){num_muons+=1; length_muon =ipart.length;}  // muons
                    if(abs(ipart.pdg)==211){return 3;}  //Charged pions
                    int iG4ID_parent;  
                    //int use_plane = ipart.plane[ipart.cryostat][2].nhit>ipart.plane[ipart.cryostat][1].nhit ? 2:1;
                    int use_plane = 2;
                    
                    if(abs(ipart.pdg)==111){            //Neutral pions - reject if any of its gamma > 25 MeV
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue2 : sr->true_particles ){
                        iG4ID_parent=itrue2.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID && abs(itrue2.pdg) == 22){if(itrue2.plane[ipart.cryostat][use_plane].visE*1000>25){return 3;};
                        }}}                        
                        }
                    if(abs(ipart.pdg) == 22){                    
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }}}
                        dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                        } 
                    if(abs(ipart.pdg)== 22 && dep_E>25.0){return 3;}   //primary photons   

                    if(abs(ipart.pdg)== 2212){                    
                    if(ipart.daughters.size()>0){
                        for ( auto const& itrue : sr->true_particles ){
                        iG4ID_parent=itrue.parent;
                        //sum depE daughters 
                        if(iG4ID_parent==ipart.G4ID ){dep_E+=itrue.plane[ipart.cryostat][use_plane].visE*1000;
                        }}}
                        dep_E += ipart.plane[ipart.cryostat][use_plane].visE*1000;
                        }
                    if(abs(ipart.pdg)== 2212 && dep_E>50.0){num_protons_above50+=1;} //protons
                    dep_E=0;  
                    }
                if(num_muons == 1 && num_protons_above50 == 1){
                    if(length_muon > 50. && all_contained_truth(sr, islc)){return 2; }
                }
                if(num_muons == 1 && num_protons_above50 > 1){
                    if(length_muon > 50. && all_contained_truth(sr, islc)){return 5; }
                }
             }//fiducial 
             }//active
         
return 3;
}

bool automatic_selection_1muNp ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_cut, int cut_baryc ){

    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z))
    {
        int ipfp_mu = -1;
        if(islc.is_clear_cosmic==0)
        {
            if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
            {
                if(all_contained(islc))
                {
                    ipfp_mu=find_muon(islc,dist_cut);
                    if(ipfp_mu!=-1)
                    {
                        int num_protons =0;
                        int num_pions =0;
                        int num_showers =0; 

                        for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
                        {
                            if(int(ipfp)==ipfp_mu)continue;
                            if(id_pfp(islc, ipfp,dist_cut)==1){num_protons+=1;}
                            if(id_pfp(islc, ipfp,dist_cut)==2){num_pions+=1;}
                            if(id_pfp(islc, ipfp,dist_cut)==3){num_showers+=1;}
                        }
                        if(num_protons>0 && num_pions==0 && num_showers==0)
                        {
                            return true;
                        }//1mu1p 
                    }//muon with conditions found
                }//all tracks of slice contained  
            }//new Barycenter match
        }//fiducial condition
    }//check no nan vertex
   
return false;
}


struct SliceDrawInfo
{
    int ipfp_mu;
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

SliceDrawInfo GetSliceDrawInfo(const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, std::vector<int> v_ipfp_pro, std::vector<int> v_ipfp_pi, std::vector<std::pair<int,std::array<double,6>>> true_p)
{
    SliceDrawInfo thislice;

    thislice.true_particles=true_p;
    std::vector<std::vector<std::array<double,3>>> particles;
    std::vector<std::vector<double>> particles_dedx;
    thislice.vertex[0] = islc.vertex.x;
    thislice.vertex[1] = islc.vertex.y;
    thislice.vertex[2] = islc.vertex.z;

    thislice.ipfp_mu=ipfp_mu;
    thislice.v_ipfp_pro= v_ipfp_pro;
    thislice.v_ipfp_pi = v_ipfp_pi;

    for(std::size_t ipfp(0); ipfp<islc.reco.npfp; ipfp++)
    {
        //int plane=-1;
        //cout << islc.reco.pfp[ipfp].trk.calo[2].points.size() << " " << islc.reco.pfp[ipfp].trk.calo[0].points.size() << " " << islc.reco.pfp[ipfp].trk.calo[1].points.size() << endl;
        //if(int(islc.reco.pfp[ipfp].trk.calo[2].points.size())==0) plane=islc.reco.pfp[ipfp].trk.bestplane;
        //else plane=2;
        //int plane=islc.reco.pfp[ipfp].trk.bestplane;
        //if(plane<0 || plane>3) continue;

        int plane=2;

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
        if(abs(prim.pdg)==2212 && dep_E<50.)continue;

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

std::vector<unsigned int> run = {};
std::vector<unsigned int> subrun = {};
std::vector<unsigned int> evt = {};

ofstream dumpInfo("logSelection.txt");
int slice_counter=0;

std::vector<std::vector<std::pair<double,double>>> dumpRM;
std::vector<double> dumpMedian;
std::vector<std::vector<double>> dedx;
std::vector<std::vector<double>> rr;
std::vector<std::vector<double>> x;

int conta_muoni=0;

const SpillMultiVar disegna([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    //bool is_true_signal=false;

    for (auto const& islc : sr->slc)
    {     

        if(sr->hdr.evt % 2 !=0)continue;

        bool is_true_signal=false;
        int ipfp_mu=-1;
        std::vector<int> v_ipfp_pro;
        std::vector<int> v_ipfp_pi;


        //if(islc.truth.index>=0 && (classification_type(sr,islc)==2 || classification_type(sr,islc)==5 ) && (isInFV(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z)) ){is_true_signal=true;}

        //EVENT SELECTION
        
        
        //if (automatic_selection_1muNp(sr,islc,10,100))//1mu1p reco
	    //{           

            //if(!is_true_signal)continue;

            slice_counter+=1;
            //IDENTIFY MUONS AND PROTONS WITH CHI2
            ipfp_mu=find_truth_muon(islc,10);

            if(ipfp_mu!=-1)conta_muoni++;

	        for(std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
		    {
                if(int(ipfp)==ipfp_mu)continue;
		        if(id_pfp_truth(islc, ipfp, 10,2)==1){v_ipfp_pro.push_back(int(ipfp));}
                if(id_pfp_truth(islc, ipfp, 10,2)==2){v_ipfp_pi.push_back(int(ipfp));}
            }
            //if(v_ipfp_pi.size()==0)continue;

            //rolling median
            std::vector<std::pair<double,double>> rm_thislice={{-1,-1}};
            if(ipfp_mu!=-1)
            {
                dumpRM.push_back(rollingMedian(islc,ipfp_mu));
                rm_thislice=rollingMedian(islc,ipfp_mu);
            }


            //median
            double median=-1;
            std::vector<double> dummy;
            std::vector<double> dedx_temp;
            std::vector<double> rr_temp;
            std::vector<double> x_temp;
            if(ipfp_mu!=-1)
            {
                for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit) 
                {
                    if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr<30)
                    {
                        dedx_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx);
                        rr_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr);
                        x_temp.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x);
                    }
                    if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr<5.) 
                    {
                        dummy.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx);
                    }
                }
                if(dummy.size()!=0) {dumpMedian.push_back(mediana(dummy)); median=mediana(dummy);}
                else {dumpMedian.push_back(-1);}
                dedx.push_back(dedx_temp);
                rr.push_back(rr_temp);
                x.push_back(x_temp);
            }

            
            int max_length_pfp=-1;
            double max_length=-1;

            
            dumpInfo << "*** RUN " << sr->hdr.run << " *** EVT " << sr->hdr.evt << " ***" << endl;

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
                dumpInfo << "is contained " << isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) << endl;
                dumpInfo << "is the same TPC as vertex " << ((islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0) << endl;
                dumpInfo << "parent is primary " << islc.reco.pfp[ipfp].parent_is_primary << endl;
                Start_mom_v_proton.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
                Start_mom_v_pion.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);
                dumpInfo << "energy lost as proton " << sqrt(pow(938.3,2)+pow(Start_mom_v_proton.Mag()*1000,2))-938.3 << " energy lost as pion " << sqrt(pow(139.570,2)+pow(Start_mom_v_pion.Mag()*1000,2))-139.570 << endl;
     
                if(int(ipfp)==ipfp_mu)
                {
                    dumpInfo << "median " << median << endl;
                    dumpInfo << "rolling median ";
                    for(int i=0; i<int(rm_thislice.size()); i++) dumpInfo << rm_thislice[i].first << " ";
                    dumpInfo << endl;
                    for(int i=0; i<int(rm_thislice.size()); i++) dumpInfo << rm_thislice[i].second << " ";
                    dumpInfo << endl;
                }
                
    
                if(int(ipfp)==ipfp_mu) dumpInfo << "MUON" << endl << endl;
                else if(id_pfp_truth(islc, ipfp, 10,2)==1) dumpInfo << "PROTON" << endl << endl;
                else if(id_pfp_truth(islc, ipfp, 10,2)==2) dumpInfo << "PION" << endl << endl;
                else dumpInfo << "OTHER" << endl << endl;

            }

            dumpInfo << endl;
            dumpInfo << "*********************************** end slice ***********************************" << endl;
            
            /*
            bool match_found = false;
            
            for (int i = 0; i < int(run.size()); i++) 
            {
                if (sr->hdr.run == run[i] &&
                sr->hdr.subrun == subrun[i] &&
                sr->hdr.evt == evt[i]) 
                {
                    match_found = true; 
                    break; // trovato, non serve continuare a cercare
                }
            }
            if(!match_found)continue;
            */

            //cout << sr->hdr.run << " " << sr->hdr.subrun << " " << sr->hdr.evt << endl;
            //cout << "muone " << islc.reco.pfp[ipfp_mu].trk.len << " " << islc.reco.pfp[ipfp_mu].trk.start.x << endl;
            //for(auto const & ipfp_pro : v_ipfp_pro)
            //{
                //cout << "protone " << islc.reco.pfp[ipfp_pro].trk.len << " " << islc.reco.pfp[ipfp_pro].trk.start.x << endl;
            //}
            //for(std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
            //{
                //if(int(ipfp)==ipfp_mu)continue;
                //bool is_a_proton=false;
                //for(auto const &ipfp_pro : v_ipfp_pro){if(int(ipfp)==ipfp_pro){is_a_proton=true;}}
                //if(!is_a_proton)
                //{
                    //cout << "other " << islc.reco.pfp[ipfp].trk.len << " " << islc.reco.pfp[ipfp].trk.start.x << endl;
                //}
            //}

            //if(median > 1 && median < 3.2) {

            //INFORMAZIONI PER DISEGNARE LA SLICE 3D
            std::vector<double> nu_int_temp;
            nu_int_temp.push_back(islc.truth.E);
            nu_int_temp.push_back(islc.truth.targetPDG);
            nu_int_temp.push_back(islc.truth.Q2);
            std::vector<std::pair<int,std::array<double,6>>> true_particles=GetTrueParticlesInfo(sr,islc);
            SliceDrawInfo drawInfo = GetSliceDrawInfo(islc,ipfp_mu,v_ipfp_pro,v_ipfp_pi,true_particles);
            drawInfo.neutrino_interaction=nu_int_temp;

            //INFORMAZIONI PER DISEGNARE LA SLICE 2D
            //std::vector<std::pair<int,std::array<double,8>>> true_particles2D=GetTrueParticlesInfo2D(sr,islc);
            //SliceDrawInfo2D drawInfo2D = GetSliceDrawInfo2D(islc,ipfp_mu,v_ipfp_pro,v_ipfp_pi/*,true_particles2D*/);
            //drawInfo2D.neutrino_interaction=nu_int_temp;

            slices3D.push_back(drawInfo);
            //slices2D.push_back(drawInfo2D);
            //}

        //}//STANDARD SELECTION

    }

    return vector_active;
});
