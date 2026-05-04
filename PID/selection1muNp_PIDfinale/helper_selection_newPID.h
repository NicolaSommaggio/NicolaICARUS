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


int find_muon ( const caf::Proxy<caf::SRSlice>& islc, int dist_mucut, int plane) { 

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
   
  //int use_plane = islc.reco.pfp[ipfp].trk.calo[2].nhit>islc.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
    //int use_plane = plane;
    int use_plane = -1;
    if(plane==-1)
    {
      if(islc.reco.pfp[ipfp].trk.bestplane==0 || islc.reco.pfp[ipfp].trk.bestplane==1 || islc.reco.pfp[ipfp].trk.bestplane==2){use_plane=islc.reco.pfp[ipfp].trk.bestplane;}
      else{use_plane=2;}
    }
    else{use_plane=plane;}
    //compute new chi2
    std::vector<double> output;
    std::vector<double> dedx;
    std::vector<double> rr;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit ){
                        dedx.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
                        rr.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
                            } // calo points
                  //input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
    //output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt} 
    output = chi2_ALG(dedx,rr,0.0,25.0);


        if(islc.reco.pfp[ipfp].trk.len>max_length && ((RecoVtx-RecoStart).Mag()<dist_mucut) && islc.reco.pfp[ipfp].trk.len>50 && 
        output[0]<30 && output[1]>60 && 
        isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) && 
        (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && islc.reco.pfp[ipfp].parent_is_primary){
        max_length=islc.reco.pfp[ipfp].trk.len;
        ipfp_mu=ipfp;
            }
        }//loop of pfp to find muon
  return ipfp_mu;
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

int id_pfp ( const caf::Proxy<caf::SRSlice>& islc, int ipfp, int dist_cut, int plane ) { 
    //return 1 PROTONS
    //return 2 PIONS
    //return 3 SHOWER
    //return 9 other -> nan, not primary, too far, below energy threshold... 


    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    //    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
    //skip secondaries
    if (!(islc.reco.pfp[ipfp].parent_is_primary ))return 9;
    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.end.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return 9;
    //if(int(ipfp)==ipfp_mu)continue;     //There is always a muon, for a 1mu1p we need 2 tracks - 1 muon = 1 only proton with threshold
    //consider only primary tracks which are 20cm close to the vertex, either vtx-start or vtx-end
    TVector3 start(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);
    TVector3 end(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z);
    //condition ? result_if_true : result_if_false
    double min_dist = ((start-RecoVtx).Mag()< (end-RecoVtx).Mag() ? (start-RecoVtx).Mag() : (end-RecoVtx).Mag());
    //if(min_dist>50.0)continue;
    //if(min_dist>10.0) return 9;
    if(min_dist>50.0) return 9;
     

    //int use_plane = islc.reco.pfp[ipfp].trk.calo[2].nhit>islc.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
    int use_plane = plane;
    
    //compute new chi2
    std::vector<double> output;
    std::vector<double> dedx;
    std::vector<double> rr;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit ){
                        dedx.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
                        rr.push_back(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
                            } // calo points
                  //input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
    //output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt} 
    output = chi2_ALG(dedx,rr,0.0,25.0);
    if(islc.reco.pfp[ipfp].trackScore>=0.5){
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
    if(islc.reco.pfp[ipfp].trackScore<0.5){
        if(islc.reco.pfp[ipfp].trackScore>=0.3 && islc.reco.pfp[ipfp].trackScore<0.5 && output[1]<100 ){
            TVector3 Start_mom_v2;
            Start_mom_v2.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
            if((sqrt(pow(938.3,2)+pow(Start_mom_v2.Mag()*1000,2))-938.3)>=50.0 && ((RecoVtx-start).Mag()<dist_cut) && islc.reco.pfp[ipfp].parent_is_primary){return 1;}
        }
    if(!(islc.reco.pfp[ipfp].trackScore>=0.3 && islc.reco.pfp[ipfp].trackScore<0.5 && output[1]<100 )){
    //int use_plane2 = islc.reco.pfp[ipfp].trk.calo[2].nhit>islc.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;  
    //int use_plane2 = 2;
    int use_plane2=use_plane;
    if(std::isnan(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy))return 9;
    if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000<25.0)return 9;
    if(islc.reco.pfp[ipfp].shw.plane[use_plane2].energy*1000>25.0){if (islc.reco.pfp[ipfp].parent_is_primary ) {return 3;}}}
        }
    
    return 9;
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


bool automatic_selection_1muNp ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_cut, int cut_baryc, int plane ){

        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/)){
        {            

        {
        int ipfp_mu = -1;
        int ipfp_pro = -1;
    if( !(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z)){

        if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
    {
            if(all_contained(islc)){
                ipfp_mu=find_muon(islc,dist_cut,plane);
                if(ipfp_mu!=-1){
                    int num_protons =0;
                    int num_pions =0;
                    int num_showers =0; 

                    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){

                        int use_plane = -1;
                        if(plane==-1)
                        {
                          if(islc.reco.pfp[ipfp].trk.bestplane==0 || islc.reco.pfp[ipfp].trk.bestplane==1 || islc.reco.pfp[ipfp].trk.bestplane==2){use_plane=islc.reco.pfp[ipfp].trk.bestplane;}
                          else{use_plane=2;}
                        }
                        else{use_plane=plane;}

                        if(int(ipfp)==ipfp_mu)continue;
                        if(id_pfp(islc, ipfp,dist_cut, use_plane)==1){num_protons+=1;}
                        if(id_pfp(islc, ipfp,dist_cut, use_plane)==2){num_pions+=1;}
                        if(id_pfp(islc, ipfp,dist_cut, use_plane)==3){num_showers+=1;}
                    }
                    if(num_protons>0 && num_pions==0 && num_showers==0){
                        return true;
                            }//1mu1p 
                        }//muon with conditions found
                    
                    }//all tracks of slice contained  
        }//new Barycenter match
    }//fiducial condition
 
    }//signal
    }//only neutrinos in active!
    }//check no nan in true info

   
  return false;
}

bool automatic_selection_1muNp_truth ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, int dist_cut, int cut_baryc, int plane ){

        if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z)/* || std::isnan(islc.truth.position.x) || std::isnan(islc.truth.position.y) || std::isnan(islc.truth.position.z)*/))
        {
          int ipfp_mu = -1;
          int ipfp_pro = -1;
          if(!(std::isnan(islc.vertex.x) || std::isnan(islc.vertex.y) || std::isnan(islc.vertex.z) || std::isnan(islc.charge_center.z)) && isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z))
          {
            if(islc.barycenterFM.deltaZ_Trigger < cut_baryc && islc.barycenterFM.deltaZ_Trigger>0)   
            {
              if(all_contained(islc))
              {

                ipfp_mu=find_truth_muon(islc,dist_cut,plane);
                if(ipfp_mu!=-1)
                {
                    int num_protons =0;
                    int num_pions =0;
                    int num_showers =0; 

                    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
                    {
                        if(int(ipfp)==ipfp_mu)continue;
                        int use_plane = -1;
                        if(plane==-1)
                        {
                          if(islc.reco.pfp[ipfp].trk.bestplane==0 || islc.reco.pfp[ipfp].trk.bestplane==1 || islc.reco.pfp[ipfp].trk.bestplane==2){use_plane=islc.reco.pfp[ipfp].trk.bestplane;}
                          else{use_plane=2;}
                        }
                        else{use_plane=plane;}

                        if(id_pfp_truth(islc, ipfp,dist_cut, use_plane)==1){num_protons+=1;}
                        if(id_pfp_truth(islc, ipfp,dist_cut, use_plane)==2){num_pions+=1;}
                        if(id_pfp_truth(islc, ipfp,dist_cut, use_plane)==3){num_showers+=1;}
                    }
                    if(num_protons>0 && num_pions==0 && num_showers==0)
                    {
                        return true;
                    }//1mu1p 
                }//muon with conditions found
              }//all tracks of slice contained  
            }//new Barycenter match
          }//fiducial condition
        }//check no nan in true info
return false;
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

std::array<std::vector<double>,2> rollingMedian(const caf::Proxy<caf::SRSlice>& islc, int ipfp, int plane)
{
    std::vector<double> dedx_temp;
    std::vector<double> rr_temp;
    std::vector<double> temp_rr_rm={-1.};
    std::vector<double> temp_rm={-1.};
    std::array<std::vector<double>,2> rm;
    for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[plane].points.size(); ++ihit)
    {
        if(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr<30)
        {
            dedx_temp.push_back(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx);
            rr_temp.push_back(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr);
        }
    }
    if(dedx_temp.size()<9)return {temp_rr_rm,temp_rm};
    
    for(int i=4; i+4<int(dedx_temp.size()); i++)
    {
        std::vector<double> dummy = {dedx_temp[i-4], dedx_temp[i-3], dedx_temp[i-2], dedx_temp[i-1], dedx_temp[i] ,dedx_temp[i+1],dedx_temp[i+2], dedx_temp[i+3], dedx_temp[i+4]};
        double media_rr = (rr_temp[i-4] + rr_temp[i-3] + rr_temp[i-2] + rr_temp[i-1] + rr_temp[i] + rr_temp[i+1] + rr_temp[i+2] + rr_temp[i+3] + rr_temp[i+4])/9.;
        rm[0].push_back(media_rr);
        rm[1].push_back(mediana(dummy));
    }
    return rm;
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


bool is_muon (const caf::Proxy<caf::SRSlice>& islc, int ipfp) 
{

    if(std::isnan(islc.reco.pfp[ipfp].trk.start.x) || std::isnan(islc.reco.pfp[ipfp].trk.len)) return false;

    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 RecoStart;
    RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);

    std::vector<double> output;
    std::vector<double> dedx;
    std::vector<double> rr;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
    {   
        dedx.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
        rr.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
    } 
    output = chi2_ALG(dedx,rr,0.0,25.0);

    if( 
        islc.reco.pfp[ipfp].trackScore>0.5 &&
        output[0]<30 && output[1]>60 &&
        //( (RecoVtx-RecoStart).Mag()<10 ) && 
        islc.reco.pfp[ipfp].trk.len>20 && 
        (islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)>0 && 
        //islc.reco.pfp[ipfp].parent_is_primary
        isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) && 
        isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,5.0)
      )
    {return true;}

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

namespace param {
        bool print = false;         // print on terminal
        bool max2Seg = true;        // consider only nuMuCC with maximum 2 segments (ONLY MC)
        bool checkVtx = true;       // check if the vertex is between the two tracks

        namespace slc {
            float Tss = 0.7;        // overlap of the tracks on the z-axis
            float theta = 35.;      // angle between the directions of the tracks
            float Dss = 0.8;        // 3D distance between the tracks (ratio)
        }
        //stitching extra slice -> disabilitato
        namespace evt {
            float Tss = 0.;        // overlap of the tracks on the z-axis
            float theta = 0.;      // angle between the directions of the tracks
            float Dss = 0.;        // 3D distance between the tracks (ratio)
        }
    } // namespace param

struct PFP {
        // first track
        caf::SRVector3D vertex;      ///< Vertex of slice
        caf::SRVector3D start;       ///< Start point of track
        caf::SRVector3D end;         ///< End point of track
        float MClen;                 ///< True length of the track (MC only)
        float MCp_muon;              ///< True momentum of the track (MC only)
        float len;                   ///< Reconstructed track length
        float p_muon;                ///< Momentum estimate from trk range (muon hypothesis)
        int slcID;                   ///< Slice ID of the track
        int slice_index;             ///< Index of the slice --> NS
        std::vector<double> pid_proba;///< predicted probabilities for each class --> NS
        double chi2mu;               ///chi2 as muon --> NS
        double chi2pro;              ///chi2 as proton --> NS
        int id;                      ///< PFP ID of the track
        int pfp_index;               
        int G4ID;                    ///< G4ID of the track (MC only)
        int pdg;                     ///< PDG of the track (MC only)
        //bool muon_1muNp = false;     ///< Muon passing the 1muNp selection (data only)
        //bool nuMuCC = false;         ///< Muon from a numuCC interaction correctly stitched (MC only)

        // second track
        caf::SRVector3D start2;      ///< Start point of track
        caf::SRVector3D end2;        ///< End point of track
        float len2;                  ///< Reconstructed track length
        float p_muon2;               ///< Momentum estimate from trk range (muon hypothesis)
        int slcID2;                  ///< Slice ID of the track
        int slice_index2;            ///< Index of the slice --> NS
        std::vector<double> pid_proba2;///< predicted probabilities for each class --> NS
        double chi2mu2;               ///chi2 as muon --> NS
        double chi2pro2;              ///chi2 as proton --> NS
        int id2;                     ///< PFP ID of the track
        int pfp_index2;
        int G4ID2;                   ///< G4ID of the track (MC only)
        int pdg2;                    ///< PDG of the track (MC only)

        // stitching
        size_t nStitch = 0;          ///< Number of stitching in the same event
        float hole;                  ///< Distance between the two track edges to be stitched (0 if overlapping)
        float Len;                   ///< Length of the track after the stitching
        float P_muon;                ///< Momentum estimate from the track range (muon hypothesis) after the stitching

        // event
        unsigned int   run;       ///< run number
        unsigned int   subrun;    ///< subrun number
        unsigned int   evt;       ///< ART event number, indexes trigger windows.
    };

    // sort start and end points of the track according to the distance from the vertex
    void Icarus202503SortPFP(PFP &a, const caf::SRVector3D &vtx) {
        float dist[2];
        dist[0] = std::hypot(a.start.x-vtx.x, a.start.y-vtx.y, a.start.z-vtx.z);
        dist[1] = std::hypot(a.end.x-vtx.x, a.end.y-vtx.y, a.end.z-vtx.z);
        if(dist[1] < dist[0]) std::swap(a.start, a.end);
    }

    // sort start and end points of the track according to 3D distance and z-axis
    void Icarus202503SortPFP(PFP &a, PFP &b) {
        std::vector<float> dist,ord;
        dist.push_back(std::hypot(a.start.x-b.start.x, a.start.y-b.start.y, a.start.z-b.start.z));
        dist.push_back(std::hypot(a.start.x-b.end.x, a.start.y-b.end.y, a.start.z-b.end.z));
        dist.push_back(std::hypot(a.end.x-b.start.x, a.end.y-b.start.y, a.end.z-b.start.z));
        dist.push_back(std::hypot(a.end.x-b.end.x, a.end.y-b.end.y, a.end.z-b.end.z));
        ord = dist;
        std::sort(ord.begin(), ord.end());
        if(ord[0]==dist[0]) std::swap(a.start, a.end);
        else if(ord[0]==dist[1]) {
            std::swap(a.start, a.end);
            std::swap(b.start, b.end);
        }
        else if(ord[0]==dist[3]) std::swap(b.start, b.end);
    }

    // stitching function in the event
    std::vector<PFP> stitch(const caf::SRSpillProxy* sr, bool mc) 
    {
        std::vector<PFP> stitch;
        std::vector<PFP> pfp;
        size_t nStitch = 0;
        int muonIdx = -1;

        // loop on the slices to import the pfp
        for(int i=0; i<sr->nslc; i++) {

            //check the truth matching
            if(!is_good_truth_matching(sr->slc.at(i)))continue;

            // check if the slice is in the fiducial volume
            if(!isInFV(sr->slc.at(i).vertex.x, sr->slc.at(i).vertex.y, sr->slc.at(i).vertex.z))continue;

            for(size_t j=0; j<sr->slc.at(i).reco.npfp; j++) 
            {
                if(!is_muon(sr->slc.at(i),int(j)))continue;

                std::vector<double> output;
                std::vector<double> dedx;
                std::vector<double> rr;
                for ( std::size_t ihit(0); ihit < sr->slc.at(i).reco.pfp[j].trk.calo[2].points.size(); ++ihit )
                {   
                    dedx.push_back(sr->slc.at(i).reco.pfp[j].trk.calo[2].points[ihit].dedx);
                    rr.push_back(sr->slc.at(i).reco.pfp[j].trk.calo[2].points[ihit].rr);
                } 
                output = chi2_ALG(dedx,rr,0.0,25.0);

                PFP P;
                P.vertex.x = sr->slc.at(i).vertex.x;
                P.vertex.y = sr->slc.at(i).vertex.y;
                P.vertex.z = sr->slc.at(i).vertex.z;
                P.start.x = sr->slc.at(i).reco.pfp.at(j).trk.start.x;
                P.start.y = sr->slc.at(i).reco.pfp.at(j).trk.start.y;
                P.start.z = sr->slc.at(i).reco.pfp.at(j).trk.start.z;
                P.end.x = sr->slc.at(i).reco.pfp.at(j).trk.end.x;
                P.end.y = sr->slc.at(i).reco.pfp.at(j).trk.end.y;
                P.end.z = sr->slc.at(i).reco.pfp.at(j).trk.end.z;
                P.len = sr->slc.at(i).reco.pfp.at(j).trk.len;
                P.p_muon = sr->slc.at(i).reco.pfp.at(j).trk.rangeP.p_muon;
                P.slcID = sr->slc.at(i).reco.pfp.at(j).slcID;
                P.slice_index = i; //-->NS
                P.pid_proba = PIDproba(sr->slc.at(i),j); //-->NS
                P.chi2mu = output[0]; //-->NS
                P.chi2pro = output[1]; //-->NS
                P.pfp_index = j;
                P.id = sr->slc.at(i).reco.pfp.at(j).id;

                // MC
                if(mc) 
                {
                    if(!std::isnan(sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.G4ID)) P.G4ID = sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.G4ID;
                    else P.G4ID = -1;
                    if(!std::isnan(sr->slc.at(i).reco.pfp.at(j).trk.truth.p.pdg)) P.pdg = sr->slc.at(i).reco.pfp.at(j).trk.truth.p.pdg;
                    else P.pdg = -1;
                    if(!std::isnan(sr->slc.at(i).reco.pfp.at(j).trk.truth.p.length)) P.MClen = sr->slc.at(i).reco.pfp.at(j).trk.truth.p.length;
                    else P.MClen = -1;
                }
                pfp.push_back(P);
            }
        }

        // sort the pfp according to the distance from the vertex
        std::sort(pfp.begin(), pfp.end(), [](const PFP &a, const PFP &b) {
            caf::SRVector3D baryA, baryB;
            baryA.x = 0.5*(a.start.x+a.end.x);
            baryA.y = 0.5*(a.start.y+a.end.y);
            baryA.z = 0.5*(a.start.z+a.end.z);
            baryB.x = 0.5*(b.start.x+b.end.x);
            baryB.y = 0.5*(b.start.y+b.end.y);
            baryB.z = 0.5*(b.start.z+b.end.z);
            float distA = std::hypot(baryA.x-a.vertex.x, baryA.y-a.vertex.y, baryA.z-a.vertex.z);
            float distB = std::hypot(baryB.x-a.vertex.x, baryB.y-a.vertex.y, baryB.z-a.vertex.z);
            return distA < distB;
        });
        // stitching
        for(size_t i=0; i<pfp.size(); i++) {
            for(size_t j=i+1; j<pfp.size(); j++) {

                std::vector<PFP> trk;
                caf::SRVector3D vtx;
                vtx.x = pfp.at(i).vertex.x;
                vtx.y = pfp.at(i).vertex.y;
                vtx.z = pfp.at(i).vertex.z;
                // sort vertex and tracks according to the baycenters
                std::vector<float> ordz;
                ordz.push_back(vtx.z);                                       // vertex
                ordz.push_back( 0.5*(pfp.at(i).start.z+pfp.at(i).end.z) );   // (start 1 + end 1)/2
                ordz.push_back( 0.5*(pfp.at(j).start.z+pfp.at(j).end.z) );   // (start 2 + end 2)/2
                std::sort(ordz.begin(), ordz.end());
                // if the first vertex is in the middle with tracks
                // of different slices, try again with the second vertex
                if(ordz[1] == vtx.z && pfp.at(i).slcID!=pfp.at(j).slcID) {
                    ordz.clear();
                    vtx.x = pfp.at(j).vertex.x;
                    vtx.y = pfp.at(j).vertex.y;
                    vtx.z = pfp.at(j).vertex.z;
                    ordz.push_back(vtx.z);                                       // vertex
                    ordz.push_back( 0.5*(pfp.at(i).start.z+pfp.at(i).end.z) );   // (start 1 + end 1)/2
                    ordz.push_back( 0.5*(pfp.at(j).start.z+pfp.at(j).end.z) );   // (start 2 + end 2)/2
                    std::sort(ordz.begin(), ordz.end());
                }
                // check vertex in the middle
                if(param::checkVtx && ordz[1] == vtx.z) {continue;}
                // check which segment is closest to the vertex
                if(ordz[1]==0.5*(pfp.at(i).start.z+pfp.at(i).end.z)) {
                    trk.push_back(pfp.at(i));
                    trk.push_back(pfp.at(j));
                    // sort start and end points of each segment according to distance from the vertex
                    for(size_t k=0; k<2; k++) Icarus202503SortPFP(trk.at(k), vtx);
                }
                else if(ordz[1]==0.5*(pfp.at(j).start.z+pfp.at(j).end.z)) {
                    trk.push_back(pfp.at(j));
                    trk.push_back(pfp.at(i));
                    // sort start and end points of each segment according to distance from the vertex
                    for(size_t k=0; k<2; k++) Icarus202503SortPFP(trk.at(k), vtx);
                }
                else if(ordz[0]==0.5*(pfp.at(i).start.z+pfp.at(i).end.z)) {
                    trk.push_back(pfp.at(i));
                    trk.push_back(pfp.at(j));
                    // sort start and end end point of each track according to the 3D distance and z-axis
                    Icarus202503SortPFP(trk.at(0), trk.at(1));
                }
                else {
                    trk.push_back(pfp.at(j));
                    trk.push_back(pfp.at(i));
                    // sort start and end end point of each track according to the 3D distance and z-axis
                    Icarus202503SortPFP(trk.at(0), trk.at(1));
                }
                // check overlap on z-axis for tracks in the same slice
                std::vector<float> zlen;
                zlen.push_back(std::abs(trk.at(0).end.z-trk.at(0).start.z));
                zlen.push_back(std::abs(trk.at(1).end.z-trk.at(1).start.z));
                std::sort(zlen.begin(),zlen.end());
                bool Zoverlap = (ordz[0]==vtx.z && trk.at(1).start.z < trk.at(0).end.z) || (ordz[2]==vtx.z && trk.at(1).start.z > trk.at(0).end.z);
                if(Zoverlap) 
                {
                    if(trk.at(0).slcID==trk.at(1).slcID) 
                    {
                        // selection
                        if(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= param::slc::Tss) {continue;}
                    }
                    else 
                    {
                        // selection
                        if(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= param::evt::Tss) {continue;}
                    }
                }
                // check angle between tracks
                float a[3],b[3];
                a[0] = trk.at(0).end.x-trk.at(0).start.x;
                a[1] = trk.at(0).end.y-trk.at(0).start.y;
                a[2] = trk.at(0).end.z-trk.at(0).start.z;
                b[0] = trk.at(1).end.x-trk.at(1).start.x;
                b[1] = trk.at(1).end.y-trk.at(1).start.y;
                b[2] = trk.at(1).end.z-trk.at(1).start.z;
                float Ra = std::hypot(a[0],a[1],a[2]);
                float Rb = std::hypot(b[0],b[1],b[2]);
                float c = (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/(Ra*Rb); // scalar product
                float phi = acos(c); // [rad]
                if(trk.at(0).slcID==trk.at(1).slcID) 
                {
                    // selection
                    if(phi*180/3.141592 >= param::slc::theta){continue;}
                }
                else {
                    // selection
                    if(phi*180/3.141592 >= param::evt::theta){continue;}
                }
                // check distance between tracks
                std::vector<float> length;
                length.push_back(Ra);
                length.push_back(Rb);
                std::sort(length.begin(), length.end());
                float dist = std::hypot(trk.at(1).start.x-trk.at(0).end.x, trk.at(1).start.y-trk.at(0).end.y, trk.at(1).start.z-trk.at(0).end.z);
                if(trk.at(0).slcID==trk.at(1).slcID) 
                {
                    // selection
                    if(dist/length[0] >= param::slc::Dss){continue;}
                }
                else 
                {
                    // selection
                    if(dist/length[0] >= param::evt::Dss){continue;}
                }

                nStitch++;
                PFP sth;
                // track 1
                sth.vertex.x = vtx.x;                             // vertex x
                sth.vertex.y = vtx.y;                             // vertex y
                sth.vertex.z = vtx.z;                             // vertex z
                sth.start.x = trk.at(0).start.x;                  // start x
                sth.start.y = trk.at(0).start.y;                  // start y
                sth.start.z = trk.at(0).start.z;                  // start z
                sth.end.x = trk.at(0).end.x;                      // end x
                sth.end.y = trk.at(0).end.y;                      // end y
                sth.end.z = trk.at(0).end.z;                      // end z
                sth.len = trk.at(0).len;                          // length
                sth.p_muon = trk.at(0).p_muon;                    // momentum
                sth.slcID = trk.at(0).slcID;                      // slcID
                sth.id = trk.at(0).id;                            // pfpID
                sth.pfp_index = trk.at(0).pfp_index;
                // track 2
                sth.start2.x = trk.at(1).start.x;                 // start x
                sth.start2.y = trk.at(1).start.y;                 // start y
                sth.start2.z = trk.at(1).start.z;                 // start z
                sth.end2.x = trk.at(1).end.x;                     // end x
                sth.end2.y = trk.at(1).end.y;                     // end y
                sth.end2.z = trk.at(1).end.z;                     // end z
                sth.len2 = trk.at(1).len;                         // length
                sth.p_muon2 = trk.at(1).p_muon;                   // momentum
                sth.slcID2 = trk.at(1).slcID;                     // slcID
                sth.id2 = trk.at(1).id;                           // pfpID
                sth.pfp_index2 = trk.at(1).pfp_index;
                // stitching
                sth.nStitch = nStitch;                            // number of stitching in the same event
                if(Zoverlap) {
                    // distance between the two track edges to be stitched
                    sth.hole = std::hypot(trk.at(1).start.x-trk.at(0).end.x,
                                            trk.at(1).start.y-trk.at(0).end.y,
                                                trk.at(1).start.z-trk.at(0).end.z);
                    // total length after the stitching
                    sth.Len = sth.len + sth.len2 + sth.hole;
                }
                else {
                    sth.hole = 0.;
                    sth.Len = sth.len + sth.len2;                 // total length after the stitching
                }
                sth.P_muon = GetTrackMomentum(sth.Len, 13);       // total momentum after the stitching
                // MC or data
                if(mc)
                {
                    sth.MClen = trk.at(0).MClen;                        // MC length
                    if(sth.MClen > 0) sth.MCp_muon = GetTrackMomentum(sth.MClen, 13);    // MC momentum
                    else sth.MCp_muon = -1;
                    sth.G4ID = trk.at(0).G4ID;                          // G4ID 1
                    sth.pdg = trk.at(0).pdg;                            // pdg 1
                    sth.G4ID2 = trk.at(1).G4ID;                         // G4ID 2
                    sth.pdg2 = trk.at(1).pdg;                           // pdg 2
                }
                // event
                sth.slice_index = trk.at(0).slice_index;
                sth.slice_index2 = trk.at(1).slice_index;
                sth.run = sr->hdr.run;                             // run number
                sth.subrun = sr->hdr.subrun;                       // subrun number
                sth.evt = sr->hdr.evt;                             // event number
                stitch.push_back(sth);
            }
        }
        return stitch;
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

//////////////////////////////////////////////////////////////////////////////////
///////////                   DEBUG 1mu0p0pi                          ////////////
//////////////////////////////////////////////////////////////////////////////////

std::vector<int> get_protons_pfp ( const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc){

  std::vector<int> v_ipfp_pro;

  int ipfp_mu = -1;
  ipfp_mu = find_muon_newPID(islc,10,-1,-1,-1);

  for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
  {
    if(int(ipfp)==ipfp_mu)continue;
    if(id_pfp_newPID(islc,ipfp,10)==1)
    { 
      v_ipfp_pro.push_back(ipfp);
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


//ofstream debug1mu0p0pi("debug1mu0p0pi_dedxmag1_newtrain_onlycoll.txt"); 
//ofstream debug1muNp("debug1muNp_dedxmag1_newtrain_onlycoll.txt");

ofstream debug1mu0p0pi("debug1mu0p0pi_features.txt");
ofstream debug1muNp("debug1muNp_feature.txt");

const SpillMultiVar fdebug1mu0p0pi([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  std::vector<double> vector_active;

  int slice_counter = -1;
  for(const auto &islc : sr->slc)
  {
    slice_counter ++;
    std::string slice_class_reco = selection_1muNp_newPID(sr,islc,10,100);

    std::string slice_class_true = classification_type_generic(sr,islc);

    if(is1muNp(slice_class_reco)==true && is1mu_true(slice_class_true)==true)
    {
      std::vector<int> v_ipfp_pro = get_protons_pfp(sr,islc);

        for(const auto &ipfp : v_ipfp_pro)
        {
            if(islc.reco.pfp[ipfp].trk.truth.p.pdg == 2212)continue;

            int bestplane = find_best_plane(islc,ipfp);

            //std::vector<double> prediction_proba;
            //prediction_proba = PIDproba(islc,ipfp);

            std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
            double depE = compute_depE(islc,ipfp,bestplane);
            std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

            std::vector<double> track_features;

            track_features.insert(track_features.end(), lr.begin(), lr.end());
            track_features.push_back(depE);
            track_features.insert(track_features.end(), dvars.begin(), dvars.end());

            debug1mu0p0pi << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
            for(const auto &feature : track_features){debug1mu0p0pi << feature << " ";}
            debug1mu0p0pi << endl;

            /*
            debug1mu0p0pi << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
            debug1mu0p0pi << bestplane << " " << get_n_valid_hit(islc,ipfp) << " " << islc.reco.pfp[ipfp].trk.len << " " << compute_depE(islc,ipfp,bestplane) << " " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " ";

            for(const auto &prob : prediction_proba)
            {
              debug1mu0p0pi << prob << "_";
            }

            debug1mu0p0pi << " ";
            
            for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
            {
              debug1mu0p0pi << islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr << "_";
            } 

            debug1mu0p0pi << " ";

            for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
            {
              debug1mu0p0pi << islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx << "_";
            } 
            
            debug1mu0p0pi << endl;
            */
        }
    }

    
    if(is1muNp(slice_class_reco)==true && is1muNp_true(slice_class_true)==true)
    {
      std::vector<int> v_ipfp_pro = get_protons_pfp(sr,islc);

        for(const auto &ipfp : v_ipfp_pro)
        {
            int bestplane = find_best_plane(islc,ipfp);

            //std::vector<double> prediction_proba;
            //prediction_proba = PIDproba(islc,ipfp);

            std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
            double depE = compute_depE(islc,ipfp,bestplane);
            std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

            std::vector<double> track_features;

            track_features.insert(track_features.end(), lr.begin(), lr.end());
            track_features.push_back(depE);
            track_features.insert(track_features.end(), dvars.begin(), dvars.end());

            debug1muNp << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
            for(const auto &feature : track_features){debug1muNp << feature << " ";}
            debug1muNp << endl;

            /*
            debug1muNp << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " ";
            debug1muNp << bestplane << " " << get_n_valid_hit(islc,ipfp) << " " << islc.reco.pfp[ipfp].trk.len << " " << compute_depE(islc,ipfp,bestplane) << " " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " ";
            
            for(const auto &prob : prediction_proba)
            {
              debug1muNp << prob << "_";
            }

            debug1muNp << " ";

            for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
            {
              debug1muNp << islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr << "_";
            } 

            debug1muNp << " ";

            for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
            {
              debug1muNp << islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx << "_";
            } 
            
            debug1muNp << endl;
            */
        }
    }
    
  }

  return vector_active;

});


///////////////////////////////////////////////////////////////////////////////////
//////////                PRED PROBA DISTRIBUTIONS                /////////////////
///////////////////////////////////////////////////////////////////////////////////

/*
struct particle
{
  int run;
  int evt;
  int slice_index;
  int bestplane;
  int nhit;
  int pdg;
  double length;
  double depE;
  //double theta_drift;
  std::array<double,3> dir;
  std::vector<double> v_pred_proba;
  std::vector<double> rr;
  std::vector<double> dedx;
};

std::vector<particle> muons;
std::vector<particle> protons;
std::vector<particle> pions;

int nmuons=0;
int npions=0;
int nprotons=0;

const SpillMultiVar fdump_pred_proba([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  std::vector<double> vector_active;

  int slice_counter = -1;
  for(const auto &islc : sr->slc)
  {
    slice_counter ++;

    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    
    if (islc.tmatch.eff <= 0.5 || (vertex_true-vertex_reco).Mag() >= 100 || islc.is_clear_cosmic)continue;

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
    {
      std::vector<double> prediction_proba;
      prediction_proba = PIDproba(islc,ipfp);

      double length = islc.reco.pfp[ipfp].trk.len;

      int bestplane = find_best_plane(islc,ipfp);

      int pdg = std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg);

      int nhit = get_n_valid_hit(islc,ipfp);

      int predicted_class = PIDclass(islc,ipfp);

      double depE_temp = -1;

      //double theta_drift = 180./M_PI*acos(abs(islc.reco.pfp[ipfp].trk.dir.x));

      std::array<double,3> dir = {islc.reco.pfp[ipfp].trk.dir.x,islc.reco.pfp[ipfp].trk.dir.y,islc.reco.pfp[ipfp].trk.dir.z};

      std::vector<double> rr_temp;
      std::vector<double> dedx_temp;

      if(bestplane !=-1)
      {
        depE_temp = compute_depE(islc,ipfp,bestplane);
      
        for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
        {
          {
            rr_temp.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr);
            dedx_temp.push_back(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx);
          }
        }
      }

      particle par;
      par.run = sr->hdr.run;
      par.evt = sr->hdr.evt;
      par.slice_index = slice_counter;
      par.bestplane = bestplane;
      par.nhit = nhit;
      par.pdg = pdg;
      par.length = length;
      par.depE = depE_temp;
      par.dir = dir;
      par.v_pred_proba = prediction_proba;
      par.rr = rr_temp;
      par.dedx = dedx_temp;

      double prob_as_muon = prediction_proba[0]+prediction_proba[1];
      double prob_as_proton = prediction_proba[2]+prediction_proba[3];
      double prob_as_pion = prediction_proba[4]+prediction_proba[5];

      //if(predicted_class == 0 || predicted_class == 1)
      if( (prob_as_muon > prob_as_proton) && (prob_as_muon > prob_as_pion) )
      {
        muons.push_back(par);
        nmuons ++;
      }
      //if(predicted_class == 4 || predicted_class == 5)
      if( (prob_as_pion > prob_as_muon) && (prob_as_pion > prob_as_proton) )
      {
        pions.push_back(par);
        npions ++;
      }
      //if(predicted_class == 2 || predicted_class == 3)
      if( (prob_as_proton > prob_as_muon) && (prob_as_proton > prob_as_pion) )
      {
        protons.push_back(par);
        nprotons ++;
      }
      
    }
  }

  return vector_active;

});
*/

///////////////////////////////////////////////////////////////////////////////////
//////////                           SELECTION                           //////////
///////////////////////////////////////////////////////////////////////////////////

/*
std::vector<std::string> slices_reco_class;
std::vector<std::string> slices_true_class;
int tot_1muNp;
//ofstream dumpNuE("dumpNuE_temp.txt");

bool ismc;

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
      double nuE = -1;
      if (sr->hdr.ismc)
      {
        nuE = islc.truth.E;
        slice_class_true = classification_type_generic(sr,islc);
      }
      slices_true_class.push_back(slice_class_true);

      std::string slice_class = selection_1muNp_newPID(sr,islc,10,100,-1,-1,-1);
      slices_reco_class.push_back(slice_class);

    }//loop over all slices

  return vector_active;
});
*/

///////////////////////////////////////////////////////////////////////////////////
//////////                     SELECTION WITH STITCHING                  //////////
///////////////////////////////////////////////////////////////////////////////////

/*
ofstream dump_stitching_AMR("dump_stitching_AMR_light.txt");

    // stitching function
    const SpillMultiVar kStitch([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<double> vector_active;
        sth = stitch(sr, sr->hdr.ismc);

        //aggiunta da nicola ////////////////
        for(size_t i=0; i<sth.size(); i++) 
        {
            dump_stitching_AMR << sth.at(i).run << " " << sth.at(i).evt << " " << sth.at(i).slice_index << " " << sth.at(i).G4ID << " " << sth.at(i).pdg << " " << sth.at(i).G4ID2 << " " << sth.at(i).pdg2 << " " << sth.at(i).len << " " << sth.at(i).Len << " " << sth.at(i).MClen << endl;
        }
        /////////////////////////////////////

        return vector_active;
    });

*/
/*
std::vector<std::string> slices_reco_class;
std::vector<std::string> slices_true_class;
int tot_1muNp;
ofstream dumpNuE("dumpNuE_temp.txt");
ofstream dump_reco_class("dump_reco_class_temp.txt");

bool ismc;

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
        dumpNuE << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " " << nu.E << " " << nu.index << endl;
        if(nu.index == 0)tot_1muNp++;
      }
    }
  }

    std::vector<double> vector_active;

    //stitching su tutte le slices dello spill
    std::vector<PFP> sth;
    //sth = stitch(sr, sr->hdr.ismc);

    //for(size_t i=0; i<sth.size(); i++) 
    //{
      //cout << sth.at(i).run << " " << sth.at(i).evt << " " << sth.at(i).slice_index << " " << sth.at(i).G4ID << " " << sth.at(i).pdg << " " << sth.at(i).G4ID2 << " " << sth.at(i).pdg2 << " " << sth.at(i).len << " " << sth.at(i).Len << " " << sth.at(i).MClen << endl;
    //}
    //cout << "end of stitching for this spill" << endl;

    int slice_counter = -1;
    for(const auto &islc : sr->slc)
    { 
      slice_counter ++;

      std::string slice_class_true = "N/D";
      double nuE = -1;
      if (sr->hdr.ismc)
      {
        nuE = islc.truth.E;
        slice_class_true = classification_type_generic(sr,islc);
      }
      slices_true_class.push_back(slice_class_true);

      int mu1_pfp = -1;
      int mu2_pfp = -1;
      double stitched_len = -1;
      for(const auto &s : sth)
      {
        if((int)s.run == (int)sr->hdr.run && (int)s.evt == (int)sr->hdr.evt && s.slice_index == slice_counter)
        {
          mu1_pfp = s.pfp_index;
          mu2_pfp = s.pfp_index2;
          stitched_len = s.Len;
        }
      }

      std::string slice_class = selection_1muNp_newPID(sr,islc,10,100,mu1_pfp,mu2_pfp,stitched_len);
      slices_reco_class.push_back(slice_class);

      //dump_reco_class << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " " << nuE << " " << slice_class << " " << slice_class_true << endl;

      if(is1muNp(slice_class))
      {
        dump_reco_class << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " " << slice_class << " " << std::string(sr->hdr.sourceName) << endl;
      }

    }//loop over all slices

  return vector_active;
});
*/


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////                                                 DAUGHTERS STITCHING                                           //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
std::vector<std::string> slices_reco_class;

std::vector<double> len_reco_stitched;
std::vector<double> len_true;
std::vector<double> len_reco_before_stitching;
std::vector<std::array<double,3>> end_reco_stitched;
std::vector<std::array<double,3>> end_true;
std::vector<std::array<double,3>> end_reco_before_stitching;
std::vector<bool> correctly_stitched; 

ofstream dump_stitching("dump_stitching_new.txt");
ofstream dump_stitching_NS("dump_stitching_NS_angle_new.txt");

const SpillMultiVar selection_newPID([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    int slice_counter = -1;
    for(const auto &islc : sr->slc)
    { 
      slice_counter ++;

      if(!isInFV(islc.vertex.x,islc.vertex.y,islc.vertex.z))continue;
      if(!is_good_truth_matching(islc))continue;

      //std::string slice_class = selection_1muNp_newPID(sr,islc,10,100);
      //slices_reco_class.push_back(slice_class);

      //if(!is1muNp(slice_class))continue;

      std::string slice_class_true = classification_type_generic(sr,islc);

      if(!is1muNp_true(slice_class_true))continue;

      int G4ID_true_muon = -1;
      std::array<double,3> end_true_muon;
      double len_true_muon=-1;
      double muon_momentum = -1;
      for(const auto & iprim : islc.truth.prim)
      {
        if(std::abs(iprim.pdg)==13)
        {
          G4ID_true_muon = iprim.G4ID;
          len_true_muon = iprim.length;
          end_true_muon = {iprim.end.x, iprim.end.y, iprim.end.z};
          muon_momentum = std::sqrt(std::pow(iprim.genp.x,2)+std::pow(iprim.genp.y,2)+std::pow(iprim.genp.z,2));
        }
      }

      for( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
      {
          if(!is_muon(islc,ipfp))continue;

          int pid_class = PIDclass(islc,ipfp);

          //asking for a muon mip at least 50 cm long
          //if(pid_class==1)
          //{
            //search for the best daughter having more hits in the last 5 cm
            std::array<double,3> track_direction_end = {islc.reco.pfp[ipfp].trk.dir_end.x, islc.reco.pfp[ipfp].trk.dir_end.y, islc.reco.pfp[ipfp].trk.dir_end.z};
            int best_daughter_pfp = -1;
            int best_daughter_hits = 0;
            double best_daughter_depE = 0;
            double best_angle = -1;
            for(int daughter=0; daughter<int(islc.reco.pfp[ipfp].ndaughters); daughter++)
            {  
              int d = islc.reco.pfp[ipfp].daughters[daughter];
              int daughter_hits = 0;
              double temp_depE = 0;

              for (std::size_t jpfp(0); jpfp < islc.reco.npfp; ++jpfp)
              {
                if(islc.reco.pfp[jpfp].id!=d)continue;

                int d_bestplane = islc.reco.pfp[jpfp].trk.bestplane;

                if(d_bestplane==-1)continue; 

                std::array<double,3> daughter_direction = {islc.reco.pfp[jpfp].trk.dir.x, islc.reco.pfp[jpfp].trk.dir.y, islc.reco.pfp[jpfp].trk.dir.z};
                double angle_end = std::acos(track_direction_end[0]*daughter_direction[0] + track_direction_end[1]*daughter_direction[1] + track_direction_end[2]*daughter_direction[2]);

                for ( std::size_t ihit(0); ihit < islc.reco.pfp[jpfp].trk.calo[d_bestplane].points.size(); ++ihit )
                {
                  if(islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].rr <= 5)
                  {  
                    daughter_hits++; 
                    temp_depE = temp_depE + islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].dedx*islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].pitch;
                  }
                }
                if(daughter_hits>best_daughter_hits)
                {
                  best_daughter_hits = daughter_hits;
                  best_daughter_pfp = jpfp;
                  best_daughter_depE = temp_depE;
                  best_angle = angle_end;
                }
              }
            }//loop over all pfp's daughters

            //searching for muon rising secondaries
            if(pid_class == 1 && best_daughter_depE < 50 && best_daughter_depE > 20 && std::abs(best_angle)*180/M_PI<45 && islc.reco.pfp[best_daughter_pfp].trk.len > 20 && isInContained(islc.reco.pfp[best_daughter_pfp].trk.end.x, islc.reco.pfp[best_daughter_pfp].trk.end.y, islc.reco.pfp[best_daughter_pfp].trk.end.z, 5))
            //if(best_daughter_pfp!=-1 && (pid_class == 1 || pid_class == 5) && PIDclass(islc,best_daughter_pfp)==0 && islc.reco.pfp[best_daughter_pfp].trk.len > 20)
            {
              //if(islc.reco.pfp[ipfp].trk.truth.p.pdg==2212){cout << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << endl; }

              if(islc.reco.pfp[ipfp].trk.truth.p.G4ID == G4ID_true_muon && islc.reco.pfp[best_daughter_pfp].trk.truth.p.G4ID == G4ID_true_muon){correctly_stitched.push_back(true);}
              else{correctly_stitched.push_back(false);}

              len_true.push_back(len_true_muon);
              end_true.push_back(end_true_muon);

              len_reco_stitched.push_back(islc.reco.pfp[best_daughter_pfp].trk.len + islc.reco.pfp[ipfp].trk.len);
              end_reco_stitched.push_back({islc.reco.pfp[best_daughter_pfp].trk.end.x, islc.reco.pfp[best_daughter_pfp].trk.end.y, islc.reco.pfp[best_daughter_pfp].trk.end.z});

              end_reco_before_stitching.push_back({islc.reco.pfp[ipfp].trk.end.x, islc.reco.pfp[ipfp].trk.end.y, islc.reco.pfp[ipfp].trk.end.z});
              len_reco_before_stitching.push_back(islc.reco.pfp[ipfp].trk.len);

              TVector3 reco_muon_momentum;
              reco_muon_momentum.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.z);


              dump_stitching << end_true_muon[0] << " " <<  end_true_muon[1] << " " << end_true_muon[2] << " ";
              dump_stitching << islc.reco.pfp[ipfp].trk.end.x << " " << islc.reco.pfp[ipfp].trk.end.y << " " << islc.reco.pfp[ipfp].trk.end.z << " ";
              dump_stitching << islc.reco.pfp[best_daughter_pfp].trk.end.x << " " << islc.reco.pfp[best_daughter_pfp].trk.end.y << " " << islc.reco.pfp[best_daughter_pfp].trk.end.z << " ";
              dump_stitching << len_true_muon << " " << islc.reco.pfp[ipfp].trk.len << " " << islc.reco.pfp[best_daughter_pfp].trk.len + islc.reco.pfp[ipfp].trk.len << " ";
              dump_stitching << muon_momentum << " " << reco_muon_momentum.Mag() << " " << GetTrackMomentum(islc.reco.pfp[best_daughter_pfp].trk.len + islc.reco.pfp[ipfp].trk.len, 13) << endl;

              dump_stitching_NS << sr->hdr.run << " " << sr->hdr.evt << " " << slice_counter << " " << islc.reco.pfp[ipfp].trk.truth.p.G4ID << " " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[best_daughter_pfp].trk.truth.p.G4ID << " " << islc.reco.pfp[best_daughter_pfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << islc.reco.pfp[best_daughter_pfp].trk.len + islc.reco.pfp[ipfp].trk.len << " " << islc.reco.pfp[ipfp].trk.truth.p.length << endl;
            }

          //}
      }

    }//loop over all slices

  return vector_active;
});
*/

///////////////////////////////////////////////////////////////////
/////////                 LR ON DATA                    ///////////
///////////////////////////////////////////////////////////////////

/*
ofstream lr_data("lr_data.txt");
ofstream dedx_data("dedx_data.txt");

const SpillMultiVar DATAlikelihood([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for(const auto &islc : sr->slc)
    {
      if(islc.is_clear_cosmic)continue;

      for( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
      {
        if(islc.reco.pfp[ipfp].trk.len > 150)
        {
          int bestplane = find_best_plane(islc,ipfp);
          if(bestplane!=2)continue;
          std::vector<double> likelihood_ratios = likelihood(islc,ipfp,prob_d_coll, prob_d_ind1, prob_d_ind2);

          if(isInContained(islc.reco.pfp[ipfp].trk.end.x, islc.reco.pfp[ipfp].trk.end.y, islc.reco.pfp[ipfp].trk.end.z,5))
          {
            lr_data << "stopping ";
            for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
            {
              dedx_data << "stopping " << islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr << " " << islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx << endl;
            }
          }
          else
          {
            lr_data << "exiting ";
            for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
            {
              dedx_data << "exiting " << islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr << " " << islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx << endl;
            }
          }
          
          for(const auto &lr : likelihood_ratios){lr_data << lr << " ";}
          lr_data << endl; 
        }
      }
      
    }

    return vector_active;
});
*/


//ofstream chi2dump("chi2dump.txt");
//ofstream depEdump("depEdump_respun.txt");
//ofstream depEdump("depEdump_highLT.txt");
//ofstream depEdump("depEdump_null.txt");

//ofstream proton_depEdump("proton_depEdump_respun.txt");

/*
ofstream proton_depEdump("proton_depEdump_highLT.txt");

const SpillMultiVar fchi2dump([](const caf::SRSpillProxy* sr)-> std::vector<double>
{

  std::vector<double> vector_active;

  for(const auto &islc : sr->slc)
  { 
    //if(!is_good_truth_matching(islc))continue;

    //bool is_true_1muNp = false;
    //std::string slice_class_true = classification_type_generic(sr,islc);

    //is_true_1muNp = is1muNp_true(slice_class_true);

    TVector3 RecoVtx;
    RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);

    for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp )
    {
      TVector3 RecoStart;
      RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);

      if(!islc.reco.pfp[ipfp].parent_is_primary)continue;

      if(!isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0))continue; 
        
      if((islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x)<0)continue;

      if((RecoVtx-RecoStart).Mag()>=10)continue;



      //muon search 

      //if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness >= 0.5 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity >= 0.5)
      //{
        if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==2212 && islc.reco.pfp[ipfp].trk.len > 10)
        {
          
          double depE_25 = 0;
          double depE_10 = 0;
          double depE_5 = 0;
          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
          {   
            if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr < 25)
            {
              depE_25 = depE_25 + islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[2].points[ihit].pitch;
            }
            if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr < 10)
            {
              depE_10 = depE_10 + islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[2].points[ihit].pitch;
            }
            if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr < 5)
            {
              depE_5 = depE_5 + islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[2].points[ihit].pitch;
            }
             
          } 

          if(depE_5 == 0 || depE_10 == 0 || depE_25 == 0)continue;

          proton_depEdump << depE_5 << " " << depE_10 << " " << depE_25 << endl;

          //chi2dump << "muon " << output[0] << " " << output[1] << endl;
        }
      //}

      
      //proton search

      if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness >= 0.3 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity >= 0.3)
      {
        TVector3 start_mom;
        start_mom.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
        double ke = sqrt(pow(938.3,2)+pow(start_mom.Mag()*1000,2))-938.3;

        if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==2212 && ke >= 50)
        {
          std::vector<double> output;
          std::vector<double> dedx;
          std::vector<double> rr;
          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
          {   
            dedx.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
            rr.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
          } 
          output = chi2_ALG(dedx,rr,0.0,25.0);

          chi2dump << "proton " << output[0] << " " << output[1] << endl;
        }
      }

      //pion search

      if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness >= 0.3 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity >= 0.3)
      {
        TVector3 start_mom;
        start_mom.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);
        double ke = sqrt(pow(139.570,2)+pow(start_mom.Mag()*1000,2))-139.570;

        if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211 && ke >= 50)
        {
          std::vector<double> output;
          std::vector<double> dedx;
          std::vector<double> rr;
          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
          {   
            dedx.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
            rr.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
          } 
          output = chi2_ALG(dedx,rr,0.0,25.0);

          chi2dump << "pion " << output[0] << " " << output[1] << endl;
        }
      }

      

    }  

  }

  return vector_active;

});

*/