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

#include "data_struct.h"

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

Double_t correction_function_25only(Double_t *x, Double_t *par)
{
    if(x[0]<=2) return 0.0127159;
    else if (x[0]>=15) return -0.0235248;
    else return exp(par[0]+par[1]*x[0])+par[2]+par[3]*exp(-0.5*pow((x[0]-par[4])/par[5],2));
}
TF1 *correction_function = new TF1("correction_function",correction_function_25only, 0., 16., 6 );


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

int true_selection(const caf::SRSpillProxy* sr, const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
    //-1 : unclassified
    //0 : muon rising
    //1 : muon mip
    //2 : proton rising
    //3 : proton interacting
    //4 : pion rising
    //5 : pion interacting 

    int bestplane = islc.reco.pfp[ipfp].trk.bestplane;
    if(bestplane==-1)return -1;

    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()==0)return -1;

    //controllo che vertice reco e true siano abbastanza vicini per eliminare la componente di cosmici
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////                                                 Data Struct                                                  //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//std::vector<RecoSlice> slices; 
int tot_1muNp=0;

std::array<std::vector<TH1D*>,6> prob_d_coll;
std::array<std::vector<TH1D*>,6> prob_d_ind1;
std::array<std::vector<TH1D*>,6> prob_d_ind2;

std::array<std::vector<TH1D*>,6> prob_d_coll_150;
std::array<std::vector<TH1D*>,6> prob_d_ind1_150;
std::array<std::vector<TH1D*>,6> prob_d_ind2_150;

std::string txt_name = "/exp/icarus/data/users/nsommagg/data_struct_multiplane/dumpNuE_provaMaria.txt";
ofstream dumpNuE(txt_name.c_str());

TTree * tree;
RecoSlice thisslice;

const SpillMultiVar DataLoader([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  std::vector<double> vector_active;
  //if(sr->hdr.run==9435){cout << "trovato" << endl;}
  if(sr->hdr.run>0){ 

  //count how many true 1muNp neutrino interactions
  int nu_multiplicity=0;
  for ( auto const& nu : sr->mc.nu )
  {
    bool is_1muNp = classification_type_MC(sr,nu);
    if(is_1muNp)
    {
      tot_1muNp++;
      dumpNuE << sr->hdr.run << " " << sr->hdr.evt << " " << nu.position.x << " " << nu.position.y << " " << nu.position.z << " " << nu.E << " " << nu.index << " " << nu_multiplicity << endl;
      nu_multiplicity = nu_multiplicity + 1;
      //if(nu.index != 0) cout << "FOUND" << endl; 
    }
  }

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
      itrack.reco.nhits_coll = islc.reco.pfp[ipfp].trk.calo[2].points.size();
      
      itrack.reco.likelihood_ratios = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
      itrack.reco.likelihood_ratios_mediumbinning = likelihood(islc,ipfp, prob_d_coll_150, prob_d_ind1_150, prob_d_ind2_150);

      

      double depE=0;
      int bestplane = find_best_plane(islc,ipfp);
      if(bestplane!=-1){depE = compute_depE(islc,ipfp,bestplane);}
      /*
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit)
      {
        if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr>5)continue;
        depE = depE + islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[2].points[ihit].pitch;
      }

      if(depE==0)
      {
	      int use_plane = -1;
   	    if(islc.reco.pfp[ipfp].trk.bestplane==0 || islc.reco.pfp[ipfp].trk.bestplane==1 || islc.reco.pfp[ipfp].trk.bestplane==2){use_plane=islc.reco.pfp[ipfp].trk.bestplane;}
   	    else{use_plane=2;}
        for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit)
          {
            if(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr>5)continue;
            depE = depE + islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].pitch;
          }
      }
      */

      itrack.reco.deposited_energy = depE;
      itrack.truth.pdg = islc.reco.pfp[ipfp].trk.truth.p.pdg;
      itrack.truth.true_class = true_selection(sr,islc,ipfp);
      //std::vector<double> l = likelihood(islc,ipfp);
      //if(true_selection(sr,islc,ipfp)==0)
      //{
        //cout << l.size();
        //for(int i=0; i<(int)l.size(); i++){cout << l[i] << " ";}
        //cout << endl;
      //}


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

          int bestplane = islc.reco.pfp[jpfp].trk.bestplane;

          if(bestplane==-1)continue;

          for ( std::size_t ihit(0); ihit < islc.reco.pfp[jpfp].trk.calo[bestplane].points.size(); ++ihit )
          {
            if(islc.reco.pfp[jpfp].trk.calo[bestplane].points[ihit].rr <= 5)
            {
              temp_rr.push_back(islc.reco.pfp[jpfp].trk.calo[bestplane].points[ihit].rr);
              temp_dedx.push_back(islc.reco.pfp[jpfp].trk.calo[bestplane].points[ihit].dedx);
              temp_depE = temp_depE + islc.reco.pfp[jpfp].trk.calo[bestplane].points[ihit].dedx*islc.reco.pfp[jpfp].trk.calo[bestplane].points[ihit].pitch;
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
        itrack.reco.daughter_depE = -1;
        itrack.reco.daughter_angle_end = -1;
      }
      else
      {
        itrack.reco.daughter_depE = daughter_depE;
        itrack.reco.daughter_angle_end = angle_end;
      }
      
      
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
