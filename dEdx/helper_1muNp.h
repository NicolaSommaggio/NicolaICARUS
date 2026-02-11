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

//#include "daughtersInfo.h"

using namespace ana;

TFile* file = TFile::Open("/storage/gpfs_data/icarus/local/users/marterop/sbnana_v09_78_06/mc_test/dEdxrestemplates.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

const SpillCut kCRTPMTNeutrino([](const caf::SRSpillProxy* spill){
  for(const auto& match: spill->crtpmt_matches) {
      //Define the interval depending on Data or MC files
    double min_time =-1; double max_time =-1;
//std::cout << " SCAN Run: " << spill->hdr.run << " Event: " << spill->hdr.evt << std::endl;
    if(spill->hdr.evt==122
|| spill->hdr.evt==634
|| spill->hdr.evt==640
|| spill->hdr.evt==778
|| spill->hdr.evt==1658
|| spill->hdr.evt==2058
|| spill->hdr.evt==2061
|| spill->hdr.evt==2083
|| spill->hdr.evt==2147
|| spill->hdr.evt==2440
|| spill->hdr.evt==2851
|| spill->hdr.evt==2880
|| spill->hdr.evt==3712
|| spill->hdr.evt==4297
|| spill->hdr.evt==5061
|| spill->hdr.evt==5261
|| spill->hdr.evt==5405
|| spill->hdr.evt==5770
|| spill->hdr.evt==6278
|| spill->hdr.evt==6544
|| spill->hdr.evt==6641
|| spill->hdr.evt==6712
|| spill->hdr.evt==6734
|| spill->hdr.evt==6807
|| spill->hdr.evt==7417
|| spill->hdr.evt==8048
|| spill->hdr.evt==8222
|| spill->hdr.evt==8301
|| spill->hdr.evt==8560
|| spill->hdr.evt==8621
|| spill->hdr.evt==8802
|| spill->hdr.evt==9307
|| spill->hdr.evt==9426
|| spill->hdr.evt==9576
|| spill->hdr.evt==9711
|| spill->hdr.evt==9799
|| spill->hdr.evt==9868
|| spill->hdr.evt==10253
|| spill->hdr.evt==10962
|| spill->hdr.evt==11546
|| spill->hdr.evt==11809
|| spill->hdr.evt==12148
|| spill->hdr.evt==12348
|| spill->hdr.evt==12426
|| spill->hdr.evt==12531
|| spill->hdr.evt==12588
|| spill->hdr.evt==12819
|| spill->hdr.evt==12858
|| spill->hdr.evt==12879
|| spill->hdr.evt==12929
|| spill->hdr.evt==13085
|| spill->hdr.evt==13679
|| spill->hdr.evt==13744
|| spill->hdr.evt==14562
|| spill->hdr.evt==15096
|| spill->hdr.evt==15180
|| spill->hdr.evt==15231
|| spill->hdr.evt==15724
|| spill->hdr.evt==15814
|| spill->hdr.evt==16111
|| spill->hdr.evt==16138
|| spill->hdr.evt==16273
|| spill->hdr.evt==16298
|| spill->hdr.evt==16867
|| spill->hdr.evt==17158
|| spill->hdr.evt==17270
|| spill->hdr.evt==17344
|| spill->hdr.evt==17438
|| spill->hdr.evt==17905
|| spill->hdr.evt==18106
|| spill->hdr.evt==18191
|| spill->hdr.evt==18598
|| spill->hdr.evt==18823
|| spill->hdr.evt==19097
|| spill->hdr.evt==19112
|| spill->hdr.evt==19286
|| spill->hdr.evt==19510
|| spill->hdr.evt==19567
|| spill->hdr.evt==19632
|| spill->hdr.evt==20274
|| spill->hdr.evt==20708
|| spill->hdr.evt==20781
|| spill->hdr.evt==20908
|| spill->hdr.evt==21034
|| spill->hdr.evt==21056
|| spill->hdr.evt==21128
|| spill->hdr.evt==21229
|| spill->hdr.evt==21455
|| spill->hdr.evt==21672
|| spill->hdr.evt==22008
|| spill->hdr.evt==22569
|| spill->hdr.evt==22638
|| spill->hdr.evt==22674
|| spill->hdr.evt==22704
|| spill->hdr.evt==22751
|| spill->hdr.evt==22884
|| spill->hdr.evt==22983
|| spill->hdr.evt==23121
|| spill->hdr.evt==24279
|| spill->hdr.evt==24514)
{
                //std::cout << " FLAT SCAN Run: " << spill->hdr.run << " Event: " << spill->hdr.evt << std::endl;
  }
    if(spill->hdr.ismc){min_time = 0.0; max_time = 1.6;}
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    if(match.flashGateTime > min_time && match.flashGateTime < max_time && match.flashClassification == 0){  return true;} 



  }
  return false;
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////                                         Event Selection                                                      //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//this macro simply find and store all the informations we need for all the selected tracks, it already computes Eint


std::vector<double> CHIde;

std::vector<double> dEdxind1;
std::vector<double> rrind1;
std::vector<double> dEdxind2;
std::vector<double> rrind2;

std::vector<double> dEdxbestplane;
std::vector<double> rrbestplane;
std::vector<double> pitchbestplane;

std::vector<double> pitch_ind1;
std::vector<double> pitch_ind2;

std::vector<std::vector<double>> pitch_ind1_MU;
std::vector<std::vector<double>> pitch_ind2_MU;

std::vector<std::vector<double>> pitch_ind1_PRO;
std::vector<std::vector<double>> pitch_ind2_PRO;

std::vector<std::vector<double>> pitch_ind1_PI;
std::vector<std::vector<double>> pitch_ind2_PI;

std::vector<std::vector<double>> dEdx_bestplane_MU;
std::vector<std::vector<double>> rr_bestplane_MU;
std::vector<std::vector<double>> dEdx_bestplane_PRO;
std::vector<std::vector<double>> rr_bestplane_PRO;
std::vector<std::vector<double>> dEdx_bestplane_PI;
std::vector<std::vector<double>> rr_bestplane_PI;
std::vector<std::vector<double>> pitch_bestplane_MU;
std::vector<std::vector<double>> pitch_bestplane_PRO;
std::vector<std::vector<double>> pitch_bestplane_PI;

std::vector<double> hitxbestplane;
std::vector<double> hitybestplane;
std::vector<double> hitzbestplane;
std::vector<std::vector<double>> hitx_bestplane_MU;
std::vector<std::vector<double>> hity_bestplane_MU;
std::vector<std::vector<double>> hitz_bestplane_MU;
std::vector<std::vector<double>> hitx_bestplane_PRO;
std::vector<std::vector<double>> hity_bestplane_PRO;
std::vector<std::vector<double>> hitz_bestplane_PRO;
std::vector<std::vector<double>> hitx_bestplane_PI;
std::vector<std::vector<double>> hity_bestplane_PI;
std::vector<std::vector<double>> hitz_bestplane_PI;



std::vector<double> CHIr;
std::vector<double> CHIpitch;
std::vector<int> CHImult;
std::vector<double> CHIEint;
std::vector<double> CHIwidth;
std::vector<double> CHIintegral;
std::vector<double> CHIsumadc;
std::vector<double> CHIx;
std::vector<double> CHIy;
std::vector<double> CHIz;
std::vector<double> CHIdQdx;

std::vector<double> wire;

std::vector<int> CHIbest_plane_mu;
std::vector<int> CHIbest_plane_pro;
std::vector<int> CHIbest_plane_pi;

std::vector<std::vector<double>> chi_quadro_coll_MU;
std::vector<std::vector<double>> chi_quadro_ind1_MU;
std::vector<std::vector<double>> chi_quadro_ind2_MU;

std::vector<std::vector<double>> chi_quadro_coll_PRO;
std::vector<std::vector<double>> chi_quadro_ind1_PRO;
std::vector<std::vector<double>> chi_quadro_ind2_PRO;

std::vector<std::vector<double>> chi_quadro_coll_PI;
std::vector<std::vector<double>> chi_quadro_ind1_PI;
std::vector<std::vector<double>> chi_quadro_ind2_PI;



std::vector<std::vector<double>> dEdx_ind1_MU;
std::vector<std::vector<double>> rr_ind1_MU;
std::vector<std::vector<double>> dEdx_ind2_MU;
std::vector<std::vector<double>> rr_ind2_MU;

std::vector<std::vector<double>> dEdx_ind1_PRO;
std::vector<std::vector<double>> rr_ind1_PRO;
std::vector<std::vector<double>> dEdx_ind2_PRO;
std::vector<std::vector<double>> rr_ind2_PRO;

std::vector<std::vector<double>> dEdx_ind1_PI;
std::vector<std::vector<double>> rr_ind1_PI;
std::vector<std::vector<double>> dEdx_ind2_PI;
std::vector<std::vector<double>> rr_ind2_PI;

std::vector<std::vector<double>> CHIwireMU;
std::vector<std::vector<double>> CHIdeMU;
std::vector<std::vector<double>> CHIrrMUninv;
std::vector<std::vector<double>> CHIEintMUallt;
std::vector<std::vector<double>> CHIrrMU;
std::vector<std::vector<int>> CHImultMU;
std::vector<std::vector<double>> CHIwidthMU;
std::vector<std::vector<double>> CHIpitchMU;
std::vector<std::vector<double>> CHIdQdxMU;

std::vector<std::vector<double>> CHIdePRO;
std::vector<std::vector<double>> CHIrrPROninv;
std::vector<std::vector<double>> CHIEintPROallt;
std::vector<std::vector<double>> CHIrrPRO;
std::vector<std::vector<int>> CHImultPRO;
std::vector<std::vector<double>> CHIwidthPRO;
std::vector<std::vector<double>> CHIpitchPRO;
std::vector<std::vector<double>> CHIdQdxPRO;

std::vector<std::vector<double>> CHIdePI;
std::vector<std::vector<double>> CHIrrPIninv;
std::vector<std::vector<double>> CHIpitchPI;

std::vector<std::vector<double>> CHIintegralMU;
std::vector<std::vector<double>> CHIintegralPRO;
std::vector<std::vector<double>> CHIsumadcMU;
std::vector<std::vector<double>> CHIsumadcPRO;

std::vector<std::vector<double>> CHIxMU;
std::vector<std::vector<double>> CHIyMU;
std::vector<std::vector<double>> CHIzMU;
std::vector<std::vector<double>> CHIxPRO;
std::vector<std::vector<double>> CHIyPRO;
std::vector<std::vector<double>> CHIzPRO;

std::vector<std::vector<double>> CHIxPI;
std::vector<std::vector<double>> CHIyPI;
std::vector<std::vector<double>> CHIzPI;

std::vector<std::vector<double>> CHIstart3DrecoMU;
std::vector<std::vector<double>> CHIstart3DmcMU;
std::vector<std::vector<double>> CHIstart3DrecoPRO;
std::vector<std::vector<double>> CHIstart3DmcPRO;
std::vector<std::vector<double>> CHIstart3DrecoPI;
std::vector<std::vector<double>> CHIstart3DmcPI;

std::vector<std::vector<double>> CHIend3DrecoMU;
std::vector<std::vector<double>> CHIend3DmcMU;
std::vector<std::vector<double>> CHIend3DrecoPRO;
std::vector<std::vector<double>> CHIend3DmcPRO;
std::vector<std::vector<double>> CHIend3DrecoPI;
std::vector<std::vector<double>> CHIend3DmcPI;

std::vector<std::vector<double>> CHIgenMomentumMU;
std::vector<std::vector<double>> CHIgenMomentumPRO;


//std::vector<std::vector<unsigned int>> MUdaughtersID;
//std::vector<std::vector<unsigned int>> PROdaughtersID;

std::vector<unsigned int> PROndaughters_reco;
std::vector<unsigned int> MUndaughters_reco;

std::vector<double> CHItlRECOmu;
std::vector<double> CHItlRECOpro;
std::vector<double> CHItlMCmu;
std::vector<double> CHItlMCpro;

std::vector<double> CHItlRECOpi;
std::vector<double> CHItlMCpi;

std::vector<std::vector<double>> CHIvtxRECOmu;
std::vector<std::vector<double>> CHIvtxMCmu;
std::vector<std::vector<double>> CHIvtxRECOpro;
std::vector<std::vector<double>> CHIvtxMCpro;
std::vector<std::vector<double>> CHIvtxRECOpi;
std::vector<std::vector<double>> CHIvtxMCpi;

std::vector<double> vdirx_mu;
std::vector<double> vdiry_mu;
std::vector<double> vdirz_mu;
std::vector<double> vdirx_pro;
std::vector<double> vdiry_pro;
std::vector<double> vdirz_pro;

std::vector<unsigned int> vrun_mu;
std::vector<unsigned int> vsubrun_mu;
std::vector<unsigned int> vevt_mu;
std::vector<bool> visMC_mu;

std::vector<unsigned int> vrun_pro;
std::vector<unsigned int> vsubrun_pro;
std::vector<unsigned int> vevt_pro;
std::vector<bool> visMC_pro;

std::vector<unsigned int> vrun_pi;
std::vector<unsigned int> vsubrun_pi;
std::vector<unsigned int> vevt_pi;
std::vector<bool> visMC_pi;

std::vector<int> end_process_mu;
std::vector<int> end_process_pro;
std::vector<int> end_process_pi;

std::vector<double> dummy;
//std::vector<int> dummy_id;

//std::vector<std::vector<daughtersInfo>> muons_daughter;
//std::vector<std::vector<daughtersInfo>> protons_daughter;

std::vector<std::vector<double>> rm_mu_coll;
std::vector<std::vector<double>> rr_rm_mu_coll;
std::vector<std::vector<double>> rm_pro_coll;
std::vector<std::vector<double>> rr_rm_pro_coll;
std::vector<std::vector<double>> rm_pi_coll;
std::vector<std::vector<double>> rr_rm_pi_coll;

std::vector<std::vector<double>> rm_mu_bestplane;
std::vector<std::vector<double>> rr_rm_mu_bestplane;
std::vector<std::vector<double>> rm_pro_bestplane;
std::vector<std::vector<double>> rr_rm_pro_bestplane;
std::vector<std::vector<double>> rm_pi_bestplane;
std::vector<std::vector<double>> rr_rm_pi_bestplane;

std::vector<std::vector<double>> rm_mu_ind1;
std::vector<std::vector<double>> rr_rm_mu_ind1;
std::vector<std::vector<double>> rm_pro_ind1;
std::vector<std::vector<double>> rr_rm_pro_ind1;
std::vector<std::vector<double>> rm_pi_ind1;
std::vector<std::vector<double>> rr_rm_pi_ind1;

std::vector<std::vector<double>> rm_mu_ind2;
std::vector<std::vector<double>> rr_rm_mu_ind2;
std::vector<std::vector<double>> rm_pro_ind2;
std::vector<std::vector<double>> rr_rm_pro_ind2;
std::vector<std::vector<double>> rm_pi_ind2;
std::vector<std::vector<double>> rr_rm_pi_ind2;

std::vector<bool> has_true_secondaries_MU;
std::vector<bool> has_true_secondaries_PRO;
std::vector<bool> has_true_secondaries_PI;

std::vector<double> visible_energy_mu;
std::vector<double> visible_energy_pro;
std::vector<double> visible_energy_pi;

std::vector<double> hit_purity_mu;
std::vector<double> hit_completeness_mu;
std::vector<double> hit_purity_pro;
std::vector<double> hit_completeness_pro;
std::vector<double> hit_purity_pi;
std::vector<double> hit_completeness_pi;

std::vector<double> energy_purity_mu;
std::vector<double> energy_completeness_mu;
std::vector<double> energy_purity_pro;
std::vector<double> energy_completeness_pro;
std::vector<double> energy_purity_pi;
std::vector<double> energy_completeness_pi;

std::vector<std::vector<int>> pdg_matches_mu;
std::vector<std::vector<double>> energy_matches_mu;
std::vector<std::vector<double>> end_points_distance_matches_mu;
std::vector<int> start_process_mu;
std::vector<double> energy_at_last_hit_mu;
std::vector<double> energy_at_first_hit_mu;

std::vector<std::vector<int>> pdg_matches_pro;
std::vector<std::vector<double>> energy_matches_pro;
std::vector<std::vector<double>> end_points_distance_matches_pro;
std::vector<int> start_process_pro;
std::vector<double> energy_at_last_hit_pro;
std::vector<double> energy_at_first_hit_pro;

std::vector<std::vector<int>> pdg_matches_pi;
std::vector<std::vector<double>> energy_matches_pi;
std::vector<std::vector<double>> end_points_distance_matches_pi;
std::vector<int> start_process_pi;
std::vector<double> energy_at_last_hit_pi;
std::vector<double> energy_at_first_hit_pi;

std::vector<std::vector<bool>> is_daughter_mu;
std::vector<std::vector<bool>> is_daughter_pro;
std::vector<std::vector<bool>> is_daughter_pi;

ofstream dumpPionSecondaries("dumpPionSecondaries.txt");
ofstream dumpMuonSecondaries("dumpMuonSecondaries.txt");
ofstream dumpProtonSecondaries("dumpProtonSecondaries.txt");

bool is_mc=true; 

int choosen_plane=2;

const SpillMultiVar DataLoader([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
  //TF1 *correction_function = new TF1("correction_function",correction_function_25only, 0., 16., 6 );
  //correction_function->SetParameters(-1.08107e+00, -1.38302e+00, -2.84635e-02, 5.87297e-02, 7.17814e+00, -3.51461e+00);

  std::vector<double> vector_active;

  int count_slices=-1;
  for (auto const& islc : sr->slc)
  {     
    count_slices+=1;

    int ipfp_mu=-1;
    std::vector<int> v_ipfp_pro; 
    std::vector<int> v_ipfp_pi; 
    double distanza;

      /*
      if(automatic_selection_1muNp_truth(sr,islc,10,100,choosen_plane))
	    {
        ipfp_mu=find_muon(islc,10,choosen_plane); //trovo l'indice del muone tra tutte le particelle nella slice

        if(ipfp_mu==-1)continue;
          
        if(is_mc) // we check for MC info only if we work on MC data 
        {
	        if(!(islc.truth.index>=0 && (classification_type(sr,islc)==2 || classification_type(sr,islc)==5 ))) //1mu1p E>50 true
	        {
            continue;
          }
        }
      */

    //MUONS
    ipfp_mu=find_truth_muon(islc,10,choosen_plane);
    //rolling median muons

    //searching all the matched true particles
    if(is_mc && ipfp_mu!=-1)
    {

      cout << islc.reco.pfp[ipfp_mu].trk.truth.p.endE << endl;
      energy_at_first_hit_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.startE);
      energy_at_last_hit_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.endE);
      start_process_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.start_process);

      std::vector<int> temp_pdg;
      std::vector<double> temp_energy;
      std::vector<bool> is_daughter;
      std::vector<double> end_points_distance;
      std::vector<double> is_primary;
      for (const auto& true_p : sr->true_particles)
      {
        for (auto const& match: islc.reco.pfp[ipfp_mu].trk.truth.matches)
        {
          if(true_p.G4ID==match.G4ID)
          {
            temp_pdg.push_back(true_p.pdg);
            temp_energy.push_back(match.energy/3.);

            TVector3 end_point_match;
            end_point_match.SetXYZ(true_p.end.x, true_p.end.y, true_p.end.z);
            TVector3 end_point_track;
            end_point_track.SetXYZ(islc.reco.pfp[ipfp_mu].trk.end.x, islc.reco.pfp[ipfp_mu].trk.end.y, islc.reco.pfp[ipfp_mu].trk.end.z);
            end_points_distance.push_back( (end_point_match-end_point_track).Mag() );

            if((int)true_p.parent == islc.reco.pfp[ipfp_mu].trk.truth.p.G4ID)is_daughter.push_back(true);
            else{is_daughter.push_back(false);}
          }
        }
      }
      pdg_matches_mu.push_back(temp_pdg);
      end_points_distance_matches_mu.push_back(end_points_distance);
      energy_matches_mu.push_back(temp_energy);
      is_daughter_mu.push_back(is_daughter);
      temp_pdg.clear();
      end_points_distance.clear();
      temp_energy.clear();
      is_daughter.clear();
    }


    if(ipfp_mu!=-1)
    {
      int plane_for_rm=-1;
      int bestplane=islc.reco.pfp[ipfp_mu].trk.bestplane;
      if(bestplane!=-1 && islc.reco.pfp[ipfp_mu].trk.calo[bestplane].nhit>0){plane_for_rm=bestplane;}
      else{plane_for_rm=2;}

      rr_rm_mu_bestplane.push_back(rollingMedian(islc,ipfp_mu,plane_for_rm)[0]);
      rr_rm_mu_ind1.push_back(rollingMedian(islc,ipfp_mu,0)[0]);
      rr_rm_mu_ind2.push_back(rollingMedian(islc,ipfp_mu,1)[0]);
      rr_rm_mu_coll.push_back(rollingMedian(islc,ipfp_mu,2)[0]);
      rm_mu_bestplane.push_back(rollingMedian(islc,ipfp_mu,plane_for_rm)[1]);
      rm_mu_ind1.push_back(rollingMedian(islc,ipfp_mu,0)[1]);
      rm_mu_ind2.push_back(rollingMedian(islc,ipfp_mu,1)[1]);
      rm_mu_coll.push_back(rollingMedian(islc,ipfp_mu,2)[1]);
    }

    //CHI2 MUONS
    if(int(ipfp_mu)!=-1)
    {
      std::vector<double> output;
      std::vector<double> dedx;
      std::vector<double> rr;

      std::array<std::vector<double>,3> chi_diff_planes_mu;
      //storing the chi2 of muon in the three wire planes
      for(int iplane=0; iplane<3; iplane++)
      {
        std::vector<double> output;
        std::vector<double> dedx;
        std::vector<double> rr;
        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[iplane].points.size(); ++ihit )
        {
          dedx.push_back(islc.reco.pfp[ipfp_mu].trk.calo[iplane].points[ihit].dedx);
          rr.push_back(islc.reco.pfp[ipfp_mu].trk.calo[iplane].points[ihit].rr);
        } // calo points
        output = chi2_ALG(dedx,rr,0.0,25.0);
        chi_diff_planes_mu[iplane]=output;
      }
      chi_quadro_ind1_MU.push_back(chi_diff_planes_mu[0]);
      chi_quadro_ind2_MU.push_back(chi_diff_planes_mu[1]);
      chi_quadro_coll_MU.push_back(chi_diff_planes_mu[2]);
    

    //RUN and EVT information
    vrun_mu.push_back(sr->hdr.run);
    vsubrun_mu.push_back(sr->hdr.subrun);
    vevt_mu.push_back(sr->hdr.evt);
    visMC_mu.push_back(sr->hdr.ismc);
    }

	  for (std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
		{
      if(int(ipfp)==ipfp_mu)continue;

      int use_plane = -1;
      if(choosen_plane==-1)
      {
        if(islc.reco.pfp[ipfp].trk.bestplane==0 || islc.reco.pfp[ipfp].trk.bestplane==1 || islc.reco.pfp[ipfp].trk.bestplane==2){use_plane=islc.reco.pfp[ipfp].trk.bestplane;}
        else{use_plane=2;}
      }
      else{use_plane=choosen_plane;}

      //PROTONS
      if(id_pfp_truth(islc, ipfp, 10,use_plane)==1)
			{
        int ipfp_pro=int(ipfp);

        //searching all the matched true particles
        if(is_mc)
        {
          energy_at_first_hit_pro.push_back(islc.reco.pfp[ipfp_pro].trk.truth.p.startE);
          energy_at_last_hit_pro.push_back(islc.reco.pfp[ipfp_pro].trk.truth.p.endE);
          start_process_pro.push_back(islc.reco.pfp[ipfp_pro].trk.truth.p.start_process);

          std::vector<int> temp_pdg;
          std::vector<double> temp_energy;
          std::vector<bool> is_daughter;
          std::vector<double> end_points_distance;
          for (const auto& true_p : sr->true_particles)
          {
            for (auto const& match: islc.reco.pfp[ipfp_pro].trk.truth.matches)
            {
              if(true_p.G4ID==match.G4ID)
              {
                temp_pdg.push_back(true_p.pdg);

                TVector3 end_point_match;
                end_point_match.SetXYZ(true_p.end.x, true_p.end.y, true_p.end.z);
                TVector3 end_point_track;
                end_point_track.SetXYZ(islc.reco.pfp[ipfp_pro].trk.end.x, islc.reco.pfp[ipfp_pro].trk.end.y, islc.reco.pfp[ipfp_pro].trk.end.z);
                end_points_distance.push_back( (end_point_match-end_point_track).Mag() );

                temp_energy.push_back(match.energy/3.);
                if((int)true_p.parent == islc.reco.pfp[ipfp_pro].trk.truth.p.G4ID)is_daughter.push_back(true);
                else{is_daughter.push_back(false);}
              }
            }
          } 
          pdg_matches_pro.push_back(temp_pdg);
          end_points_distance_matches_pro.push_back(end_points_distance);
          energy_matches_pro.push_back(temp_energy);
          is_daughter_pro.push_back(is_daughter);
          temp_pdg.clear();
          end_points_distance.clear();
          temp_energy.clear();
          is_daughter.clear();
        }

        //rolling median protons

        int plane_for_rm=-1;
        int bestplane=islc.reco.pfp[ipfp_pro].trk.bestplane;
        if(bestplane!=-1 && islc.reco.pfp[ipfp_pro].trk.calo[bestplane].nhit>0){plane_for_rm=bestplane;}
        else{plane_for_rm=2;}

        rr_rm_pro_bestplane.push_back(rollingMedian(islc,ipfp_pro,plane_for_rm)[0]);
        rr_rm_pro_ind1.push_back(rollingMedian(islc,ipfp_pro,0)[0]);
        rr_rm_pro_ind2.push_back(rollingMedian(islc,ipfp_pro,1)[0]);
        rr_rm_pro_coll.push_back(rollingMedian(islc,ipfp_pro,2)[0]);
        rm_pro_bestplane.push_back(rollingMedian(islc,ipfp_pro,plane_for_rm)[1]);
        rm_pro_ind1.push_back(rollingMedian(islc,ipfp_pro,0)[1]);
        rm_pro_ind2.push_back(rollingMedian(islc,ipfp_pro,1)[1]);
        rm_pro_coll.push_back(rollingMedian(islc,ipfp_pro,2)[1]);

        //CHI2 PROTONS
        std::array<std::vector<double>,3> chi_diff_planes;          
        for(int iplane=0; iplane<3; iplane++)
        {
          std::vector<double> output;
          std::vector<double> dedx;
          std::vector<double> rr;
          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_pro].trk.calo[iplane].points.size(); ++ihit )
          {
            dedx.push_back(islc.reco.pfp[ipfp_pro].trk.calo[iplane].points[ihit].dedx);
            rr.push_back(islc.reco.pfp[ipfp_pro].trk.calo[iplane].points[ihit].rr);
          } // calo points
          output = chi2_ALG(dedx,rr,0.0,25.0);
          chi_diff_planes[iplane]=output;
        }
        chi_quadro_ind1_PRO.push_back(chi_diff_planes[0]);
        chi_quadro_ind2_PRO.push_back(chi_diff_planes[1]);
        chi_quadro_coll_PRO.push_back(chi_diff_planes[2]);

        //RUN and EVT informations
        v_ipfp_pro.push_back(ipfp_pro);
        vrun_pro.push_back(sr->hdr.run);
        vsubrun_pro.push_back(sr->hdr.subrun);
        vevt_pro.push_back(sr->hdr.evt);
        visMC_pro.push_back(sr->hdr.ismc);
      }//if it is a proton

      //PIONS
      if(id_pfp_truth(islc, ipfp, 10,use_plane)==2)
      {

        //searching all the matched true particles
        if(is_mc)
        {
          energy_at_first_hit_pi.push_back(islc.reco.pfp[ipfp].trk.truth.p.startE);
          energy_at_last_hit_pi.push_back(islc.reco.pfp[ipfp].trk.truth.p.endE);
          start_process_pi.push_back(islc.reco.pfp[ipfp].trk.truth.p.start_process);

          std::vector<int> temp_pdg;
          std::vector<double> temp_energy;
          std::vector<bool> is_daughter;
          std::vector<double> end_points_distance;
          for (const auto& true_p : sr->true_particles)
          {
            for (auto const& match: islc.reco.pfp[ipfp].trk.truth.matches)
            {
              if(true_p.G4ID==match.G4ID)
              {
                temp_pdg.push_back(true_p.pdg);

                TVector3 end_point_match;
                end_point_match.SetXYZ(true_p.end.x, true_p.end.y, true_p.end.z);
                TVector3 end_point_track;
                end_point_track.SetXYZ(islc.reco.pfp[ipfp].trk.end.x, islc.reco.pfp[ipfp].trk.end.y, islc.reco.pfp[ipfp].trk.end.z);
                end_points_distance.push_back( (end_point_match-end_point_track).Mag() );

                temp_energy.push_back(match.energy/3.);
                if((int)true_p.parent == islc.reco.pfp[ipfp].trk.truth.p.G4ID)is_daughter.push_back(true);
                else{is_daughter.push_back(false);}
              }
            }
          } 
          pdg_matches_pi.push_back(temp_pdg);
          end_points_distance_matches_pi.push_back(end_points_distance);
          energy_matches_pi.push_back(temp_energy);
          is_daughter_pi.push_back(is_daughter);
          temp_pdg.clear();
          end_points_distance.clear();
          temp_energy.clear();
          is_daughter.clear();
        }

      int plane_for_rm=-1;
      int bestplane=islc.reco.pfp[ipfp].trk.bestplane;
      if(bestplane!=-1 && islc.reco.pfp[ipfp].trk.calo[bestplane].nhit>0){plane_for_rm=bestplane;}
      else{plane_for_rm=2;}

      //rolling median pions
      rr_rm_pi_bestplane.push_back(rollingMedian(islc,ipfp,plane_for_rm)[0]);
      rr_rm_pi_ind1.push_back(rollingMedian(islc,ipfp,0)[0]);
      rr_rm_pi_ind2.push_back(rollingMedian(islc,ipfp,1)[0]);
      rr_rm_pi_coll.push_back(rollingMedian(islc,ipfp,2)[0]);
      rm_pi_bestplane.push_back(rollingMedian(islc,ipfp,plane_for_rm)[1]);
      rm_pi_ind1.push_back(rollingMedian(islc,ipfp,0)[1]);
      rm_pi_ind2.push_back(rollingMedian(islc,ipfp,1)[1]);
      rm_pi_coll.push_back(rollingMedian(islc,ipfp,2)[1]);
        

        std::array<std::vector<double>,3> chi_diff_planes;
        //storing the chi2 of the pion in the three wire planes
        for(int iplane=0; iplane<3; iplane++)
        {
          std::vector<double> output;
          std::vector<double> dedx;
          std::vector<double> rr;
          for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[iplane].points.size(); ++ihit )
          {
            dedx.push_back(islc.reco.pfp[ipfp].trk.calo[iplane].points[ihit].dedx);
            rr.push_back(islc.reco.pfp[ipfp].trk.calo[iplane].points[ihit].rr);
          } // calo points
          output = chi2_ALG(dedx,rr,0.0,25.0);
          chi_diff_planes[iplane]=output;
        }
        chi_quadro_ind1_PI.push_back(chi_diff_planes[0]);
        chi_quadro_ind2_PI.push_back(chi_diff_planes[1]);
        chi_quadro_coll_PI.push_back(chi_diff_planes[2]);

        //RUN and EVT informations
        v_ipfp_pi.push_back(ipfp);
        vrun_pi.push_back(sr->hdr.run);
        vsubrun_pi.push_back(sr->hdr.subrun);
        vevt_pi.push_back(sr->hdr.evt);
        visMC_pi.push_back(sr->hdr.ismc);
    }
  }//loop over all particles in the slice

          /*
          //fill secondaries dEdx and RR of MUON and PROTONS
          //MUONS
          std::vector<daughtersInfo> dummy_muons_daughtersInfo;
          bool has_daughters = false;
          for(int daughter=0; daughter<int(islc.reco.pfp[ipfp_mu].ndaughters); daughter++)
          {
            has_daughters=true;
            daughtersInfo dummy_daughter;
            int daughters_id = islc.reco.pfp[ipfp_mu].daughters[daughter];
            for(std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp)
            {
              std::vector<double> dEdxdummy;
              std::vector<double> rrdummy;
              if(islc.reco.pfp[ipfp].id==daughters_id)
              {
                for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
                {
                  if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx<0.5 || islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx>500.)continue;

                  dEdxdummy.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
                  rrdummy.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
                }
                if(is_mc){dummy_daughter.pdg=islc.reco.pfp[ipfp].trk.truth.p.pdg;}
                dummy_daughter.dEdx=dEdxdummy;
                dummy_daughter.rr=rrdummy;
              }
            }
            dummy_muons_daughtersInfo.push_back(dummy_daughter);
          }
          if(!has_daughters)
          {
            daughtersInfo dummy_daughter;
            dummy_daughter.dEdx={-1.};
            dummy_daughter.rr={-1.};
            dummy_daughter.pdg={-1};
            dummy_muons_daughtersInfo.push_back(dummy_daughter);
          }
          muons_daughter.push_back(dummy_muons_daughtersInfo);
          

          //for PROTONS
          for(auto const &ipfp_pro : v_ipfp_pro)
          {
            has_daughters=false;
            std::vector<daughtersInfo> dummy_protons_daughtersInfo;
            for(int daughter=0; daughter<int(islc.reco.pfp[ipfp_pro].ndaughters); daughter++)
            {
              has_daughters=true;
              daughtersInfo dummy_daughter;
              int daughters_id = islc.reco.pfp[ipfp_pro].daughters[daughter];
              for(std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp)
              {
                std::vector<double> dEdxdummy;
                std::vector<double> rrdummy;
                if(islc.reco.pfp[ipfp].id==daughters_id)
                {
                  for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit )
                  {
                    if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx<0.5 || islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx>500.)continue;

                    dEdxdummy.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
                    rrdummy.push_back(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
                  }
                  if(is_mc){dummy_daughter.pdg=islc.reco.pfp[ipfp].trk.truth.p.pdg;}
                  dummy_daughter.dEdx=dEdxdummy;
                  dummy_daughter.rr=rrdummy;
                }
              } 
              dummy_protons_daughtersInfo.push_back(dummy_daughter);
            }
            if(!has_daughters)
            {
              daughtersInfo dummy_daughter;
              dummy_daughter.dEdx={-1.};
              dummy_daughter.rr={-1.};
              dummy_daughter.pdg={-1};
              dummy_protons_daughtersInfo.push_back(dummy_daughter);
            } 
            protons_daughter.push_back(dummy_protons_daughtersInfo);

          }
            */

  //MC informations storage
  if(is_mc)
  {
    //for(int id=0; id<int(islc.reco.pfp[ipfp_mu].trk.truth.p.daughters.size()); id++)
    //{ dummy_id.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.daughters[id]); //}
    //MUdaughtersID.push_back(dummy_id);
    //dummy_id.clear();

    //MUON TRACK TRUE INFO
    if(ipfp_mu!=-1)
    {
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.genp.x);
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.genp.y);
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.genp.z);
      CHIgenMomentumMU.push_back(dummy);
      dummy.clear();

      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.start.x);
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.start.y);
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.start.z);
      CHIstart3DmcMU.push_back(dummy);
      dummy.clear();
            
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.end.x);
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.end.y);
      dummy.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.end.z);
      CHIend3DmcMU.push_back(dummy);
      dummy.clear();

      if(islc.reco.pfp[ipfp_mu].trk.truth.p.cryostat!=-1)visible_energy_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.plane[islc.reco.pfp[ipfp_mu].trk.truth.p.cryostat][2].visE);

      CHItlMCmu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.length);

      end_process_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.end_process);

      hit_purity_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.bestmatch.hit_purity);
      hit_completeness_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.bestmatch.hit_completeness);
      energy_purity_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.bestmatch.energy_purity);
      energy_completeness_mu.push_back(islc.reco.pfp[ipfp_mu].trk.truth.bestmatch.energy_completeness);

      dumpMuonSecondaries << "----------------------------------------------------------" << endl;
      dumpMuonSecondaries << "RUN= " << sr->hdr.run << " EVT= " << sr->hdr.evt << " | " << "particleEP= " << islc.reco.pfp[ipfp_mu].trk.truth.p.end_process << " | secondaries: " ;
      bool has_true_secondaries = false;
      int G4ID = islc.reco.pfp[ipfp_mu].trk.truth.p.G4ID;
      for(const auto& true_particle : sr->true_particles)
      {
        if((int)true_particle.parent==G4ID)
        {
          has_true_secondaries=true;
          dumpMuonSecondaries << true_particle.pdg << " ";
        }
      }
      dumpMuonSecondaries << endl;
      dumpMuonSecondaries << "INDUCTION 1 ";
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[0].points.size(); ++ihit){dumpMuonSecondaries << islc.reco.pfp[ipfp_mu].trk.calo[0].points[ihit].x << " ";}
      dumpMuonSecondaries << endl;
      dumpMuonSecondaries << "INDUCTION 2 ";
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[1].points.size(); ++ihit){dumpMuonSecondaries << islc.reco.pfp[ipfp_mu].trk.calo[1].points[ihit].x << " ";}
      dumpMuonSecondaries << endl;
      dumpMuonSecondaries << "COLLECTION "; 
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit){dumpMuonSecondaries << islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x << " ";}
      dumpMuonSecondaries << endl;
      has_true_secondaries_MU.push_back(has_true_secondaries);
      //searching all the matched true particles
      for (const auto& true_p : sr->true_particles)
      {
        for (auto const& match: islc.reco.pfp[ipfp_mu].trk.truth.matches)
        {
          if(true_p.G4ID==match.G4ID)
          {
            dumpMuonSecondaries << "matched_energy= " <<match.energy/3 << " pdg= " << true_p.pdg << " | " << endl; 
          }
        }
      }
      dumpMuonSecondaries << endl;

      // VERTEX TRUE INFO            
      dummy.push_back(islc.truth.position.x);
      dummy.push_back(islc.truth.position.y);
      dummy.push_back(islc.truth.position.z);
      CHIvtxMCmu.push_back(dummy);
      dummy.clear();
    }
    
    //PROTONS TRACK TRUE INFO
    for(int i=0; i<int(v_ipfp_pro.size()); i++)
    {
      int ind_ipfp=v_ipfp_pro[i];

      //for(int id=0; id<int(islc.reco.pfp[ind_ipfp].trk.truth.p.daughters.size()); id++)
      //{ dummy_id.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.daughters[id]);}
      //PROdaughtersID.push_back(dummy_id);
      //dummy_id.clear();

      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.genp.x);
      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.genp.y);
      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.genp.z);
      CHIgenMomentumPRO.push_back(dummy);
      dummy.clear();  

      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.start.x);
      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.start.y);
      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.start.z);
      CHIstart3DmcPRO.push_back(dummy);
      dummy.clear();

      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.end.x);
      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.end.y);
      dummy.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.end.z);
      CHIend3DmcPRO.push_back(dummy);
      dummy.clear();

      if(islc.reco.pfp[ind_ipfp].trk.truth.p.cryostat!=-1)visible_energy_pro.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.plane[islc.reco.pfp[ind_ipfp].trk.truth.p.cryostat][2].visE);

      CHItlMCpro.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.length);

      end_process_pro.push_back(islc.reco.pfp[ind_ipfp].trk.truth.p.end_process);

      hit_purity_pro.push_back(islc.reco.pfp[ind_ipfp].trk.truth.bestmatch.hit_purity);
      hit_completeness_pro.push_back(islc.reco.pfp[ind_ipfp].trk.truth.bestmatch.hit_completeness);
      energy_purity_pro.push_back(islc.reco.pfp[ind_ipfp].trk.truth.bestmatch.energy_purity);
      energy_completeness_pro.push_back(islc.reco.pfp[ind_ipfp].trk.truth.bestmatch.energy_completeness);

      dumpProtonSecondaries << "----------------------------------------------------------" << endl;
      dumpProtonSecondaries << "RUN= " << sr->hdr.run << " EVT= " << sr->hdr.evt << " | " << "particleEP= " << islc.reco.pfp[ind_ipfp].trk.truth.p.end_process << " | secondaries: " ;
      
      bool has_true_secondaries = false;
      int G4ID = islc.reco.pfp[ind_ipfp].trk.truth.p.G4ID;
      for(const auto& true_particle : sr->true_particles)
      {
        if((int)true_particle.parent==G4ID)
        {
          has_true_secondaries=true;
          dumpProtonSecondaries << true_particle.pdg << " " << true_particle.start.x << " " << true_particle.end.x << endl;
        }
      }
      has_true_secondaries_PRO.push_back(has_true_secondaries);

      dumpProtonSecondaries << endl;
      dumpProtonSecondaries << "INDUCTION 1 ";
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ind_ipfp].trk.calo[0].points.size(); ++ihit){dumpProtonSecondaries << islc.reco.pfp[ind_ipfp].trk.calo[0].points[ihit].x << " ";}
      dumpProtonSecondaries << endl;
      dumpProtonSecondaries << "INDUCTION 2 ";
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ind_ipfp].trk.calo[1].points.size(); ++ihit){dumpProtonSecondaries << islc.reco.pfp[ind_ipfp].trk.calo[1].points[ihit].x << " ";}
      dumpProtonSecondaries << endl;
      dumpProtonSecondaries << "COLLECTION "; 
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ind_ipfp].trk.calo[2].points.size(); ++ihit){dumpProtonSecondaries << islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].x << " ";}
      dumpProtonSecondaries << endl;
      
      //searching all the matched true particles
      for (const auto& true_p : sr->true_particles)
      {
        for (auto const& match: islc.reco.pfp[ind_ipfp].trk.truth.matches)
        {
          if(true_p.G4ID==match.G4ID)
          {
            dumpProtonSecondaries << "matched_energy= " <<match.energy/3 << " pdg= " << true_p.pdg << " | " << endl; 
          }
        }
      }
      dumpProtonSecondaries << endl;

      //HAS SECONDARIES AT TRUTH LEVEL
      //bool has_true_secondaries = false;
      //int G4ID = islc.reco.pfp[ind_ipfp].trk.truth.p.G4ID;
      //for(const auto& true_particle : sr->true_particles)
      //{if((int)true_particle.parent==G4ID)has_true_secondaries=true;}
      //has_true_secondaries_PRO.push_back(has_true_secondaries);

      // VERTEX TRUE INFO
      dummy.push_back(islc.truth.position.x);
      dummy.push_back(islc.truth.position.y);
      dummy.push_back(islc.truth.position.z);
      CHIvtxMCpro.push_back(dummy);
      dummy.clear();
    }

    //PIONS TRACK TRUE INFO
    for(int i=0; i<int(v_ipfp_pi.size()); i++)
    {
      int pfp_index_pi=v_ipfp_pi[i];

      dummy.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.start.x);
      dummy.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.start.y);
      dummy.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.start.z);
      CHIstart3DmcPI.push_back(dummy);
      dummy.clear();

      dummy.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.end.x);
      dummy.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.end.y);
      dummy.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.end.z);
      CHIend3DmcPI.push_back(dummy);
      dummy.clear();

      if(islc.reco.pfp[pfp_index_pi].trk.truth.p.cryostat!=-1)visible_energy_pi.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.plane[islc.reco.pfp[pfp_index_pi].trk.truth.p.cryostat][2].visE);

      CHItlMCpi.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.length);

      end_process_pi.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.p.end_process);

      hit_purity_pi.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.bestmatch.hit_purity);
      hit_completeness_pi.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.bestmatch.hit_completeness);
      energy_purity_pi.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.bestmatch.energy_purity);
      energy_completeness_pi.push_back(islc.reco.pfp[pfp_index_pi].trk.truth.bestmatch.energy_completeness);

      //HAS SECONDARIES AT TRUTH LEVEL
      //endprocess codes: 
      //3 : DECAY
      //9 : pipInelastic
      //10 : pimInelastic
      //45 : BertiniCaptureAtRest
      dumpPionSecondaries << "----------------------------------------------------------" << endl;
      dumpPionSecondaries << "RUN= " << sr->hdr.run << " EVT= " << sr->hdr.evt << " | " << "particleEP= " << islc.reco.pfp[pfp_index_pi].trk.truth.p.end_process << " | secondaries: " ;
      bool has_true_secondaries = false;
      int G4ID = islc.reco.pfp[pfp_index_pi].trk.truth.p.G4ID;
      for(const auto& true_particle : sr->true_particles)
      {
        if((int)true_particle.parent==G4ID)
        {
          has_true_secondaries=true;
          dumpPionSecondaries << true_particle.pdg << " ";
        }
      }
      dumpPionSecondaries << endl;
      dumpPionSecondaries << "INDUCTION 1 ";
      for(std::size_t ihit(0); ihit < islc.reco.pfp[pfp_index_pi].trk.calo[0].points.size(); ++ihit){dumpPionSecondaries << islc.reco.pfp[pfp_index_pi].trk.calo[0].points[ihit].x << " ";}
      dumpPionSecondaries << endl;
      dumpPionSecondaries << "INDUCTION 2 ";
      for(std::size_t ihit(0); ihit < islc.reco.pfp[pfp_index_pi].trk.calo[1].points.size(); ++ihit){dumpPionSecondaries << islc.reco.pfp[pfp_index_pi].trk.calo[1].points[ihit].x << " ";}
      dumpPionSecondaries << endl;
      dumpPionSecondaries << "COLLECTION "; 
      for(std::size_t ihit(0); ihit < islc.reco.pfp[pfp_index_pi].trk.calo[2].points.size(); ++ihit){dumpPionSecondaries << islc.reco.pfp[pfp_index_pi].trk.calo[2].points[ihit].x << " ";}
      dumpPionSecondaries << endl;
      has_true_secondaries_PI.push_back(has_true_secondaries);

      // VERTEX TRUE INFO
      dummy.push_back(islc.truth.position.x);
      dummy.push_back(islc.truth.position.y);
      dummy.push_back(islc.truth.position.z);
      CHIvtxMCpi.push_back(dummy);
      dummy.clear();
    }
  }//MC informations storage

  if(ipfp_mu!=-1)
  {
    //RECO MUON TRACK INFO
    CHIbest_plane_mu.push_back(islc.reco.pfp[ipfp_mu].trk.bestplane);

    MUndaughters_reco.push_back(islc.reco.pfp[ipfp_mu].ndaughters);

    dummy.push_back(islc.reco.pfp[ipfp_mu].trk.start.x);
    dummy.push_back(islc.reco.pfp[ipfp_mu].trk.start.y);
    dummy.push_back(islc.reco.pfp[ipfp_mu].trk.start.z);
    CHIstart3DrecoMU.push_back(dummy);
    dummy.clear();

    dummy.push_back(islc.reco.pfp[ipfp_mu].trk.end.x);
    dummy.push_back(islc.reco.pfp[ipfp_mu].trk.end.y);
    dummy.push_back(islc.reco.pfp[ipfp_mu].trk.end.z);
    CHIend3DrecoMU.push_back(dummy);
    dummy.clear();

    CHItlRECOmu.push_back(islc.reco.pfp[ipfp_mu].trk.len);

    vdirx_mu.push_back(islc.reco.pfp[ipfp_mu].trk.dir.x);
    vdiry_mu.push_back(islc.reco.pfp[ipfp_mu].trk.dir.y);
    vdirz_mu.push_back(islc.reco.pfp[ipfp_mu].trk.dir.z);

    // VERTEX TRUE INFO
    dummy.push_back(islc.vertex.x);
    dummy.push_back(islc.vertex.y);
    dummy.push_back(islc.vertex.z);
    CHIvtxRECOmu.push_back(dummy);
    dummy.clear();

  }
          
  //RECO PROTON TRACK INFO
  for(int i=0; i<int(v_ipfp_pro.size()); i++)
  {
    int ind_ipfp=v_ipfp_pro[i];

    CHIbest_plane_pro.push_back(islc.reco.pfp[ind_ipfp].trk.bestplane);

    PROndaughters_reco.push_back(islc.reco.pfp[ind_ipfp].ndaughters);
          
    dummy.push_back(islc.reco.pfp[ind_ipfp].trk.start.x);
    dummy.push_back(islc.reco.pfp[ind_ipfp].trk.start.y);
    dummy.push_back(islc.reco.pfp[ind_ipfp].trk.start.z);
    CHIstart3DrecoPRO.push_back(dummy);
    dummy.clear();

    dummy.push_back(islc.reco.pfp[ind_ipfp].trk.end.x);
    dummy.push_back(islc.reco.pfp[ind_ipfp].trk.end.y);
    dummy.push_back(islc.reco.pfp[ind_ipfp].trk.end.z);
    CHIend3DrecoPRO.push_back(dummy);
    dummy.clear();
          
    CHItlRECOpro.push_back(islc.reco.pfp[ind_ipfp].trk.len);

    vdirx_pro.push_back(islc.reco.pfp[ind_ipfp].trk.dir.x);
    vdiry_pro.push_back(islc.reco.pfp[ind_ipfp].trk.dir.y);
    vdirz_pro.push_back(islc.reco.pfp[ind_ipfp].trk.dir.z);

    // VERTEX TRUE INFO
    dummy.push_back(islc.vertex.x);
    dummy.push_back(islc.vertex.y);
    dummy.push_back(islc.vertex.z);
    CHIvtxRECOpro.push_back(dummy);
    dummy.clear();
  }

  //RECO PION TRACK INFO
  for(int i=0; i<int(v_ipfp_pi.size()); i++)
  {
    int index_pfp_pi=v_ipfp_pi[i];

    CHIbest_plane_pi.push_back(islc.reco.pfp[index_pfp_pi].trk.bestplane);
          
    dummy.push_back(islc.reco.pfp[index_pfp_pi].trk.start.x);
    dummy.push_back(islc.reco.pfp[index_pfp_pi].trk.start.y);
    dummy.push_back(islc.reco.pfp[index_pfp_pi].trk.start.z);
    CHIstart3DrecoPI.push_back(dummy);
    dummy.clear();

    dummy.push_back(islc.reco.pfp[index_pfp_pi].trk.end.x);
    dummy.push_back(islc.reco.pfp[index_pfp_pi].trk.end.y);
    dummy.push_back(islc.reco.pfp[index_pfp_pi].trk.end.z);
    CHIend3DrecoPI.push_back(dummy);
    dummy.clear();
          
    CHItlRECOpi.push_back(islc.reco.pfp[index_pfp_pi].trk.len);

    // VERTEX TRUE INFO
    dummy.push_back(islc.vertex.x);
    dummy.push_back(islc.vertex.y);
    dummy.push_back(islc.vertex.z);
    CHIvtxRECOpi.push_back(dummy);
    dummy.clear(); 
  }
          
  //CALORIMETRIA MUONI
  if(ipfp_mu!=-1)
  {
    //CALORIMETRIA BESTPLANE
    int pplane=-1;
    int bestplane=islc.reco.pfp[ipfp_mu].trk.bestplane;
    if(bestplane!=-1 && islc.reco.pfp[ipfp_mu].trk.calo[bestplane].nhit>0){pplane=bestplane;}
    else{pplane=2;}
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[pplane].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].dedx<500.)
      {
        dEdxbestplane.push_back(islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].dedx);
        rrbestplane.push_back(islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].rr);
        pitchbestplane.push_back(islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].pitch);
        hitxbestplane.push_back(islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].x);
        hitybestplane.push_back(islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].y);
        hitzbestplane.push_back(islc.reco.pfp[ipfp_mu].trk.calo[pplane].points[ihit].z);
      }
    }
    dEdx_bestplane_MU.push_back(dEdxbestplane);
    dEdxbestplane.clear();
    rr_bestplane_MU.push_back(rrbestplane);
    rrbestplane.clear();
    pitch_bestplane_MU.push_back(pitchbestplane);
    pitchbestplane.clear();
    hitx_bestplane_MU.push_back(hitxbestplane);
    hitxbestplane.clear();
    hity_bestplane_MU.push_back(hitybestplane);
    hitybestplane.clear();
    hitz_bestplane_MU.push_back(hitzbestplane);
    hitzbestplane.clear();

    //CALORIMETRIA IND1
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[0].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ipfp_mu].trk.calo[0].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[0].points[ihit].dedx<500.)
      {
        dEdxind1.push_back(islc.reco.pfp[ipfp_mu].trk.calo[0].points[ihit].dedx);
        rrind1.push_back(islc.reco.pfp[ipfp_mu].trk.calo[0].points[ihit].rr);
        pitch_ind1.push_back(islc.reco.pfp[ipfp_mu].trk.calo[0].points[ihit].pitch);
      }
    }
    dEdx_ind1_MU.push_back(dEdxind1);
    dEdxind1.clear();
    rr_ind1_MU.push_back(rrind1);
    rrind1.clear();
    pitch_ind1_MU.push_back(pitch_ind1);
    pitch_ind1.clear();

    //CALORIMETRIA IND2
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[1].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ipfp_mu].trk.calo[1].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[1].points[ihit].dedx<500.)
      {
        dEdxind2.push_back(islc.reco.pfp[ipfp_mu].trk.calo[1].points[ihit].dedx);
        rrind2.push_back(islc.reco.pfp[ipfp_mu].trk.calo[1].points[ihit].rr);
        pitch_ind2.push_back(islc.reco.pfp[ipfp_mu].trk.calo[1].points[ihit].pitch);
      }
    }
    dEdx_ind2_MU.push_back(dEdxind2);
    dEdxind2.clear();
    rr_ind2_MU.push_back(rrind2);
    rrind2.clear();
    pitch_ind2_MU.push_back(pitch_ind2);
    pitch_ind2.clear();

    //CALORIMENTRIA COLLECTION
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx<500.)
      {    
        wire.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].wire);                                  
        CHIdQdx.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dqdx);
        CHIde.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx);           
        CHIpitch.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].pitch);
        CHIr.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr);
        CHImult.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].mult);
        CHIwidth.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].width);
        CHIintegral.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].integral);
        CHIsumadc.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].sumadc);
        CHIx.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x);
        CHIy.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].y);
        CHIz.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z);
      }
    }//loop over all hits of that particle

    CHIwireMU.push_back(wire);
    CHIdQdxMU.push_back(CHIdQdx);
    CHIdeMU.push_back(CHIde);
    CHIrrMUninv.push_back(CHIr);
    CHIwidthMU.push_back(CHIwidth);
    CHIpitchMU.push_back(CHIpitch);
    CHIintegralMU.push_back(CHIintegral);
    CHIsumadcMU.push_back(CHIsumadc);
    CHImultMU.push_back(CHImult);
    CHIxMU.push_back(CHIx);
    CHIyMU.push_back(CHIy);
    CHIzMU.push_back(CHIz);

    std::reverse(CHIde.begin(),CHIde.end());
    std::reverse(CHIr.begin(), CHIr.end());
    std::reverse(CHIpitch.begin(), CHIpitch.end());

    double Eint=0;
    for(int i=0; i<int(CHIde.size()); i++){ Eint = Eint + CHIde[i]*CHIpitch[i]; CHIEint.push_back(Eint); }     

    CHIEintMUallt.push_back(CHIEint);    
    CHIrrMU.push_back(CHIr);
            
    wire.clear();
    CHIdQdx.clear();
    CHIx.clear();
    CHIy.clear();
    CHIz.clear();
    CHImult.clear();
    CHIpitch.clear();
    CHIEint.clear();
    CHIde.clear();
    CHIr.clear();
    CHIwidth.clear();
    CHIintegral.clear();
    CHIsumadc.clear();

  }

  //CALORIMETRIA PROTONI

  for(int i=0; i<int(v_ipfp_pro.size()); i++)
  {
    int ind_ipfp=v_ipfp_pro[i];

    //CALORIMETRIA BESTPLANE
    int pplane=-1;
    int bestplane=islc.reco.pfp[ind_ipfp].trk.bestplane;
    if(bestplane!=-1 && islc.reco.pfp[ind_ipfp].trk.calo[bestplane].nhit>0){pplane=bestplane;}
    else{pplane=2;}
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ind_ipfp].trk.calo[pplane].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].dedx>0.5 && islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].dedx<500.)
      {
        dEdxbestplane.push_back(islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].dedx);
        rrbestplane.push_back(islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].rr);
        pitchbestplane.push_back(islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].pitch);
        hitxbestplane.push_back(islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].x);
        hitybestplane.push_back(islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].y);
        hitzbestplane.push_back(islc.reco.pfp[ind_ipfp].trk.calo[pplane].points[ihit].z);
      }
    }
    dEdx_bestplane_PRO.push_back(dEdxbestplane);
    dEdxbestplane.clear();
    rr_bestplane_PRO.push_back(rrbestplane);
    rrbestplane.clear();
    pitch_bestplane_PRO.push_back(pitchbestplane);
    pitchbestplane.clear();
    hitx_bestplane_PRO.push_back(hitxbestplane);
    hitxbestplane.clear();
    hity_bestplane_PRO.push_back(hitybestplane);
    hitybestplane.clear();
    hitz_bestplane_PRO.push_back(hitzbestplane);
    hitzbestplane.clear();

    //CALORIMETRIA IND1
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ind_ipfp].trk.calo[0].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ind_ipfp].trk.calo[0].points[ihit].dedx>0.5 && islc.reco.pfp[ind_ipfp].trk.calo[0].points[ihit].dedx<500.)
      {
        dEdxind1.push_back(islc.reco.pfp[ind_ipfp].trk.calo[0].points[ihit].dedx);
        rrind1.push_back(islc.reco.pfp[ind_ipfp].trk.calo[0].points[ihit].rr);
        pitch_ind1.push_back(islc.reco.pfp[ind_ipfp].trk.calo[0].points[ihit].pitch);
      }
    }
    dEdx_ind1_PRO.push_back(dEdxind1);
    dEdxind1.clear();
    rr_ind1_PRO.push_back(rrind1);
    rrind1.clear();
    pitch_ind1_PRO.push_back(pitch_ind1);
    pitch_ind1.clear();

    //CALORIMETRIA IND2
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ind_ipfp].trk.calo[1].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ind_ipfp].trk.calo[1].points[ihit].dedx>0.5 && islc.reco.pfp[ind_ipfp].trk.calo[1].points[ihit].dedx<500.)
      {
        dEdxind2.push_back(islc.reco.pfp[ind_ipfp].trk.calo[1].points[ihit].dedx);
        rrind2.push_back(islc.reco.pfp[ind_ipfp].trk.calo[1].points[ihit].rr);
        pitch_ind2.push_back(islc.reco.pfp[ind_ipfp].trk.calo[1].points[ihit].pitch);
      }
    }
    dEdx_ind2_PRO.push_back(dEdxind2);
    dEdxind2.clear();
    rr_ind2_PRO.push_back(rrind2);
    rrind2.clear();
    pitch_ind2_PRO.push_back(pitch_ind2);
    pitch_ind2.clear();

    //CALORIMETRIA COLLECTION
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ind_ipfp].trk.calo[2].points.size(); ++ihit )
    {
      if(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].dedx>0.5 && islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].dedx<500.)
      {   
        CHIdQdx.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].dqdx);
        CHIde.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].dedx);           
        vector_active.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].dedx);
        CHIpitch.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].pitch);
        CHIr.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].rr);
        CHImult.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].mult);
        CHIwidth.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].width);
        CHIintegral.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].integral);
        CHIsumadc.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].sumadc);
        CHIx.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].x);
        CHIy.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].y);
        CHIz.push_back(islc.reco.pfp[ind_ipfp].trk.calo[2].points[ihit].z);
      }
    }//loop over all hits of that particle

    CHIdQdxPRO.push_back(CHIdQdx);
    CHIdePRO.push_back(CHIde);
    CHIrrPROninv.push_back(CHIr);
    CHIwidthPRO.push_back(CHIwidth);
    CHIpitchPRO.push_back(CHIpitch);
    CHIintegralPRO.push_back(CHIintegral);
    CHImultPRO.push_back(CHImult);
    CHIsumadcPRO.push_back(CHIsumadc);
    CHIxPRO.push_back(CHIx);
    CHIyPRO.push_back(CHIy);
    CHIzPRO.push_back(CHIz);

    std::reverse(CHIde.begin(),CHIde.end());
    std::reverse(CHIr.begin(), CHIr.end());
    std::reverse(CHIpitch.begin(), CHIpitch.end());

    double Eint=0;
    for(int i=0; i<int(CHIde.size()); i++){ Eint = Eint + CHIde[i]*CHIpitch[i]; CHIEint.push_back(Eint);}     

    CHIEintPROallt.push_back(CHIEint);    
    CHIrrPRO.push_back(CHIr);
            
    CHIdQdx.clear();
    CHIx.clear();
    CHIy.clear();
    CHIz.clear();
    CHImult.clear();
    CHIpitch.clear();
    CHIEint.clear();
    CHIde.clear();
    CHIr.clear();
    CHIwidth.clear();
    CHIintegral.clear();
    CHIsumadc.clear();
  }//loop over all proton track in that slice

  //CALORIMETRIA PIONI
  for(int i=0; i<int(v_ipfp_pi.size()); i++)
  {
    int index_pfp_pi=v_ipfp_pi[i];

    //CALORIMETRIA BESTPLANE
    int pplane=-1;
    int bestplane=islc.reco.pfp[index_pfp_pi].trk.bestplane;
    if(bestplane!=-1 && islc.reco.pfp[index_pfp_pi].trk.calo[bestplane].nhit>0){pplane=bestplane;}
    else{pplane=2;}
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points.size(); ++ihit )
    {
      if(islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].dedx>0.5 && islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].dedx<500.)
      {
        dEdxbestplane.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].dedx);
        rrbestplane.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].rr);
        pitchbestplane.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].pitch);
        hitxbestplane.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].x);
        hitybestplane.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].y);
        hitzbestplane.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[pplane].points[ihit].z);
      }
    }
    dEdx_bestplane_PI.push_back(dEdxbestplane);
    dEdxbestplane.clear();
    rr_bestplane_PI.push_back(rrbestplane);
    rrbestplane.clear();
    pitch_bestplane_PI.push_back(pitchbestplane);
    pitchbestplane.clear();
    hitx_bestplane_PI.push_back(hitxbestplane);
    hitxbestplane.clear();
    hity_bestplane_PI.push_back(hitybestplane);
    hitybestplane.clear();
    hitz_bestplane_PI.push_back(hitzbestplane);
    hitzbestplane.clear();

    //CALORIMETRIA IND1
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[index_pfp_pi].trk.calo[0].points.size(); ++ihit )
    {
      if(islc.reco.pfp[index_pfp_pi].trk.calo[0].points[ihit].dedx>0.5 && islc.reco.pfp[index_pfp_pi].trk.calo[0].points[ihit].dedx<500.)
      {
        dEdxind1.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[0].points[ihit].dedx);
        rrind1.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[0].points[ihit].rr);
        pitch_ind1.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[0].points[ihit].pitch);
      }
    }
    dEdx_ind1_PI.push_back(dEdxind1);
    dEdxind1.clear();
    rr_ind1_PI.push_back(rrind1);
    rrind1.clear();
    pitch_ind1_PI.push_back(pitch_ind1);
    pitch_ind1.clear();
    

    //CALORIMETRIA IND2
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[index_pfp_pi].trk.calo[1].points.size(); ++ihit )
    {
      if(islc.reco.pfp[index_pfp_pi].trk.calo[1].points[ihit].dedx>0.5 && islc.reco.pfp[index_pfp_pi].trk.calo[1].points[ihit].dedx<500.)
      {
        dEdxind2.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[1].points[ihit].dedx);
        rrind2.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[1].points[ihit].rr);
        pitch_ind2.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[1].points[ihit].pitch);
      }
    }
    dEdx_ind2_PI.push_back(dEdxind2);
    dEdxind2.clear();
    rr_ind2_PI.push_back(rrind2);
    rrind2.clear();
    pitch_ind2_PI.push_back(pitch_ind2);
    pitch_ind2.clear();


    //CALORIMETRIA COLLECTION
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[index_pfp_pi].trk.calo[2].points.size(); ++ihit )
    {
      if(islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].dedx>0.5 && islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].dedx<500.)
      {   
        CHIde.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].dedx);           
        CHIr.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].rr);
        CHIx.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].x);
        CHIy.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].y);
        CHIz.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].z);
        CHIpitch.push_back(islc.reco.pfp[index_pfp_pi].trk.calo[2].points[ihit].pitch);
      }
    }//loop over all hits of that particle


    CHIdePI.push_back(CHIde);
    CHIrrPIninv.push_back(CHIr);
    CHIxPI.push_back(CHIx);
    CHIyPI.push_back(CHIy);
    CHIzPI.push_back(CHIz);
    CHIpitchPI.push_back(CHIpitch);
 
    CHIx.clear();
    CHIy.clear();
    CHIz.clear();
    CHIde.clear();
    CHIr.clear();
    CHIpitch.clear();

  }//loop over all pions track in that slice

}//loop over all slices

  return vector_active;
});
