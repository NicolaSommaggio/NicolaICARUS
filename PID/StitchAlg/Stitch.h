////////////////////////////////////////////////////////////////////
// Definition of the stitching function
//
// Author: Alessandro Maria Ricci
////////////////////////////////////////////////////////////////////

#ifndef STITCH_H
#define STITCH_H

#include "Selection.h"
#include "TrackMomentumCalculator.h"
#include <fstream>


bool isInFV (const caf::Proxy<caf::SRSlice>& islc)
{
  double x = islc.vertex.x;
  double y = islc.vertex.y;
  double z = islc.vertex.z;
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
  //need to add a check to avoid having the vtx slice in one cryo and the end track in the other cryo

    //fiducial volume for the dangling cable 
  if(x>210.0 && y > 60.0 && z> 290.0 && z< 390.0 ) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
			( x >  61.94 + 25 && x <  358.49 - 25 )) &&
		  ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
		  ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
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

bool all_contained ( const caf::Proxy<caf::SRSlice>& islc ) 
{ 

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

//FNAL
TFile* file = TFile::Open("/exp/icarus/data/users/nsommagg/RefCurvesChi2.root");
auto dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
auto dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
auto dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
auto dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

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
        islc.reco.pfp[ipfp].trackScore>0.4 &&
        //(output[0]<30 && output[1]>60) &&
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

namespace ana {

    namespace param {
        bool print = false;         // print on terminal
        bool max2Seg = true;        // consider only nuMuCC with maximum 2 segments (ONLY MC)
        bool checkVtx = true;       // check if the vertex is between the two tracks
        float mcLmin = 50.;         // minimal length of the track on mc (ONLY MC)
        float trkLmin = 20.;        // minimal length of the reconstructed track
        float trkScore = 0.4;       // the MVA score that determines how track/shower like a PFP is
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

    struct muon {
        std::vector<int> slcID;
        int G4ID = -2;
        size_t Nseg = 0;
        int run; //-->NS
        int evt; //-->NS
        bool calib_Tss_found = false;
        bool calib_theta_found = false;
        bool calib_Dss_found = false;
        bool split_intra_slice = false;
        bool split_extra_slice = false;
        bool correct_stitched_intra = false;
        bool wrong_stitched_intra = false;
        bool correct_stitched_extra = false;
        bool wrong_stitched_extra = false;
        bool vertex_failure = false;
        bool Tss_failure = false;
        bool theta_failure = false;
        bool Dss_failure = false;
    };

    struct calibration {
        std::vector<double> Tss_MuMu_IS;        // IS = intra-slice, ES = extra-slice
        std::vector<double> Tss_MuOther_IS;
        std::vector<double> Tss_Data_IS;        // Data = distribution in the real data
        std::vector<double> Tss_MuMu_ES;
        std::vector<double> Tss_MuOther_ES;
        std::vector<double> Tss_Data_ES;
        std::vector<double> Theta_MuMu_IS;
        std::vector<double> Theta_MuOther_IS;
        std::vector<double> Theta_Data_IS;
        std::vector<double> Theta_MuMu_ES;
        std::vector<double> Theta_MuOther_ES;
        std::vector<double> Theta_Data_ES;
        std::vector<double> Dss_MuMu_IS;
        std::vector<double> Dss_MuOther_IS;
        std::vector<double> Dss_Data_IS;
        std::vector<double> Dss_MuMu_ES;
        std::vector<double> Dss_MuOther_ES;
        std::vector<double> Dss_Data_ES;
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

    std::vector<std::array<int,3>> NS_evts;

    // stitching function in the event
    std::vector<PFP> stitch(const caf::SRSpillProxy* sr, std::vector<muon> &MCmuons, calibration &calib, bool calibration, bool mc) {
        std::vector<PFP> stitch;
        std::vector<PFP> pfp;
        size_t nStitch = 0;
        int muonIdx = -1;
        if(mc) {
            // select mu primary daugthers of numuCC interactions, contained in a single cryostat
            // and with a minimal length of the track
            for(size_t i=0; i<sr->mc.nnu; i++) {
                for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                    if(sr->mc.nu.at(i).prim.at(j).pdg != std::abs(13)) continue;
                    //if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                    //if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                    muon mu;
                    mu.G4ID = sr->mc.nu.at(i).prim.at(j).G4ID;
                    mu.run = sr->hdr.run; //-->NS
                    mu.evt = sr->hdr.evt; //-->NS
                    MCmuons.push_back(mu);
                }
            }
            if(MCmuons.empty()) return stitch;
        }

        // loop on the slices to import the pfp
        for(int i=0; i<sr->nslc; i++) {

            //bool stitched_NS = false;
            //for(const auto &NS_evt : NS_evts)
            //{
            //  if( ((int)sr->hdr.run == NS_evt[0]) && ((int)sr->hdr.evt == NS_evt[1]) )
            //  {
            //    cout << "match" << endl;
            //    stitched_NS = true;
            //  }
            //}

            //if(stitched_NS == false)continue;

            //cout << (int)sr->hdr.run << " " << (int)sr->hdr.evt << " ";

            //check the truth matching
            if(!is_good_truth_matching(sr->slc.at(i)))continue;

            // check if the slice is in the fiducial volume
            if(!isInFV(sr->slc.at(i)))continue;
            //-> if(!kIcarus202412RecoFiducial(sr->slc.at(i))) continue;
            // apply the 1muNp selection if data
            //-> if(!mc && !kIcarus202507Contained1muNp(sr->slc.at(i))) continue;
            //->nicola if(!all_contained(sr->slc.at(i)))continue;

            //cout << "passed ";

            //for data
            if(!mc) muonIdx = kIcarus202507MuonIdx(sr->slc.at(i));

            for(size_t j=0; j<sr->slc.at(i).reco.npfp; j++) {
                // check if the pfp is a muon
                //if(!kIcarus202401MuonTrack(sr->slc.at(i).reco.pfp.at(j), param::trkScore, param::trkLmin)) continue;

                if(!is_muon(sr->slc.at(i),int(j)))continue;

                //cout << " pfp found " << endl;
                //double pid_class = PIDclass(sr->slc.at(i),j);
                //if(!(pid_class == 0 || pid_class == 1))continue;

                std::vector<double> output;
                std::vector<double> dedx;
                std::vector<double> rr;
                for ( std::size_t ihit(0); ihit < sr->slc.at(i).reco.pfp[j].trk.calo[2].points.size(); ++ihit )
                {   
                    dedx.push_back(sr->slc.at(i).reco.pfp[j].trk.calo[2].points[ihit].dedx);
                    rr.push_back(sr->slc.at(i).reco.pfp[j].trk.calo[2].points[ihit].rr);
                } 
                output = chi2_ALG(dedx,rr,0.0,25.0);

                if(!(output[0]<30 && output[1]>60))continue;

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
                P.id = sr->slc.at(i).reco.pfp.at(j).id;
                // data
                if(!mc) {
                    int idx = j;
                    if(idx==muonIdx) P.muon_1muNp = true;
                }
                // MC
                if(mc) {
                    if(!std::isnan(sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.G4ID)) P.G4ID = sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.G4ID;
                    else P.G4ID = -1;
                    if(!std::isnan(sr->slc.at(i).reco.pfp.at(j).trk.truth.p.pdg)) P.pdg = sr->slc.at(i).reco.pfp.at(j).trk.truth.p.pdg;
                    else P.pdg = -1;
                    if(!std::isnan(sr->slc.at(i).reco.pfp.at(j).trk.truth.p.length)) P.MClen = sr->slc.at(i).reco.pfp.at(j).trk.truth.p.length;
                    else P.MClen = -1;
                    for(size_t k=0; k<MCmuons.size(); k++) {
                        if(MCmuons.at(k).G4ID == P.G4ID) {
                            MCmuons.at(k).Nseg++;
                            MCmuons.at(k).slcID.push_back(P.slcID);
                        }
                    }
                }
                pfp.push_back(P);
            }
        }
        
        //for(size_t i=0; i<MCmuons.size(); i++) 
        //{
            //if(MCmuons.at(i).Nseg > 2){cout << MCmuons.at(i).run << " " << MCmuons.at(i).evt << endl;}
        //}

        // Remove muons broken with more than 2 segments
        if(mc && param::max2Seg) {
            // delete all elements that satisfy the condition NSeg>2
            MCmuons.erase(std::remove_if(MCmuons.begin(), MCmuons.end(), [](const muon& mu) {
                return mu.Nseg > 2;
            }), MCmuons.end());
        }
        for(size_t i=0; i<MCmuons.size(); i++) {
            if(MCmuons.at(i).Nseg<2) continue;
            for(size_t j=0; j<MCmuons.at(i).slcID.size(); j++) {
                if(MCmuons.at(i).slcID.at(j)!=MCmuons.at(i).slcID.at(0)) MCmuons.at(i).split_extra_slice = true;
                for(size_t k=0; k<MCmuons.at(i).slcID.size(); k++) {
                    if(j==k) continue;
                    if(MCmuons.at(i).slcID.at(j)==MCmuons.at(i).slcID.at(k)) {
                        MCmuons.at(i).split_intra_slice = true;
                        break;
                    }
                }
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


            //bool stitched_NS = false;
            //for(const auto &NS_evt : NS_evts)
            //{
              //if( ((int)pfp.at(i).run == NS_evt[0]) && ((int)pfp.at(i).evt== NS_evt[1]) && ((int)pfp.at(i).slice_index == NS_evt[2]) &&
                  //((int)pfp.at(j).run == NS_evt[0]) && ((int)pfp.at(j).evt== NS_evt[1]) && ((int)pfp.at(j).slice_index == NS_evt[2]) 
                //)
              //{
                //stitched_NS = true;
              //}
            //}

            //if(stitched_NS==false)continue;

                // skip the evaluation if both pfps do not pass the 1muNp selection in the data
                if(!mc && !pfp.at(i).muon_1muNp && !pfp.at(j).muon_1muNp) continue;
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
                if(param::checkVtx && ordz[1] == vtx.z) {
                    // fill register of the MC muons
                    if(mc) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(pfp.at(i).G4ID==MCmuons.at(k).G4ID && pfp.at(j).G4ID==MCmuons.at(k).G4ID)
                                MCmuons.at(k).vertex_failure = true;
                        }
                    }
                    continue;
                }
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
                if(Zoverlap) {
                    if(trk.at(0).slcID==trk.at(1).slcID) {
                        // Tss intra-slice calibration
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).calib_Tss_found == false) {
                                    calib.Tss_MuMu_IS.push_back( std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] );
                                    MCmuons.at(k).calib_Tss_found = true;
                                }
                                else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                         (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) )
                                    calib.Tss_MuOther_IS.push_back( std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] );
                            }
                        }
                        else calib.Tss_Data_IS.push_back( std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] );
                        // selection
                        if(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= param::slc::Tss) {
                            // fill register of the MC muons
                            if(mc)
                                for(size_t k=0; k<MCmuons.size(); k++)
                                    if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID)
                                        MCmuons.at(k).Tss_failure = true;
                            if(!calibration) continue;
                        }
                    }
                    else {
                        // Tss extra-slice calibration
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).calib_Tss_found == false) {
                                    calib.Tss_MuMu_ES.push_back( std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] );
                                    MCmuons.at(k).calib_Tss_found = true;
                                }
                                else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                         (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) )
                                    calib.Tss_MuOther_ES.push_back( std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] );
                            }
                        }
                        else calib.Tss_Data_ES.push_back( std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] );
                        // selection
                        if(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= param::evt::Tss) {
                            // fill register of the MC muons
                            if(mc)
                                for(size_t k=0; k<MCmuons.size(); k++)
                                    if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID)
                                        MCmuons.at(k).Tss_failure = true;
                            if(!calibration) continue;
                        }
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
                if(trk.at(0).slcID==trk.at(1).slcID) {
                    // theta intra-slice calibration
                    if(mc) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).calib_theta_found == false) {
                                calib.Theta_MuMu_IS.push_back( phi*180/3.141592 );
                                MCmuons.at(k).calib_theta_found = true;
                            }
                            else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                     (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) )
                                calib.Theta_MuOther_IS.push_back( phi*180/3.141592 );
                        }
                    }
                    else calib.Theta_Data_IS.push_back( phi*180/3.141592 );
                    // selection
                    if(phi*180/3.141592 >= param::slc::theta) {
                        // fill register of the MC muons
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID)
                                    MCmuons.at(k).theta_failure = true;
                            }
                        }
                        if(!calibration) continue;
                    }
                }
                else {
                    // theta extra-slice calibration
                    if(mc) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).calib_theta_found == false) {
                                calib.Theta_MuMu_ES.push_back( phi*180/3.141592 );
                                MCmuons.at(k).calib_theta_found = true;
                            }
                            else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                     (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) )
                                calib.Theta_MuOther_ES.push_back( phi*180/3.141592 );
                        }
                    }
                    else calib.Theta_Data_ES.push_back( phi*180/3.141592 );
                    // selection
                    if(phi*180/3.141592 >= param::evt::theta) {
                        // fill register of the MC muons
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID)
                                    MCmuons.at(k).theta_failure = true;
                            }
                        }
                        if(!calibration) continue;
                    }
                }
                // check distance between tracks
                std::vector<float> length;
                length.push_back(Ra);
                length.push_back(Rb);
                std::sort(length.begin(), length.end());
                float dist = std::hypot(trk.at(1).start.x-trk.at(0).end.x, trk.at(1).start.y-trk.at(0).end.y, trk.at(1).start.z-trk.at(0).end.z);
                if(trk.at(0).slcID==trk.at(1).slcID) {
                    // Dss intra-slice calibration
                    if(mc) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).calib_Dss_found == false) {
                                calib.Dss_MuMu_IS.push_back( dist/length[0] );
                                MCmuons.at(k).calib_Dss_found = true;
                            }
                            else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                     (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) )
                                calib.Dss_MuOther_IS.push_back( dist/length[0] );
                        }
                    }
                    else calib.Dss_Data_IS.push_back( dist/length[0] );
                    // selection
                    if(dist/length[0] >= param::slc::Dss) {
                        // fill register of the MC muons
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID)
                                    MCmuons.at(k).Dss_failure = true;
                            }
                        }
                        if(!calibration) continue;
                    }
                }
                else {
                    // Dss extra-slice calibration
                    if(mc) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).calib_Dss_found == false) {
                                calib.Dss_MuMu_ES.push_back( dist/length[0] );
                                MCmuons.at(k).calib_Dss_found = true;
                            }
                            else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                     (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) )
                                calib.Dss_MuOther_ES.push_back( dist/length[0] );
                        }
                    }
                    else calib.Dss_Data_ES.push_back( dist/length[0] );
                    // selection
                    if(dist/length[0] >= param::evt::Dss) {
                        // fill register of the MC muons
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID)
                                    MCmuons.at(k).Dss_failure = true;
                            }
                        }
                        if(!calibration) continue;
                    }
                }
                // fill register of the MC muons
                if(mc) {
                    for(size_t k=0; k<MCmuons.size(); k++) {
                        if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) {
                            if(trk.at(0).slcID==trk.at(1).slcID) MCmuons.at(k).correct_stitched_intra = true;
                            else MCmuons.at(k).correct_stitched_extra = true;
                        }
                        else if( (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) ||
                                 (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ) {
                            if(trk.at(0).slcID==trk.at(1).slcID) MCmuons.at(k).wrong_stitched_intra = true;
                            else MCmuons.at(k).wrong_stitched_extra = true;
                        }
                    }
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
                if(!mc) sth.muon_1muNp = true;
                else {
                    for(size_t k=0; k<MCmuons.size(); k++) {
                        if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) sth.nuMuCC = true;
                    }
                    sth.MClen = trk.at(0).MClen;                        // MC length
                    if(sth.MClen > 0) sth.MCp_muon = GetTrackMomentum(sth.MClen, 13);     // MC momentum
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


    // fill the histogram with the output of stitch
    void FillStitchHist(std::vector<muon> &MCmuons, std::vector<double> &hist) {
        for(size_t i=0; i<MCmuons.size(); i++) {
            //cout << MCmuons.at(i).run << " " << MCmuons.at(i).evt << endl;
            // only intra-slice split muons
            if(MCmuons.at(i).split_intra_slice && !MCmuons.at(i).split_extra_slice) {
                if(MCmuons.at(i).correct_stitched_intra && !MCmuons.at(i).wrong_stitched_intra) hist.push_back(0.);        // correct stitch
                else if(!MCmuons.at(i).correct_stitched_intra && MCmuons.at(i).wrong_stitched_intra) hist.push_back(1.);   // wrong stitch
                else if(MCmuons.at(i).correct_stitched_intra && MCmuons.at(i).wrong_stitched_intra) hist.push_back(2.);    // correct + wrong stitch
                else {
                    hist.push_back(3.);                                         // no stitch
                    if(MCmuons.at(i).vertex_failure) hist.push_back(5.);        // no stitch due to vertex
                    else if(MCmuons.at(i).Tss_failure) hist.push_back(6.);      // no stitch due to Tss
                    else if(MCmuons.at(i).theta_failure) hist.push_back(7.);    // no stitch due to theta
                    else if(MCmuons.at(i).Dss_failure) hist.push_back(8.);      // no stitch due to Dss
                }
            }
            // only extra-slice split muons
            else if(!MCmuons.at(i).split_intra_slice && MCmuons.at(i).split_extra_slice) {
                if(MCmuons.at(i).correct_stitched_extra && !MCmuons.at(i).wrong_stitched_extra) hist.push_back(10.);        // correct stitch
                else if(!MCmuons.at(i).correct_stitched_extra && MCmuons.at(i).wrong_stitched_extra) hist.push_back(11.);   // wrong stitch
                else if(MCmuons.at(i).correct_stitched_extra && MCmuons.at(i).wrong_stitched_extra) hist.push_back(12.);    // correct + wrong stitch
                else {
                    hist.push_back(13.);                                        // no stitch
                    if(MCmuons.at(i).vertex_failure) hist.push_back(15.);       // no stitch due to vertex
                    else if(MCmuons.at(i).Tss_failure) hist.push_back(16.);     // no stitch due to Tss
                    else if(MCmuons.at(i).theta_failure) hist.push_back(17.);   // no stitch due to theta
                    else if(MCmuons.at(i).Dss_failure) hist.push_back(18.);     // no stitch due to Dss
                }
            }
            // both intra and extra-slice split muons
            else if(MCmuons.at(i).split_intra_slice && MCmuons.at(i).split_extra_slice) {
                if( MCmuons.at(i).correct_stitched_intra && MCmuons.at(i).correct_stitched_extra &&
                   !MCmuons.at(i).wrong_stitched_intra && !MCmuons.at(i).wrong_stitched_extra) hist.push_back(20.);        // correct stitch
                else if(!MCmuons.at(i).correct_stitched_intra && !MCmuons.at(i).correct_stitched_extra &&
                         MCmuons.at(i).wrong_stitched_intra && MCmuons.at(i).wrong_stitched_extra) hist.push_back(21.);   // wrong stitch
                else if(MCmuons.at(i).correct_stitched_intra && MCmuons.at(i).correct_stitched_extra &&
                        MCmuons.at(i).wrong_stitched_intra && MCmuons.at(i).wrong_stitched_extra) hist.push_back(22.);     // correct + wrong stitch
                else {
                    hist.push_back(23.);                                        // no stitch
                    if(MCmuons.at(i).vertex_failure) hist.push_back(25.);       // no stitch due to vertex
                    else if(MCmuons.at(i).Tss_failure) hist.push_back(26.);     // no stitch due to Tss
                    else if(MCmuons.at(i).theta_failure) hist.push_back(27.);   // no stitch due to theta
                    else if(MCmuons.at(i).Dss_failure) hist.push_back(28.);     // no stitch due to Dss
                }
            }
            // non-split muons
            else {
                if(MCmuons.at(i).wrong_stitched_intra && !MCmuons.at(i).wrong_stitched_extra) hist.push_back(30.);          // wrong stitch intra-slice
                else if(!MCmuons.at(i).wrong_stitched_intra && MCmuons.at(i).wrong_stitched_extra) hist.push_back(31.);     // wrong stitch extra-slice
                else if(MCmuons.at(i).wrong_stitched_intra && MCmuons.at(i).wrong_stitched_extra) hist.push_back(32.);      // wrong stitch intra + extra-slice
            }
        }
    }


    // print on terminal the output of stitch
    void PrintStitch(std::vector<PFP> &stitch) {
        for(size_t i=0; i<stitch.size(); i++) {
            std::cout << " " << std::endl;
            std::cout << " " << std::endl;
            std::cout << "==============" << std::endl;
            std::cout << "Run: " << stitch.at(i).run << std::endl;
            std::cout << "Subrun: " << stitch.at(i).subrun << std::endl;
            std::cout << "Event: " << stitch.at(i).evt << std::endl;
            std::cout << "==============" << std::endl;
            std::cout << " " << std::endl;
            std::cout << "===================" << std::endl;
            std::cout << "Stitching No: " << stitch.at(i).nStitch << std::endl;
            std::cout << " " << std::endl;
            std::cout << "slcID 1 = " << stitch.at(i).slcID << std::endl;
            std::cout << "slcID 2 = " << stitch.at(i).slcID2 << std::endl;
            std::cout << "pfpID 1 = " << stitch.at(i).id << std::endl;
            std::cout << "pfpID 2 = " << stitch.at(i).id2 << std::endl;
            std::cout << std::endl;
            //aggiunto da Nicola **** -> looking if the track is stitched correctly
            std::cout << "G4ID 1 = " << stitch.at(i).G4ID << std::endl;
            std::cout << "pdg 1 = " << stitch.at(i).pdg << std::endl;
            std::cout << "G4ID 2 = " << stitch.at(i).G4ID2 << std::endl;
            std::cout << "pdg 2 = " << stitch.at(i).pdg2 << std::endl;
            std::cout << std::endl;
            //****************
            std::cout << "geom = {v,s1,e1,s2,e2}" << std::endl;
            std::cout << "x[" << 0 << "] = " << stitch.at(i).vertex.x <<
                "\t y[" << 0 << "] = " << stitch.at(i).vertex.y <<
                "\t z[" << 0 << "] = " << stitch.at(i).vertex.z << std::endl;
            std::cout << "x[" << 1 << "] = " << stitch.at(i).start.x <<
                "\t y[" << 1 << "] = " << stitch.at(i).start.y <<
                "\t z[" << 1 << "] = " << stitch.at(i).start.z << std::endl;
            std::cout << "x[" << 2 << "] = " << stitch.at(i).end.x <<
                "\t y[" << 2 << "] = " << stitch.at(i).end.y <<
                "\t z[" << 2 << "] = " << stitch.at(i).end.z << std::endl;
            std::cout << "x[" << 3 << "] = " << stitch.at(i).start2.x <<
                "\t y[" << 3 << "] = " << stitch.at(i).start2.y <<
                "\t z[" << 3 << "] = " << stitch.at(i).start2.z << std::endl;
            std::cout << "x[" << 4 << "] = " << stitch.at(i).end2.x <<
                "\t y[" << 4 << "] = " << stitch.at(i).end2.y <<
                "\t z[" << 4 << "] = " << stitch.at(i).end2.z << std::endl;
            std::cout << " " << std::endl;
            std::cout << "L1 = " << stitch.at(i).len << std::endl;
            std::cout << "L1+L2 = " << stitch.at(i).Len << std::endl;
            std::cout << "p1 = " << stitch.at(i).p_muon << std::endl;
            std::cout << "p1+p2 = " << stitch.at(i).P_muon << std::endl;
            std::cout << "E1 = " << std::sqrt( stitch.at(i).p_muon*stitch.at(i).p_muon + mmu*mmu) << std::endl;
            std::cout << "E1+E2 = " << std::sqrt(stitch.at(i).P_muon*stitch.at(i).P_muon + mmu*mmu) << std::endl;
            std::cout << " " << std::endl;
        }
    }

ofstream dump_stitching_AMR("dump_stitching_AMR_chi2_ts4.txt");

    // stitching function
    const SpillMultiVar kStitch([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        std::vector<double> hist;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);

        //aggiunta da nicola ////////////////
        for(size_t i=0; i<sth.size(); i++) 
        {
            dump_stitching_AMR << sth.at(i).run << " " << sth.at(i).evt << " " << sth.at(i).slice_index << " " << sth.at(i).G4ID << " " << sth.at(i).pdg << " " << sth.at(i).G4ID2 << " " << sth.at(i).pdg2 << " " << sth.at(i).len << " " << sth.at(i).Len << " " << sth.at(i).MClen << endl;
        }
        /////////////////////////////////////

        if(sr->hdr.ismc) FillStitchHist(MCmuons, hist);
        if(!sth.empty() && param::print) PrintStitch(sth);
        return hist;
    });

    

    // MC only: return the values of the parameter Tss for intra-slice muon-muon stitching in MC sample
    const SpillMultiVar kTssMuMuIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Tss_MuMu_IS;
    });

    // MC only: return the values of the parameter Tss for intra-slice muon-other stitching in MC sample
    const SpillMultiVar kTssMuOtherIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Tss_MuOther_IS;
    });

    // Data only: return the values of the parameter Tss for intra-slice stitching in the real data
    const SpillMultiVar kTssDataIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Tss_Data_IS;
    });

    // MC only: return the values of the parameter Tss for extra-slice muon-muon stitching in MC sample
    const SpillMultiVar kTssMuMuES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Tss_MuMu_ES;
    });

    // MC only:return the values of the parameter Tss for extra-slice muon-other stitching in MC sample
    const SpillMultiVar kTssMuOtherES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Tss_MuOther_ES;
    });

    // Data only: return the values of the parameter Tss for extra-slice stitching in the real data
    const SpillMultiVar kTssDataES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Tss_Data_ES;
    });

    // MC only: return the values of the parameter Theta for intra-slice muon-muon stitching in MC sample
    const SpillMultiVar kThetaMuMuIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Theta_MuMu_IS;
    });

    // MC only: return the values of the parameter Theta for intra-slice muon-other stitching in MC sample
    const SpillMultiVar kThetaMuOtherIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Theta_MuOther_IS;
    });

    // Data only: return the values of the parameter Theta for intra-slice stitching in the real data
    const SpillMultiVar kThetaDataIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Theta_Data_IS;
    });

    // MC only: return the values of the parameter Theta for extra-slice muon-muon stitching in MC sample
    const SpillMultiVar kThetaMuMuES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Theta_MuMu_ES;
    });

    // MC only: return the values of the parameter Theta for extra-slice muon-other stitching in MC sample
    const SpillMultiVar kThetaMuOtherES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Theta_MuOther_ES;
    });

    // Data only: return the values of the parameter Theta for extra-slice stitching in the real data
    const SpillMultiVar kThetaDataES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Theta_Data_ES;
    });

    // MC only: return the values of the parameter Dss for intra-slice muon-muon stitching in MC sample
    const SpillMultiVar kDssMuMuIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true,sr->hdr.ismc);
        return calib.Dss_MuMu_IS;
    });

    // MC only: return the values of the parameter Dss for intra-slice muon-other stitching in MC sample
    const SpillMultiVar kDssMuOtherIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Dss_MuOther_IS;
    });

    // Data only: return the values of the parameter Dss for intra-slice stitching in the real data
    const SpillMultiVar kDssDataIS([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Dss_Data_IS;
    });

    // MC only: return the values of the parameter Dss for extra-slice muon-muon stitching in MC sample
    const SpillMultiVar kDssMuMuES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Dss_MuMu_ES;
    });

    // MC only: return the values of the parameter Dss for extra-slice muon-other stitching in MC sample
    const SpillMultiVar kDssMuOtherES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Dss_MuOther_ES;
    });

    // Data only: return the values of the parameter Dss for extra-slice stitching in the real data
    const SpillMultiVar kDssDataES([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, true, sr->hdr.ismc);
        return calib.Dss_Data_ES;
    });


    // adapted from kIcarus202401MuonIdx in NumuVarsIcarus202401.cxx
    bool kIcarus202503MuonTrack(const caf::SRPFPProxy &pfp, float trkScore, float trkLmin, size_t &counter) {
        // The (dis)qualification of a slice is based upon the track level information.
        bool muTrack = false;
        if(pfp.trackScore < trkScore) return muTrack;
        if(counter<2) counter=2;
        if(std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x)) return muTrack;
        auto const& trk = pfp.trk;

        //int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
        int plane = 2; // Hard code collection plane for now since induction 2 has peak at higher chi2

        //float Chi2Proton = trk.chi2pid[plane].chi2_proton;
        //float Chi2Muon = trk.chi2pid[plane].chi2_muon;
        auto chi2 = chi2pid::chi2_calculator.calculate_chi2(trk.calo[plane]);
        float Chi2Proton = chi2.chi2_proton;
        float Chi2Muon = chi2.chi2_muon;

        const bool Contained = ( !isnan(trk.end.x) &&
            (( trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5 ) ||
             ( trk.end.x > 61.94 + 5 && trk.end.x < +358.49 - 5 )) &&
            !isnan(trk.end.y) &&
            ( trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5 ) &&
            !isnan(trk.end.z) &&
            ( trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5 ) );
        if(Contained && counter<3) counter=3;
        if(trk.len >= trkLmin && counter<4) counter=4;
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len >= trkLmin );
        if(MaybeMuonContained) muTrack = true;
        return muTrack;
    }


    // MC only: return the number of mu, primary daugher of nuMuCC interactions, contained in a single TPC volume
    const SpillMultiVar kMu([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> mu;
        if(!sr->hdr.ismc) return mu;
        // select mu primary daugthers of numuCC interactions, contained in a single cryostat
        // and with a minimal length of the track
        for(size_t i=0; i<sr->mc.nnu; i++) {
            for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                if(sr->mc.nu.at(i).prim.at(j).pdg != std::abs(13)) continue;
                if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                int MCG4ID = sr->mc.nu.at(i).prim.at(j).G4ID;
                size_t counter = 0;
                // loop on the slices
                for(int k=0; k<sr->nslc; k++) {
                    // check if the slice is in the fiducial volume
                    if(!kIcarus202412RecoFiducial(sr->slc.at(k))) continue;
                    // loop on the pfp particles
                    for(size_t m=0; m<sr->slc.at(k).reco.npfp; m++) {
                        int pfpG4ID = sr->slc.at(k).reco.pfp.at(m).trk.truth.bestmatch.G4ID;
                        if(MCG4ID != pfpG4ID) continue;
                        if(counter<1) counter=1; // muons that passed the RecoFiducial cut
                        // check if the pfp is a muon
                        if(!kIcarus202503MuonTrack(sr->slc.at(k).reco.pfp.at(m), param::trkScore, param::trkLmin, counter)) continue;
                        if(counter<5) counter=5; // tracks that passed the cut on the fit
                    }
                }
                if(counter == 0) mu.push_back(0);
                else if(counter==1) {
                    mu.push_back(0);
                    mu.push_back(1);
                }
                else if(counter==2) {
                    mu.push_back(0);
                    mu.push_back(1);
                    mu.push_back(2);
                }
                else if(counter==3) {
                    mu.push_back(0);
                    mu.push_back(1);
                    mu.push_back(2);
                    mu.push_back(3);
                }
                else if(counter==4) {
                    mu.push_back(0);
                    mu.push_back(1);
                    mu.push_back(2);
                    mu.push_back(3);
                    mu.push_back(4);
                }
                else {
                    mu.push_back(0);
                    mu.push_back(1);
                    mu.push_back(2);
                    mu.push_back(3);
                    mu.push_back(4);
                    mu.push_back(5);
                }
            }
        }
        return mu;
    });

    

} // namespace ana

#endif // STITCH_H