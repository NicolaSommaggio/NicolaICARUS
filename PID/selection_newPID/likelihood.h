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

#include "../selection1muNp_PIDfinale/helper_newPID.h"
#include "helper_generic.h"

using namespace ana;


void write(std::ofstream& file, int clas, std::vector<double> &lr, double depE, double &daughter_depE, double &daughter_angle, int &plane, double &ek_range_mu, double &ek_range_pi, double &ek_range_pro, double &ek_calo, double &ek_true)
{
    file << clas << " ";
    for(int i=0; i<(int)lr.size(); i++)
    {
        file << lr[i] << " ";
    }
    file << depE << " " << daughter_depE << " " << daughter_angle << " " << ek_range_mu << " " << ek_range_pi << " " << ek_range_pro << " " << plane << " " << ek_calo << " " << ek_true << endl;
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

ofstream dump_lkl_ratios_gbdt_withdepE("/exp/icarus/data/users/nsommagg/LIKELIHOOD_FOR_TRAINING_NEW_MC.txt");

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

  int slice_number = -1;
  for(auto const &islc : sr->slc)
  {
 
    slice_number ++;
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    
    if (islc.tmatch.eff <= 0.5 || (vertex_true-vertex_reco).Mag() >= 100 || islc.is_clear_cosmic)continue;

    int ipfp_mu=-1;
    ipfp_mu=find_truth_muon(islc,10);
    
    for (std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
    {   

        /*
        if(sr->hdr.run == 9732 && sr->hdr.evt==20347 && slice_number == 1)
        {

          TVector3 RecoVtx;
          RecoVtx.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
          TVector3 RecoStart;
          RecoStart.SetXYZ(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z);

          cout << true_selection(sr,islc,ipfp) << " ";
          cout << islc.reco.pfp[ipfp].trk.truth.p.pdg << " ";
          cout << islc.reco.pfp[ipfp].trackScore << " ";
          cout << (RecoVtx-RecoStart).Mag() << " ";
          cout << islc.reco.pfp[ipfp].trk.len << " ";
          cout << isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,5.0) << " ";
          cout << isInContained(islc.reco.pfp[ipfp].trk.start.x,islc.reco.pfp[ipfp].trk.start.y,islc.reco.pfp[ipfp].trk.start.z,5.0) << " ";
          cout << islc.reco.pfp[ipfp].trk.end.x*islc.vertex.x << " ";
          cout << islc.reco.pfp[ipfp].parent_is_primary << endl;
        
        }
        */

        int bestplane = find_best_plane(islc,ipfp);

        if(bestplane == -1)continue; //no valid hits (RR < 25 cm and RR > 1 cm and dEdx < 30 MeV/cm and dEdx > 1 MeV/cm) in neither of the 3 planes

        if(bestplane != 2)continue; //<-- to look only collection

        std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
        double dep_E = compute_depE(islc,ipfp,bestplane);
        std::vector<double> dvars = compute_daughter_vars(islc,ipfp);

        double daughter_depE = dvars[0];
        double angle_end = dvars[1];

        double ke_range_mu = -1;
        double ke_range_pi = -1;
        double ke_range_pro = -1;
        double ke_calo = islc.reco.pfp[ipfp].trk.calo[bestplane].ke;
        double ke_true = -1;

        TVector3 p_from_range_mu;
        p_from_range_mu.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.z);
        ke_range_mu = sqrt(pow(105.658,2)+pow(p_from_range_mu.Mag()*1000,2))-105.658;

        TVector3 p_from_range_pi;
        p_from_range_pi.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x, (islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y, (islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);
        ke_range_pi = sqrt(pow(139.570,2)+pow(p_from_range_pi.Mag()*1000,2))-139.570;

        TVector3 p_from_range_pro;
        p_from_range_pro.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
        ke_range_pro = sqrt(pow(938.3,2)+pow(p_from_range_pro.Mag()*1000,2))-938.3;

        double ke_difference_mu_hyp = (ke_range_mu - ke_calo)/ke_calo;
        double ke_difference_pi_hyp = (ke_range_pi - ke_calo)/ke_calo;
        double ke_difference_pro_hyp = (ke_range_pro - ke_calo)/ke_calo;

        if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==13)
        {
          TVector3 true_p;
          true_p.SetXYZ(islc.reco.pfp[ipfp].trk.truth.p.genp.x,islc.reco.pfp[ipfp].trk.truth.p.genp.y,islc.reco.pfp[ipfp].trk.truth.p.genp.z);

          ke_true = sqrt(pow(105.658,2)+pow(true_p.Mag()*1000,2))-105.658;
        }
        if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211)
        {
          TVector3 true_p;
          true_p.SetXYZ(islc.reco.pfp[ipfp].trk.truth.p.genp.x,islc.reco.pfp[ipfp].trk.truth.p.genp.y,islc.reco.pfp[ipfp].trk.truth.p.genp.z);

          ke_true = sqrt(pow(139.570,2)+pow(true_p.Mag()*1000,2))-139.570;
        }
        if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==2212)
        {
          TVector3 true_p;
          true_p.SetXYZ(islc.reco.pfp[ipfp].trk.truth.p.genp.x,islc.reco.pfp[ipfp].trk.truth.p.genp.y,islc.reco.pfp[ipfp].trk.truth.p.genp.z);

          ke_true = sqrt(pow(938.3,2)+pow(true_p.Mag()*1000,2))-938.3;
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

        int pfp_id = id_pfp_truth(islc, ipfp, 10);

        if(int(ipfp)==ipfp_mu)
        //if(is_truth_muon(islc,ipfp,10))
        {
            n_mu ++;
            if(true_selection(sr,islc,ipfp)==0)
            { 
              n_c0 ++;
              write(dump_lkl_ratios_gbdt_withdepE,0,lr,dep_E,daughter_depE,angle_end,bestplane,ke_difference_mu_hyp,ke_difference_pi_hyp,ke_difference_pro_hyp,ke_calo,ke_true);
              //write_dedx(dump_dedx_range,0,rr,dedx,pitch,dirx,diry,dirz);
            }
            if(true_selection(sr,islc,ipfp)==1)
            {
              n_c1 ++;
              write(dump_lkl_ratios_gbdt_withdepE,1,lr,dep_E,daughter_depE,angle_end,bestplane,ke_difference_mu_hyp,ke_difference_pi_hyp,ke_difference_pro_hyp,ke_calo,ke_true);
              //write_dedx(dump_dedx_range,1,rr,dedx,pitch,dirx,diry,dirz);
            }
            //if(true_selection(sr,islc,ipfp)==-1)
            //{
            //  cout << "muon class -1 " << islc.reco.pfp[ipfp].trk.truth.p.pdg << " " << islc.reco.pfp[ipfp].trk.len << " " << momentum_module << endl;
            //}
            continue;
        }
        //else if(pfp_id==1)
        if(pfp_id==1)
        {
            n_pro ++;
            if(true_selection(sr,islc,ipfp)==2)
            {
              n_c2 ++;
              write(dump_lkl_ratios_gbdt_withdepE,2,lr,dep_E,daughter_depE,angle_end,bestplane,ke_difference_mu_hyp,ke_difference_pi_hyp,ke_difference_pro_hyp,ke_calo,ke_true);
              //write_dedx(dump_dedx_range,2,rr,dedx,pitch,dirx,diry,dirz);
            }
            if(true_selection(sr,islc,ipfp)==3)
            {
              n_c3 ++;
              write(dump_lkl_ratios_gbdt_withdepE,3,lr,dep_E,daughter_depE,angle_end,bestplane,ke_difference_mu_hyp,ke_difference_pi_hyp,ke_difference_pro_hyp,ke_calo,ke_true);
              //write_dedx(dump_dedx_range,3,rr,dedx,pitch,dirx,diry,dirz);
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
              n_c4 ++;
              write(dump_lkl_ratios_gbdt_withdepE,4,lr,dep_E,daughter_depE,angle_end,bestplane,ke_difference_mu_hyp,ke_difference_pi_hyp,ke_difference_pro_hyp,ke_calo,ke_true);
              //write_dedx(dump_dedx_range,4,rr,dedx,pitch,dirx,diry,dirz);
            }
            if(true_selection(sr,islc,ipfp)==5)
            {
              n_c5 ++;
              write(dump_lkl_ratios_gbdt_withdepE,5,lr,dep_E,daughter_depE,angle_end,bestplane,ke_difference_mu_hyp,ke_difference_pi_hyp,ke_difference_pro_hyp,ke_calo,ke_true);
              //write_dedx(dump_dedx_range,5,rr,dedx,pitch,dirx,diry,dirz);
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


//5901 14802 1897

//3388 muon rising
//1540 muon mip
//10132 proton rising
//2642 proton interacting
//455 pion rising
//1309 pion interacting

//5843 12719 1788

//3377 muon rising
//1523 muon mip
//8963 proton rising
//2439 proton interacting
//433 pion rising
//1251 pion interacting