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

#include "helper_generic.h"
#include "../selection1muNp_PIDfinale/helper_newPID.h"

using namespace ana;


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
    ipfp_mu=find_truth_muon(islc,10);
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

        if(id_pfp_truth(islc, ipfp, 10)==1)
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

        if(id_pfp_truth(islc, ipfp, 10)==2)
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