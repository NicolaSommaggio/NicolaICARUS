#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/SBNOnOffSysts.h"
#include "sbnana/CAFAna/Systs/UniverseOracle.h"
#include "sbnana/CAFAna/Systs/IcarusRun2DetectorSysts.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202401.h"
#include "sbnana/SBNAna/Cuts/ICARUSDataQualityCuts.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/SBNAna/Vars/Vars.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "helper_icarus_variables.h"

#include "TDirectory.h"
#include "TFile.h"

#include <string>
#include <utility>
#include <vector>

using namespace ana;


// --------------------------------------------------------------------------
// A systematic whose weight at each integer sigma point is given explicitly,
// rather than being pulled from GENIE universes.
class SBNFixedValueSyst : public ISyst
{
public:
  SBNFixedValueSyst(const std::string& systName,
                     const std::map<int, double>& sigmaToWeight)
    : ISyst(systName, systName),
      fWeights(sigmaToWeight)
  {
  }

  void Shift(double x, caf::SRSliceProxy* sr, double& weight) const override
  {
    if(sr->truth.index < 0) return;
    weight *= GetWeight(x);
  }

  void Shift(double x, caf::SRTrueInteractionProxy* nu, double& weight) const override
  {
    if(nu->index < 0) return;
    weight *= GetWeight(x);
  }

protected:
  double GetWeight(double x) const
  {
    const int isig = std::lround(x);

    auto it = fWeights.find(isig);
    if(it == fWeights.end()){
      static bool once = true;
      if(once){
        once = false;
        std::cout << "SBNFixedValueSyst: WARNING no weight defined for "
                  << ShortName() << " at sigma = " << isig
                  << ". Returning weight 1." << std::endl;
      }
      return 1;
    }

    return it->second;
  }

  std::map<int, double> fWeights;
};

//void make_tree_newAr_SBND(std::string outname = "/exp/sbnd/data/users/marterop/Ar23+_iterE/SBND_sbnana_fix_.root")
//void make_tree_newAr_SBND(std::string outname = "/exp/sbnd/data/users/marterop/Ar23+_iterE/test_nolight/SBND_test_bary_14_vars.root")
//void make_tree_newAr_SBND(std::string outname = "/exp/sbnd/data/users/marterop/Ar23+_iterE/no_crossers/SBND_syst_nocrosser_9_vars.root")
void make_tree_newAr(std::string outname = "/exp/icarus/app/users/marterop/dev_areas/sbn_syst/syst_macros/nolight_nocrt/ICARUS_syst_vars_newxsec_9.root")

{
  //SpectrumLoader mc("/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_Overlays_BNB_MC_RUN2/summer_2025/v10_06_00_06p01/flatcaf/60/*/*.root");
  //SpectrumLoader mc("/pnfs/sbnd/scratch/users/gputnam/Ar23+_iterE/SBNDSpringMC-2xjobs/27536650_[3]*/*flat.caf.root");

  //SpectrumLoader mc("/pnfs/sbnd/scratch/users/mueller/AR23p/20260706/mcp202b_bnbrockbox_1e20/28945606_8*/*.ar23p.flat.caf.root");
  //SpectrumLoader mc("/pnfs/sbn/persistent/users/mueller/osc/sbnd/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_sbnd/*/9*/*.ar23p.flat.caf.root");

//SBND official MC 
//  SpectrumLoader mc("/pnfs/sbn/data_add/sbn_nd/poms_production/mc/SBND2026A/v1_00_01/Ar23plus_cafRepro_BNBLight/flatcaf/6*/*.flat.caf.root");

  SpectrumLoader mc("/pnfs/sbn/scratch/users/twester/Run2_ReCAF2026/flatcaf/0000*9/*/*.flat.caf.root");





//SBND DIRT
  //SpectrumLoader mc("/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025B_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd/CV/caf/[8,9]*/*/caf.flat.caf*.root");

  //SpectrumLoader mc("/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_*/*.root");
  //SpectrumLoader mc("production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_06p01_flatcaf");

  
  //SpectrumLoader offbeam("Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_v10_06_00_01p05_offbeambnbmajority_flatcaf_unblind");
  //Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_v10_06_00_01p05_offbeambnbmajority_flatcaf_unblind


  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kTrueL = SIMPLEVAR(truth.baseline);
  const Var kTruePDG = SIMPLEVAR(truth.pdg);
  const Var kTrueCC = SIMPLEVAR(truth.iscc);
  const Var kSlcVX = SIMPLEVAR(vertex.x);
  const Var kSlcVY = SIMPLEVAR(vertex.y);
  const Var kSlcVZ = SIMPLEVAR(vertex.z);

  const Cut kTrueNu = SIMPLEVAR(truth.index) >= 0;
  const SpillCut kSpillSelection = kNoSpillCut;

  const Cut kSliceSelection = Selection_1muNp_slice; //kNoCut


  std::vector<std::string> nu_branch_names = {
    "trueE", "trueL", "truePDG", "CC", "recoE","deltaPt","Muon_length","Muon_endx","Muon_endy","Muon_endz","Proton_length_leading","Proton_kinetic_leading","Proton_chi2mu","Proton_chi2pro",
    "Muon_chi2mu","Muon_chi2pro","Vertex_x","Vertex_y","Vertex_z","Proton_endx","Proton_endy","Proton_endz","Transverse_angle","T3D_angle_mup",
    "Transverse_mom_reco_mu","Transverse_mom_reco_pro","Number_protons","Muon_trackScore","Proton_trackScore",
  };

  std::vector<Var> nu_vars = {
    kTrueE, kTrueL, kTruePDG, kTrueCC,kReco_energy,kTransverse_mom_reco,kMuon_length,kMuon_endx,kMuon_endy,kMuon_endz,kProton_length_leading,kProton_kinetic_leading,kProton_chi2mu,kProton_chi2pro,
    kMuon_chi2mu,kMuon_chi2pro,kVertex_x,kVertex_y,kVertex_z,kProton_endx,kProton_endy,kProton_endz,kTransverse_angle,kT3D_angle_mup,
    kTransverse_mom_reco_mu,kTransverse_mom_reco_pro,kNumber_protons,kMuon_trackScore,kProton_trackScore,
  };

  std::vector<std::string> branch_names = {
    "recoE","deltaPt","Muon_length","Muon_endx","Muon_endy","Muon_endz","Proton_length_leading","Proton_kinetic_leading","Proton_chi2mu","Proton_chi2pro",
    "Muon_chi2mu","Muon_chi2pro","Vertex_x","Vertex_y","Vertex_z","Proton_endx","Proton_endy","Proton_endz","Transverse_angle","T3D_angle_mup",
    "Transverse_mom_reco_mu","Transverse_mom_reco_pro","Number_protons","Muon_trackScore","Proton_trackScore",
  };

  std::vector<Var> vars = {
    kReco_energy,kTransverse_mom_reco,kMuon_length,kMuon_endx,kMuon_endy,kMuon_endz,kProton_length_leading,kProton_kinetic_leading,kProton_chi2mu,kProton_chi2pro,
    kMuon_chi2mu,kMuon_chi2pro,kVertex_x,kVertex_y,kVertex_z,kProton_endx,kProton_endy,kProton_endz,kTransverse_angle,kT3D_angle_mup,
    kTransverse_mom_reco_mu,kTransverse_mom_reco_pro,kNumber_protons,kMuon_trackScore,kProton_trackScore,
  };

  Spectrum dummy_spec2("", Binning::Simple(2,0,2), mc, kFillVars, kNoSpillCut); 
  Tree nutree("selectedNu", nu_branch_names, mc, nu_vars, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  Tree costree("selectedCos", branch_names, mc, vars, kSpillSelection, kSliceSelection && !kTrueNu, kNoShift, true, true);
  //Tree offbeamtree("selectedOffbeam", branch_names, offbeam, vars, kSpillSelection, kSliceSelection, kNoShift, true, true);
  //Spectrum dummy_spec("", Binning::Simple(2,0,2), offbeam, kOffbeamLivetime, kNoSpillCut); 


  std::vector<std::string> genie_names;

  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_CoulombCCQE");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_NormCCMEC");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_NormNCMEC");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_DecayAngMEC");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad");

  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_MaCCRES");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_MaNCRES");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_MvCCRES");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_MvNCRES");

  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_RDecBR1eta");

  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_MFP_pi");


  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_FrCEx_pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_FrInel_pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_FrAbs_pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi");

  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi");

  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_AhtBY");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_BhtBY");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_CV1uBY");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_CV2uBY");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_NormCCCOH");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_NormNCCOH");

  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_MaNCEL");
  genie_names.push_back("GENIEReWeight_SBN_v1_multisigma_EtaNCEL");

  genie_names.push_back("ZExpPCAWeighter_SBN_v3_MvA_b1");
  genie_names.push_back("ZExpPCAWeighter_SBN_v3_MvA_b3");
  genie_names.push_back("ZExpPCAWeighter_SBN_v3_MvA_b2");
  genie_names.push_back("ZExpPCAWeighter_SBN_v3_MvA_b4");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToHF_q0bin0");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToHF_q0bin1");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToHF_q0bin2");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToHF_q0bin3");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToHF_q0bin4");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToSF_q0bin0");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToSF_q0bin1");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToSF_q0bin2");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToSF_q0bin3");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_LFGToSF_q0bin4");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin0");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin1");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin2");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin3");
  genie_names.push_back("CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin4");
  genie_names.push_back("QEInterference_SBN_v3_QEIntf_dial_0");
  genie_names.push_back("QEInterference_SBN_v3_QEIntf_dial_1");
  genie_names.push_back("QEInterference_SBN_v3_QEIntf_dial_2");
  genie_names.push_back("QEInterference_SBN_v3_QEIntf_dial_3");
  genie_names.push_back("QEInterference_SBN_v3_QEIntf_dial_4");
  genie_names.push_back("QEInterference_SBN_v3_QEIntf_dial_5");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrG4_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrINCL_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrG4LoE_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrINCLLoE_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrG4M1E_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrINCLM1E_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrG4M2E_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrINCLM2E_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrG4HiE_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrINCLHiE_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_MFPLoE_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_MFPM1E_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_MFPM2E_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_MFPHiE_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrKin_PiProFix_N");
  genie_names.push_back("GENIEReWeight_SBN_v3_FrKin_PiProBias_N");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin0");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin1");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin2");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin3");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin0");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin1");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin2");
  genie_names.push_back("MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin3");

  genie_names.push_back("PionAbsWeighter_SBN_v3_QuasiDeuteronFraction");
  genie_names.push_back("CCQEXSecCorr_SBN_v3_CCQEXSecCorr");

  
  std::vector<const ISyst*> genie_systs;

for (const auto& name : genie_names) {
    genie_systs.push_back(new SBNWeightSyst(name));
}

  std::vector<std::string> nsigma_names = genie_names;
  std::vector<const ISyst*> nsigma_systs = genie_systs;
  std::vector<std::pair<int, int>> min_max;
  for(size_t i = 0; i < genie_names.size(); ++i) min_max.emplace_back(-3, 3);


  std::map<int, double> myvals = {
        {-3, 0.982714}, {-2, 0.9887274}, {-1, 0.99474195},
        { 0, 1.00},
        { 1, 1.005}, { 2, 1.01}, { 3, 1.015}
    };

    auto* myFixedSyst = new SBNFixedValueSyst("new_POT_norm", myvals);

    genie_systs.push_back(myFixedSyst);
    nsigma_names.push_back("new_POT_norm");
    nsigma_systs.push_back(myFixedSyst);
    min_max.emplace_back(-3, 3);


  NSigmasTree nsigtree("multisigmaTree", nsigma_names, mc, nsigma_systs, min_max, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);

  const std::vector<std::string> flux_names{ "expskin_Flux", "horncurrent_Flux", "nucleoninexsec_Flux", "nucleonqexsec_Flux", "nucleontotxsec_Flux", "pioninexsec_Flux", "pionqexsec_Flux", "piontotxsec_Flux", "piplus_Flux", "piminus_Flux", "kplus_Flux", "kminus_Flux", "kzero_Flux" };

  std::vector<std::string> multisim_names;
  std::vector<std::vector<Var>> univsKnobs;
  std::vector<unsigned int> nuniverses;
  for(const auto& name: flux_names) {
    multisim_names.push_back(name);
    size_t nuniv = 1000;
    nuniverses.push_back(nuniv);
    univsKnobs.emplace_back();
    for(size_t i = 0; i < nuniv; ++i) 
      univsKnobs.back().push_back(GetUniverseWeight(name, i));
  }


  NUniversesTree nunivtree("multisimTree", multisim_names, mc, univsKnobs, nuniverses, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);

  mc.Go();
  //offbeam.Go();

  //offbeamtree.OverrideLivetime(offbeam_livetime);


  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  costree.SaveTo(dir);
  nsigtree.SaveTo(dir);
  nunivtree.SaveTo(dir);
  //detsysttree.SaveTo(dir);
  //TDirectory *offdir = fout.mkdir("offbeam");
  //offbeamtree.SaveTo(offdir);
}
