#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/Var.h"


#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

//#include "helper_eff.h"
#include "helper_variables.h"

#include "TDirectory.h"
#include "TFile.h"

#include <string>
#include <utility>
#include <vector>

using namespace ana;

//ofstream MyFile("offbeam_run4.txt");


void make_tree_eff_sel()
{


//const std::string fdata = "/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_Overlays_BNB_MC_RUN2/summer_2025/v10_06_00_04p04/flatcaf/*/*/overlay_*flat.caf*.root";

//const std::string fdata = "/pnfs/icarus/persistent/users/jlarkin/cafs/icarus_numudis_cv_*.flat.caf.root";

//const std::string fdata = "/pnfs/icarus/persistent/users/rtriozzi/nugraph/nugraphreco_MoreVars/numinom_noyzsim_2_NuGraphReco_HIPTagger.flat.caf.root";

//const std::string fdata_unblind = "/pnfs/sbn/data/sbn_fd/poms_production/data/Reproc_Run2_SBN_wCalib/reconstructed/icaruscode_v10_06_00_04p03/bnbmajority/flatcaf_unblind/*/*/*.root";


 //SpectrumLoader data("Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_wCalib_v10_06_00_04p03_bnbmajority_flatcaf_unblind");
 //SpectrumLoader data("production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_06p01_flatcaf");
  //const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_*.flat.caf.root";

  //const std::string fdata = "/pnfs/icarus/scratch/users/petrillo/jobOutput/NeutrinoOverlayFix/run9583/cafmakerjob_icarus_detsim2d_overlay/out/*/*.flat.caf.root";

//ICARUS MC
  const std::string fdata = "/pnfs/icarus/scratch/users/gputnam/Ar23+_iterE/ICARUSSpringMC/*/*flat.caf.root";
  const std::string fdata_nicola = "/exp/icarus/data/users/nsommagg/*/84005943*_out1.flat.caf.root";


//ICARUS Offbeam
  const std::string fdata_offbeam = "Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_v10_06_00_01p05_offbeambnbmajority_flatcaf_unblind";

//ICARUS Dirt
  const std::string fdata_dirt = "/pnfs/sbnd/scratch/users/twester/run4_dirt_test/*/*/*flat.caf.root";

//ICARUS Prescaled
  const std::string fdata_prescaled = "Icaruspro_2025_wcdnn_production_Reproc_Run2_SBN_2_v10_06_00_06p03_bnbmajority_flatcaf_prescaled";

//ICARUS WireMod θxw-X Variation
  const std::string fdata_WM = "/pnfs/sbn/scratch/users/twester/Run2_WM_ReCAF2026/flatcaf/*/*/*.flat.caf.root";

//ICARUS SCE/Field Variation
  const std::string fdata_SCE = "/pnfs/sbn/data/sbn_fd/aurora/mc/v10_06_00_13/run2_sce/caf/*/*/*.flat.caf.root";

//Sample prescaled fixed light 
  const std::string fdata_fixedlight = "/pnfs/icarus/persistent/users/fpoppi/LightRun2/data/data/job_*/*.Prescaled.OKTOLOOK.flat.caf.root";

//sample MC fixed light 
  const std::string fdata_fixedlightMC = "/pnfs/icarus/persistent/users/fpoppi/LightMC/mc/job_*/*.flat.caf.root";

//sample offbeam unblind fixed light 
  const std::string fdata_fixedlightoff = "/pnfs/icarus/persistent/users/fpoppi/LightRun2/offbeam/offbeam/job_*/*.Unblind.DONOTLOOK.flat.caf.root";

// New Run 2 official 
  const std::string fMC_run2 = "/pnfs/sbn/scratch/users/twester/Run2_ReCAF2026/flatcaf/0000*[5,6,7,8,9]/*/*.flat.caf.root";


//==========================================================================================================================

//Run 4 Prescaled
  const std::string fdata_prescaledrun4 = "/pnfs/sbn/persistent/users/gputnam/icarus-recaf-2026/Run4BeamOnWF/*/*.flat.caf.root";

//Offbeam 20% Run4
  const std::string fdata_offbeamrun4 = "/pnfs/sbn/persistent/users/gputnam/icarus-recaf-2026/Run4BeamOffWF/70479293_*/*.flat.caf.root";

//Dirt run 4 
  const std::string fdata_dirtrun4 = "/pnfs/sbn/scratch/users/twester/Run4_Dirt_ReCAF2026/flatcaf/*/*/*.flat.caf.root";

//MC Overlays run4
  const std::string fdata_MCrun4 = "/pnfs/sbn/scratch/users/twester/Run4_ReCAF2026/flatcaf/*/*/*.flat.caf.root";

//MAPLE-GUMP-STUDIES RUN2
  const std::string fdata_maple_gump_RUN2 = "/pnfs/sbn/scratch/users/twester/Run2_ReCAF2026/flatcaf/000000/00000*/*.flat.caf.root";

//MAPLE-GUMP-STUDIES RUN4
  const std::string fdata_maple_gump_RUN4 = "/pnfs/sbn/scratch/users/twester/Run4_ReCAF2026/flatcaf/00007*/*/*.flat.caf.root";

//MAPLE-GUMP-STUDIES RUN4 SMALLER
  const std::string fdata_maple_gump_RUN4_smaller = "/pnfs/sbn/scratch/users/twester/Run4_ReCAF2026/flatcaf/000071/*/*.flat.caf.root";

//MAPLE-GUMP-STUDIES RUN4 CNAF
  const std::string fdata_maple_gump_RUN4_cnaf = "/storage/gpfs_data/icarus/local/users/cfarnese/files_test_Nicola_RUN4*/file*.flat.caf.root";

//MAPLE-GUMP-STUDIES RUN4 CNAF OFFBEAM
  const std::string fdata_maple_gump_RUN4_cnaf_offbeam = "/storage/gpfs_data/icarus/local/users/cfarnese/files_test_Nicola_RUN4_offbeam/offbeam_4.flat.caf.root";

//MAPLE-GUMP-STUDIES RUN4 CONCAT
  const std::string fdata_maple_gump_RUN4_concat = "/pnfs/icarus/persistent/users/cfarnese/for_Nicola/file*.flat.caf.root";

//MAPLE-GUMP-STUDIES RUN4 CONCAT OFFBEAM 
  const std::string fdata_maple_gump_RUN4_concat_offbeam = "/pnfs/icarus/persistent/users/cfarnese/for_Nicola/offbeam_4.flat.caf.root";


/*
  SpectrumLoader loader(fdata_maple_gump_RUN4_concat);         

  const Binning kBinz = Binning::Simple(300,0,30);

  Spectrum s1("", kBinz, loader, kPrint_selection ,kNoSpillCut );
  //Spectrum s1("", kBinz, loader, check_variations ,kNoSpillCut );

  loader.Go();

  double factor = s1.POT();  

  TH1D* h1 = s1.ToTH1(factor);  
*/



 SpectrumLoader data(fdata_maple_gump_RUN4_concat);


    const SpillMultiVar True_Enu = kTrue_Enu;
    const SpillMultiVar True_visible_Enu = kTrue_visible_Enu;
    const SpillMultiVar Genie_mode = kGenie_mode;
    const SpillMultiVar Eff_cuts = kEff_cuts;
    const SpillMultiVar Eff_raster_angle = kEff_raster_angle;
    const SpillMultiVar Eff_raster_nuscore = kEff_raster_nuscore;
    const SpillMultiVar Eff_cryosel = kEff_cryosel;

    

    const SpillMultiVar True_protons = kTrue_protons;



    //const SpillMultiVar True_Enu_allnumuCC = kTrue_Enu_allnumuCC;
    //const SpillMultiVar Eff_cuts_allnumuCC = kEff_cuts_allnumuCC;


    
    
    //const SpillMultiVar Eff_cuts_TEST = kEff_cuts_TEST;

    //const SpillMultiVar Crthit_distance = kCrthit_distance;
    //const SpillMultiVar Crthit_sys = kCrthit_distance;

    
    const SpillMultiVar Reco_energy = kReco_energy;
    const SpillMultiVar Reco_class = kReco_class;
    const SpillMultiVar T3D_angle_mup = kT3D_angle_mup;


    const SpillMultiVar Muon_endx = kMuon_endx;
    const SpillMultiVar Muon_endy = kMuon_endy;
    const SpillMultiVar Muon_endz = kMuon_endz;

    const SpillMultiVar Proton_endx = kProton_endx;
    const SpillMultiVar Proton_endy = kProton_endy;
    const SpillMultiVar Proton_endz = kProton_endz;


    //const SpillMultiVar Count_entries = kCount_entries;


    const SpillMultiVar Vertex_x = kVertex_x;
    const SpillMultiVar Vertex_y = kVertex_y;
    const SpillMultiVar Vertex_z = kVertex_z;

    
    
    const SpillMultiVar Muon_length = kMuon_length;
    const SpillMultiVar Proton_length_leading = kProton_length_leading;
    const SpillMultiVar Proton_kinetic_leading = kProton_kinetic_leading;
    const SpillMultiVar Proton_chi2mu = kProton_chi2mu;
    const SpillMultiVar Proton_chi2pro = kProton_chi2pro;
    const SpillMultiVar Muon_chi2mu = kMuon_chi2mu;
    const SpillMultiVar Muon_chi2pro = kMuon_chi2pro;



    const SpillMultiVar Transverse_angle = kTransverse_angle;
    const SpillMultiVar Transverse_mom_reco = kTransverse_mom_reco;
    const SpillMultiVar Transverse_mom_reco_mu = kTransverse_mom_reco_mu;
    const SpillMultiVar Transverse_mom_reco_pro = kTransverse_mom_reco_pro;
    const SpillMultiVar Number_protons = kNumber_protons;
    const SpillMultiVar Barycenter_delta = kBarycenter_delta;

    const SpillMultiVar Muon_trackScore = kMuon_trackScore;
    const SpillMultiVar Proton_trackScore = kProton_trackScore;


    const SpillMultiVar Nu_score = kNu_score;
    const SpillMultiVar CryoSel = kCryoSel;

    const SpillMultiVar E_nu_true = k_E_nu_true;
    const SpillMultiVar slice_index = k_slice_indedx;
    const SpillMultiVar flashPE = k_pe;

  
    //const SpillMultiVar Reco_energy_resolution = kReco_energy_resolution;
    //const SpillMultiVar Reco_tmatch_eff = kReco_tmatch_eff;

    //const SpillMultiVar Barycenter_delta = kBarycenter_delta;
    //const SpillMultiVar True_visible_Enu_reco = kTrue_visible_Enu_reco;

    

    
    //const SpillCut kSpillSelection = kSelectRun_M == 9435 && kSelectEvt_2d;


    const SpillCut kSpillSelection = kNoSpillCut;

    Tree nutree(
      "selectedNu", 
      {
        "true_Enu",
        "True_visible_Enu", 
        "Genie_mode",
        "Pass_cut",
        "Eff_raster_angle",
        "True_protons",
        "Eff_raster_nuscore",
        "Eff_cryosel"
      }, 
      data, 
      {
        True_Enu,
        True_visible_Enu,
        Genie_mode,
        Eff_cuts,
        Eff_raster_angle,
        True_protons,
        Eff_raster_nuscore,
        Eff_cryosel
      }, 
      kSpillSelection, 
      true
    );
    //Tree nutree_reco("selectedReco", {"recoE","Reco_class","T3D_angle_mup","Muon_endx","Muon_endy","Muon_endz","Proton_endx","Proton_endy","Proton_endz"}, data, {Reco_energy,Reco_class,T3D_angle_mup,Muon_endx,Muon_endy,Muon_endz,Proton_endx,Proton_endy,Proton_endz}, kSpillSelection, true);
    //Tree counttree("selectedReco2", {"Count_entries"}, data, {Count_entries}, kSpillSelection, true);

    Tree nutree_reco(
      "selectedReco", 
      {
        "recoE",
        "Reco_class",
        "Muon_length",
        "Proton_length_leading",
        "Proton_kinetic_leading",
        "Proton_chi2mu",
        "Proton_chi2pro",
        "Muon_chi2mu",
        "Muon_chi2pro",
        "Vertex_x",
        "Vertex_y",
        "Vertex_z",
        "Muon_endx",
        "Muon_endy",
        "Muon_endz",
        "Proton_endx",
        "Proton_endy",
        "Proton_endz",
        "Transverse_angle",
        "T3D_angle_mup",
        "deltaPt",
        "Transverse_mom_reco_mu",
        "Transverse_mom_reco_pro",
        "Number_protons",
        "Barycenter_delta",
        "Muon_trackScore",
        "Proton_trackScore",
        "Nu_score",
        "CryoSel",
        "E_nu_true",
        "Slice_index",
        "FlashPE"
      }, 
      data, 
      {
        Reco_energy,
        Reco_class,
        Muon_length,
        Proton_length_leading,
        Proton_kinetic_leading,
        Proton_chi2mu,
        Proton_chi2pro,
        Muon_chi2mu,
        Muon_chi2pro,
        Vertex_x,
        Vertex_y,
        Vertex_z,
        Muon_endx,
        Muon_endy,
        Muon_endz,
        Proton_endx,
        Proton_endy,
        Proton_endz,
        Transverse_angle,
        T3D_angle_mup,
        Transverse_mom_reco,
        Transverse_mom_reco_mu,
        Transverse_mom_reco_pro,
        Number_protons,
        Barycenter_delta,
        Muon_trackScore,
        Proton_trackScore,
        Nu_score,
        CryoSel,
        E_nu_true,
        slice_index,
        flashPE
      }, 
      kSpillSelection, 
      true
    );



    //const Binning kBinPot     = Binning::Simple(1000,0,1.);
    //Spectrum spot("Pot",kBinPot,data,gate_delta, kNoSpillCut);

  data.Go();

//MyFile << " Gate_delta " << final_gate_delta << " Livetime " << final_livetime << endl;



//  TFile fout("/exp/sbnd/data/users/marterop/sel_dirt/selected_icarus_standard_eff_sel_MC_tmatch05_dirt_nocrt.root", "RECREATE");
  //TFile fout("/exp/sbnd/data/users/marterop/sel_dirt/selected_icarus_standard_eff_sel_MC_tmatch05_nocrt_nolight_MC_vars_dirt.root", "RECREATE");
  //TFile fout("selected_icarus_standard_eff_sel_MC_tmatch05_nocrt_nolight_MC_vars_offbeam_nuscore.root", "RECREATE");

  TFile fout("sbruce_trees/CUT_TEST/MC_CONCAT/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_50_CHI2_VAR_TRACKSCORE_VAR_NOVISIBLE_VAR_CRTveto_Lmu_FV_Trigger_GUMP_FINAL_MOD.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_50_PRO_KE_50_CHI2_VAR_TRACKSCORE_VAR_NOVISIBLE_VAR_CRTveto_NOLmu_FV_Trigger.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR_TRACKSCORE_VAR_NOVISIBLE_VAR_CRTveto_Lmu_FV_Trigger.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_50_PRO_KE_50_CHI2_VAR_TRACKSCORE_VAR_VISIBLE_VAR_CRTveto_NOLmu_FV_Trigger.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_50_PRO_KE_50_CHI2_VAR_TRACKSCORE_VAR_VISIBLE_VAR_CRTveto_Lmu_FV_Trigger.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR_TRACKSCORE_VAR_VISIBLE_VAR_CRTveto_Lmu_FV_Trigger.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR_TRACKSCORE_VAR_VISIBLE_VAR_CRTveto_Lmu_FV.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR_TRACKSCORE_VAR_VISIBLE_VAR_CRTveto_Lmu.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR_TRACKSCORE_VAR_VISIBLE_VAR_CRTveto.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR_TRACKSCORE_VAR_VISIBLE_VAR_newvars.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR_TRACKSCORE_VAR.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40_CHI2_VAR.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_MU_L_40_PRO_KE_40.root", "RECREATE");
  //TFile fout("sbruce_trees/SBRUCE_TREE_MAPLE_NICOLA_RUN4_newvars.root", "RECREATE");

//  TFile fout("/exp/sbnd/data/users/marterop/sel_dirt/selected_icarus_standard_eff_sel_MC_tmatch05_dirt.root", "RECREATE");

  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  nutree_reco.SaveTo(dir); 
  //counttree.SaveTo(dir); 


}

