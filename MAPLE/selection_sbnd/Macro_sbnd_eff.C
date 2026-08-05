//#include "helper_new.h"

//Run from: /exp/icarus/app/users/marterop/dev_areas/sbn_ana/sbnd_tests


#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/Var.h"


#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

//#include "helper_eff.h"
#include "helper_variables_sbnd.h"

#include "TDirectory.h"
#include "TFile.h"

#include <string>
#include <utility>
#include <vector>

using namespace ana;



#include <set>
#include <tuple>

const SpillMultiVar n_record_spill_data([](const caf::SRSpillProxy* spill)-> std::vector<double>
{
    static std::set<std::tuple<unsigned int, unsigned int, unsigned int>> seen;

    std::vector<double> vector_active;

    unsigned int run    = spill->hdr.run;
    unsigned int subrun = spill->hdr.subrun;
    unsigned int event  = spill->hdr.evt;

    std::tuple<unsigned int, unsigned int, unsigned int> key = std::make_tuple(run, subrun, event);

    if (seen.insert(key).second) {
        vector_active.push_back(1.0); // new unique (run, subrun, evt)
    } else {
        vector_active.push_back(0.0); // already seen
    }

    return vector_active;
});


const SpillMultiVar n_pot([](const caf::SRSpillProxy* spill)-> std::vector<double>
{
    static std::set<std::tuple<unsigned int, unsigned int, unsigned int>> seen;

    std::vector<double> vector_active;

    unsigned int run    = spill->hdr.run;
    unsigned int subrun = spill->hdr.subrun;
    unsigned int event  = spill->hdr.evt;

    std::tuple<unsigned int, unsigned int, unsigned int> key = std::make_tuple(run, subrun, event);

    for(const auto &spillinfo: spill->hdr.bnbinfo) {
        if(spillinfo.TOR875>0){
            vector_active.push_back(1.0);
        }
        else {
        vector_active.push_back(0.0); // already seen
        }

        }

    return vector_active;
});

const SpillMultiVar n_noffbeambnb([](const caf::SRSpillProxy* spill)-> std::vector<double>
{
    static std::set<std::tuple<unsigned int, unsigned int, unsigned int>> seen;

    std::vector<double> vector_active;

    unsigned int run    = spill->hdr.run;
    unsigned int subrun = spill->hdr.subrun;
    unsigned int event  = spill->hdr.evt;

    std::tuple<unsigned int, unsigned int, unsigned int> key = std::make_tuple(run, subrun, event);

    if (seen.insert(key).second && spill->hdr.first_in_subrun) {
        vector_active.push_back(spill->hdr.noffbeambnb); // new unique (run, subrun, evt)
    } else {
        vector_active.push_back(0.0); // already seen
    }

    return vector_active;
});


void Macro_sbnd_eff()
{


  //const std::string fdata = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_668X.flat.caf.root";
  //const std::string fdata = "/pnfs/sbnd/scratch/users/sungbino/sbnd_2025_prod/v10_06_00_09/flatcaf/mc_1e20/flatcaf_790.root";
  //const std::string fdata = "/pnfs/icarus/scratch/users/marterop/sbnd_files/flatcaf_*.root";
  //const std::string fdata = "/pnfs/sbn/scratch/users/sungbino/sbnd/v10_06_00_09/mc_1e20/flatcaf_12*.root";

  //const std::string fdata = "/pnfs/sbnd/scratch/users/gputnam/Ar23+_iterE/SBNDSpringMC-2xjobs/27536650_[0,1,2]/*flat.caf.root";
  const std::string fdata = "/pnfs/sbnd/scratch/users/mueller/AR23p/2*/7*/*.ar23p.flat.caf.root";


  const std::string fdata_offbeam = "/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025C/v10_06_00_09/Spring25_reprocess_Intime_offbeamlight/flatcaf/offbeamlight/*/*flat.caf.root";

  const std::string fdata_dirt = "/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025B_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd/CV/caf/*/*/caf.flat.caf-*.root";


  const std::string fdata_sbnd ="/pnfs/sbn/persistent/users/mueller/osc/sbnd/mc_MCP2025C_1e20_v10_06_00_09_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_sbnd/*/[8]*/*.ar23p.flat.caf.root";
  
  const std::string fMC_sbnd ="/pnfs/sbn/data_add/sbn_nd/poms_production/mc/SBND2026A/v1_00_01/Ar23plus_cafRepro_BNBLight/flatcaf/[0]*/*.flat.caf.root";

  //SBND DATA 
    const std::string fdata_run1 ="/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025C/v10_06_00_09/Spring25_reprocess_FixedDev_bnblight/flatcaf/bnblight/*/*.flat.caf.root";



  //const std::string fdata = "/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/6*/*/caf.flat*.root";
  //const std::string fdata = "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbnd/scratch/users/sungbino/sbnd_2025_prod/v10_06_00_09/flatcaf/mc_1e20/flatcaf_1060.root";
  SpectrumLoader data(fdata_sbnd);


    const SpillMultiVar True_Enu = kTrue_Enu;
    const SpillMultiVar True_visible_Enu = kTrue_visible_Enu;
    const SpillMultiVar Genie_mode = kGenie_mode;

    const SpillMultiVar Eff_cuts = kEff_cuts;
    const SpillMultiVar Eff_raster_angle = kEff_raster_angle;

    const SpillMultiVar Reco_energy = kReco_energy;
    const SpillMultiVar Reco_class = kReco_class;
    const SpillMultiVar True_protons = kTrue_protons;


    //const SpillCut kSpillSelection = kSelectRun_M == 9435 /*&& kSelectEvt_2d*/;


    const SpillMultiVar Muon_endx = kMuon_endx;
    const SpillMultiVar Muon_endy = kMuon_endy;
    const SpillMultiVar Muon_endz = kMuon_endz;
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


    const SpillMultiVar Proton_endx = kProton_endx;
    const SpillMultiVar Proton_endy = kProton_endy;
    const SpillMultiVar Proton_endz = kProton_endz;
    const SpillMultiVar Transverse_angle = kTransverse_angle;
    const SpillMultiVar T3D_angle_mup = kT3D_angle_mup;
    const SpillMultiVar Transverse_mom_reco = kTransverse_mom_reco;
    const SpillMultiVar Transverse_mom_reco_mu = kTransverse_mom_reco_mu;
    const SpillMultiVar Transverse_mom_reco_pro = kTransverse_mom_reco_pro;
    const SpillMultiVar Number_protons = kNumber_protons;

    const SpillMultiVar Muon_trackScore = kMuon_trackScore;
    const SpillMultiVar Proton_trackScore = kProton_trackScore;
    const SpillMultiVar Flashmatch_deltaZ = kdeltaZ;
    const SpillMultiVar Nu_score = kNu_score;

    
    


    const SpillCut kSpillSelection = kNoSpillCut;

    const SpillMultiVar Num_record_spill_data = n_record_spill_data;
    const SpillMultiVar Num_pot = n_pot;
    const SpillMultiVar Num_noffbeambnb = n_noffbeambnb;

    const SpillMultiVar True_Enu_allnumuCC = kTrue_Enu_allnumuCC;
    const SpillMultiVar True_Enu_allnumuCC_nocrossers = kTrue_Enu_allnumuCC_nocrossers;


    const SpillMultiVar Eff_cuts_allnumuCC = kEff_cuts_allnumuCC;





    //Tree nutree("selectedNumuCC", {"True_Enu_allnumuCC","True_Enu_allnumuCC_nocrossers","Eff_cuts_allnumuCC"}, data, {True_Enu_allnumuCC,True_Enu_allnumuCC_nocrossers,Eff_cuts_allnumuCC}, kSpillSelection, true);

    //Tree nutree_reco("selectedReco", {"Muon_endx_all","Muon_endy_all","Muon_endz_all","Muon_endx_all_true","Muon_endy_all_true","Muon_endz_all_true"}, data, {Muon_endx_all,Muon_endy_all,Muon_endz_all,Muon_endx_all_true,Muon_endy_all_true,Muon_endz_all_true}, kSpillSelection, true);
    //Tree nutree_reco3("selectedReco3", {"Num_record_spill_data","Num_noffbeambnb"}, data, {Num_record_spill_data,Num_noffbeambnb}, kSpillSelection, true);
    //Tree nutree_reco2("selectedReco2", {"Num_pot"}, data, {Num_pot}, kSpillSelection, true);

//This one!
    Tree nutree("selectedNu", {"true_Enu","True_visible_Enu", "Genie_mode","Pass_cut","Eff_raster_angle","True_protons"}, data, {True_Enu,True_visible_Enu,Genie_mode,Eff_cuts,Eff_raster_angle,True_protons}, kSpillSelection, true);

    //Tree nutree_reco("selectedReco", {"Reco_energy","Reco_class","Reco_energy_resolution","Reco_tmatch_eff","Barycenter_delta","True_visible_Enu_reco","Crthit_distance","Crthit_sys"}, data, {Reco_energy,Reco_class,Reco_energy_resolution,Reco_tmatch_eff,Barycenter_delta,True_visible_Enu_reco,Crthit_distance,Crthit_sys}, kSpillSelection, true);
    //Tree nutree_reco("selectedReco", {"Reco_energy"}, data, {Reco_energy}, kSpillSelection, true);
    //Tree nutree_reco("selectedReco", {"recoE","Reco_class"}, data, {Reco_energy,Reco_class}, kSpillSelection, true);

    Tree nutree_reco("selectedReco", {"recoE","Reco_class","Muon_length","Proton_length_leading","Proton_kinetic_leading","Proton_chi2mu","Proton_chi2pro","Muon_chi2mu","Muon_chi2pro","Vertex_x","Vertex_y","Vertex_z","Muon_endx","Muon_endy","Muon_endz","Proton_endx","Proton_endy","Proton_endz","Transverse_angle","T3D_angle_mup","deltaPt","Transverse_mom_reco_mu","Transverse_mom_reco_pro","Number_protons","Muon_trackScore","Proton_trackScore","Flashmatch_deltaZ","Nu_score"}, data, {Reco_energy,Reco_class,Muon_length,Proton_length_leading,Proton_kinetic_leading,Proton_chi2mu,Proton_chi2pro,Muon_chi2mu,Muon_chi2pro,Vertex_x,Vertex_y,Vertex_z,Muon_endx,Muon_endy,Muon_endz,Proton_endx,Proton_endy,Proton_endz,Transverse_angle,T3D_angle_mup,Transverse_mom_reco,Transverse_mom_reco_mu,Transverse_mom_reco_pro,Number_protons,Muon_trackScore,Proton_trackScore,Flashmatch_deltaZ,Nu_score}, kSpillSelection, true);



//const Binning kBinPot     = Binning::Simple(1000,0,1.);

//Spectrum spot2("Pot",kBinPot,data,pot_original, kSpillSelection);


//dirx,y,z e muon length, proton leading length, number proton, vertex x,y,z 

  data.Go();

//cout << " POT " << final_pot << " Livetime " << final_livetime << endl;




  //TFile fout("/exp/sbnd/data/users/marterop/Ar23+_iterE/test_nolight/selected_sbnd_nolight_eff_pur_raster_vars_offbeam.root", "RECREATE");
  //TFile fout("/exp/sbnd/data/users/marterop/Ar23+_iterE/dent_nolight/selected_sbnd_eff_J_truep_v2.root", "RECREATE");
  TFile fout("/exp/sbnd/data/users/marterop/Ar23+_iterE/dent_nolight/selected_sbnd_eff_1e20_8_truep.root", "RECREATE");

  
  //TFile fout("selected_events_sBruce_10DATA_v2.root", "RECREATE");

  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  nutree_reco.SaveTo(dir); 
  //nutree_reco2.SaveTo(dir); 
  //nutree_reco3.SaveTo(dir); 


}

