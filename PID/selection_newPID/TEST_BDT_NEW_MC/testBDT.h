#include "../helper_generic.h"
#include "../../selection1muNp_PIDfinale/helper_newPID.h"



ofstream dump_test_bdt("DUMP_TEST_BDT_CUT_OVERLAYS.txt");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////                                                 testBDT                                                      //////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const SpillMultiVar testBDT([](const caf::SRSpillProxy* sr)-> std::vector<double>
{

    std::vector<double> vector_active;

    int nslice = -1;
    for(const auto &islc : sr->slc)
    {
        nslice = nslice + 1;
        
        TVector3 vertex_reco;
        TVector3 vertex_true;
        vertex_reco.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
        vertex_true.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    
        if (islc.tmatch.eff <= 0.5 || (vertex_true-vertex_reco).Mag() >= 100 || islc.is_clear_cosmic)continue;

        for (std::size_t ipfp(0); ipfp < islc.reco.npfp; ++ipfp)
        {
            int bestplane = find_best_plane(islc,ipfp);

            if(bestplane == -1){continue;}

            int true_class_pfp = true_selection(sr,islc,ipfp);
            

            double depE = compute_depE(islc,ipfp,bestplane);

            int prediction = PIDclass(islc,ipfp);
            std::vector<double> pred_porba = PIDproba(islc,ipfp);

            TVector3 Start_mom_v_proton; 
            Start_mom_v_proton.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
            double Ek_proton = sqrt(pow(938.3,2)+pow(Start_mom_v_proton.Mag()*1000,2))-938.3;

            TVector3 Start_mom_v_pion; 
            Start_mom_v_pion.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);
            double Ek_pion = sqrt(pow(139.570,2)+pow(Start_mom_v_pion.Mag()*1000,2))-139.570;

            double Ek = -1;
            if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==2212){Ek=Ek_proton;}
            else if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211){Ek=Ek_pion;}

            double start_distance = std::sqrt(std::pow(islc.reco.pfp[ipfp].trk.start.x-islc.vertex.x,2)+std::pow(islc.reco.pfp[ipfp].trk.start.y-islc.vertex.y,2)+std::pow(islc.reco.pfp[ipfp].trk.start.z-islc.vertex.z,2));
            double end_distance = std::sqrt(std::pow(islc.reco.pfp[ipfp].trk.end.x-islc.vertex.x,2)+std::pow(islc.reco.pfp[ipfp].trk.end.y-islc.vertex.y,2)+std::pow(islc.reco.pfp[ipfp].trk.end.z-islc.vertex.z,2));

            double _bar_flash = bar_flash_x(sr,islc);
            double _bar_flash_x = bar_flash_x(sr,islc);

            std::string true_class = classification_type_generic(sr,islc);


            dump_test_bdt << islc.reco.pfp[ipfp].trk.truth.p.pdg << " "; 
            dump_test_bdt << true_class_pfp << " ";
            dump_test_bdt << islc.reco.pfp[ipfp].trk.len << " ";
            dump_test_bdt << Ek << " ";
            dump_test_bdt << prediction << " ";
            for(const auto &p : pred_porba){dump_test_bdt << p << " ";}
            dump_test_bdt << islc.reco.pfp[ipfp].trackScore << " ";
            dump_test_bdt << islc.reco.pfp[ipfp].parent_is_primary << " ";
            dump_test_bdt << start_distance << " ";
            dump_test_bdt << end_distance << " ";
            dump_test_bdt << isInContained(islc.reco.pfp[ipfp].trk.end.x,islc.reco.pfp[ipfp].trk.end.y,islc.reco.pfp[ipfp].trk.end.z,10) << " ";
            dump_test_bdt << islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness << " ";
            dump_test_bdt << islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity << " ";
            dump_test_bdt << _bar_flash << " ";
            dump_test_bdt << _bar_flash_x << " ";
            dump_test_bdt << true_class << " ";
            dump_test_bdt << depE << endl;


        }//loop over all pfps   


    }//loop over all slices

  return vector_active;
});
