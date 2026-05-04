////////////////////////////////////////////////////////////////////
// Definition of the 1muNp selection at spill level
//
// Author: Alessandro Maria Ricci
////////////////////////////////////////////////////////////////////

#ifndef SELECTION_H
#define SELECTION_H

#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202412.h"

namespace ana {

    static constexpr float mp  = 0.9383f;   // proton mass in GeV/c^2
    static constexpr float mpi = 0.13957f;  // pion mass in GeV/c^2
    static constexpr float mmu = 0.10566f;  // muon mass in GeV/c^2

    // Mostly just copied from NumuVarsIcarus202106.cxx
    int kIcarus202507MuonIdx(const caf::SRSliceProxy &slc) {
        // The (dis)qualification of a slice is based upon the track level information.
        float Longest(0);
        int PTrackInd(-1);
        for (std::size_t i(0); i < slc.reco.npfp; ++i) {
            auto const& pfp = slc.reco.pfp.at(i);
            if (pfp.trackScore < 0.5) { continue; }
	        //if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x)){ continue; }
            auto const& trk = pfp.trk;
            if(std::isnan(slc.vertex.x) || std::isnan(slc.vertex.y) || std::isnan(slc.vertex.z)) continue;
            if(std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.start.y) || std::isnan(pfp.trk.start.z)) continue;
            const float Atslc = std::hypot(slc.vertex.x - trk.start.x,
                                            slc.vertex.y - trk.start.y,
                                            slc.vertex.z - trk.start.z);
            const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

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
            const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
            const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
            if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest ) {
                Longest = trk.len;
                PTrackInd = i;
            }
        }
        return PTrackInd;
    }

    static bool Icarus202401_proton_cut(const caf::SRTrackProxy& trk) {
        //int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
        int plane = 2; // Hard code collection plane for now since induction 2 has peak at higher chi2
        auto chi2 = chi2pid::chi2_calculator.calculate_chi2(trk.calo[plane]);
        return chi2.chi2_proton < 100 && chi2.chi2_proton > 0;
    }

    static bool Icarus202401_proton_cut_exist(const caf::SRTrackProxy& trk) { 
        //int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
        int plane = 2; // Hard code collection plane for now since induction 2 has peak at higher chi2
        auto chi2 = chi2pid::chi2_calculator.calculate_chi2(trk.calo[plane]);
        return chi2.chi2_proton > 0;
    }

    int kIcarus202507NumProtons(const caf::SRSliceProxy &slc) {
        int count = 0;
        auto idx = kIcarus202507MuonIdx(slc);
        int muID = -1;
        if (idx >= 0) muID = slc.reco.pfp.at(idx).id;
        for(auto& pfp: slc.reco.pfp) {
            if (pfp.trackScore < 0.4) { continue; }
            //if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
            auto const& trk = pfp.trk;
            if(std::isnan(slc.vertex.x) || std::isnan(slc.vertex.y) || std::isnan(slc.vertex.z)) continue;
            if(std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.start.y) || std::isnan(pfp.trk.start.z)) continue;
            const float Atslc = std::hypot(slc.vertex.x - trk.start.x,
                                            slc.vertex.y - trk.start.y,
                                            slc.vertex.z - trk.start.z);
            const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
            if(pfp.id != muID && Icarus202401_proton_cut(trk) && AtSlice && std::hypot(pfp.trk.rangeP.p_proton, mp) - mp > 0.05)
            ++count;
        }
        return count;
    }

    int kIcarus202507NumShowers(const caf::SRSliceProxy &slc) {
        int count = 0;
        for(const auto& pfp: slc.reco.pfp) {
            if(std::isnan(slc.vertex.x) || std::isnan(slc.vertex.y) || std::isnan(slc.vertex.z)) continue;
            if(std::isnan(pfp.shw.start.x) || std::isnan(pfp.shw.start.y) || std::isnan(pfp.shw.start.z)) continue;
            const float Atslc = std::hypot(slc.vertex.x - pfp.shw.start.x,
                                            slc.vertex.y - pfp.shw.start.y,
                                            slc.vertex.z - pfp.shw.start.z);
            const bool AtSlice = ( Atslc < 50.0 && pfp.parent_is_primary);
            const bool isShower = pfp.trackScore > 0 && (pfp.trackScore < 0.4 || (pfp.trackScore < 0.5 && !Icarus202401_proton_cut(pfp.trk)));
            if(isShower && AtSlice && pfp.shw.plane[2].energy > 0.025) count++;
        }
        return count;
    }

    int kIcarus202507NumPions(const caf::SRSliceProxy &slc) {
        int count = 0;
        auto idx = kIcarus202507MuonIdx(slc);
        int muID = -1;
        if (idx >= 0) muID = slc.reco.pfp.at(idx).id;
        for(auto& pfp: slc.reco.pfp) {
            if (pfp.trackScore < 0.5) { continue; }
            auto const& trk = pfp.trk;
            if(std::isnan(slc.vertex.x) || std::isnan(slc.vertex.y) || std::isnan(slc.vertex.z)) continue;
            if(std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.start.y) || std::isnan(pfp.trk.start.z)) continue;
            const float Atslc = std::hypot(slc.vertex.x - trk.start.x,
                                            slc.vertex.y - trk.start.y,
                                            slc.vertex.z - trk.start.z);
            const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
            //if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
            if(pfp.id != muID && !Icarus202401_proton_cut(trk) && Icarus202401_proton_cut_exist(trk) && AtSlice && std::hypot(pfp.trk.rangeP.p_pion, mpi) - mpi > 0.025)
            ++count;
        }
        return count;
    }

    static bool Icarus202401contained(const caf::SRTrackProxy& trk) {
        return ( !isnan(trk.end.x) &&
                (( trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5 ) ||
                 ( trk.end.x > 61.94 + 5 && trk.end.x < +358.49 - 5 )) &&
                !isnan(trk.end.y) &&
                ( trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5 ) &&
                !isnan(trk.end.z) &&
                ( trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5 ) );
    }

    bool kIcarus202507ContainedHadrons(const caf::SRSliceProxy &slc) {
        auto idx = kIcarus202507MuonIdx(slc);
        int muID = -1;
        if (idx >= 0) muID = slc.reco.pfp.at(idx).id;
        for(auto& pfp: slc.reco.pfp) {
            //if (pfp.trackScore < 0.5) { continue; }

            if(std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x) || std::isnan(pfp.trk.len)) continue;
            auto const& trk = pfp.trk;
            //if(pfp.id != muID && pfp.parent_is_primary)
            if(pfp.id != muID)
            if(!Icarus202401contained(trk)) return false;
        }
        return true;
    }

    bool kIcarus202507FoundMuon(const caf::SRSliceProxy &slc) {
        return kIcarus202507MuonIdx(slc) >= 0;
    }

    bool kIcarus202507ContainedMuon(const caf::SRSliceProxy &slc) {
        return kIcarus202507FoundMuon(slc) && Icarus202401contained(slc.reco.pfp.at(kIcarus202507MuonIdx(slc)).trk);
    }

    bool kIcarus202507NumuSelection(const caf::SRSliceProxy &slc) {
        return kIcarus202412RecoFiducial(slc) && kIcarus202401BaryFMCut(slc) && kIcarus202507FoundMuon(slc);
    }

    bool kIcarus202507NoPion(const caf::SRSliceProxy &slc) {
        return kIcarus202507NumPions(slc) == 0;
    }

    bool kIcarus202507ContainedNumuCCInclusive(const caf::SRSliceProxy &slc) {
        return kIcarus202507NumuSelection(slc) && kIcarus202507ContainedHadrons(slc) && kIcarus202507ContainedMuon(slc);
    }

    bool kIcarus202507NumuCC0pi(const caf::SRSliceProxy &slc) {
        return kIcarus202507ContainedNumuCCInclusive(slc) && kIcarus202507NoPion(slc);
    }

    bool kIcarus202507Contained1muNp(const caf::SRSliceProxy &slc) {
        return kIcarus202507NumuCC0pi(slc) && kIcarus202507NumProtons(slc) > 0 && kIcarus202507NumShowers(slc) == 0;
    }

}

#endif // SELECTION_H