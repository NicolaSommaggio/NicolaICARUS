////////////////////////////////////////////////////////////////////
// Definition of the spectra used for the stitching analysis
//
// Author: Alessandro Maria Ricci
////////////////////////////////////////////////////////////////////

#include "Stitch.h"

namespace ana {

    // return the end positions on the x-axis of nuMuCC tracks before the stitching
    // NOTE: for the split tracks take the segment closest to the vertex
    const SpillMultiVar kXendB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(sr->hdr.ismc) {
            // MC: loop on the slices
            // select mu primary daugthers of numuCC interactions, contained
            // in a single cryostat and with a minimal length of the track
            std::vector<int> nuMuCC;
            for(size_t i=0; i<sr->mc.nnu; i++) {
                for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                    if(sr->mc.nu.at(i).prim.at(j).pdg != std::abs(13)) continue;
                    if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                    if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                    nuMuCC.push_back(sr->mc.nu.at(i).prim.at(j).G4ID);
                }
            }
            if(nuMuCC.empty()) return Trk;
            for(size_t i=0; i<nuMuCC.size(); i++) {
                // take all segments associated to the same muon
                std::vector<PFP> pfp;
                for(int j=0; j<sr->nslc; j++) {
                    if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                    for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                        if(nuMuCC.at(i)!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                        if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                        PFP P;
                        P.vertex.x = sr->slc.at(j).vertex.x;
                        P.vertex.y = sr->slc.at(j).vertex.y;
                        P.vertex.z = sr->slc.at(j).vertex.z;
                        P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                        P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                        P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                        P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                        P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                        P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                        P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                        pfp.push_back(P);
                    }
                }
                if(pfp.empty() || pfp.size()>2) continue;
                if(pfp.size() == 1) {
                    Trk.push_back(pfp.at(0).end.x);
                    continue;
                }
                // sort the segments along the z-axis according to the distance from the vertex
                for(size_t j=0; j<pfp.size(); j++) {
                    caf::SRVector3D vtx;
                    vtx.x = pfp.at(j).vertex.x;
                    vtx.y = pfp.at(j).vertex.y;
                    vtx.z = pfp.at(j).vertex.z;
                    std::vector<float> ordz;
                    ordz.push_back(vtx.z);
                    for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                    std::sort(ordz.begin(), ordz.end());
                    // if the vertex is in the middle try again with the vertices associated to other segments
                    if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                    // sort start and end points of the segments according to the distance from the vertex
                    for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                    // fill the histogram
                    for(size_t k=0; k<pfp.size(); k++) {
                        if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                            (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                                Trk.push_back(pfp.at(k).end.x);
                                break;
                        }
                    }
                    break;
                }
            }
        }
        else {
            // DATA: loop on the slices
            for(int i=0; i<sr->nslc; i++) {
                // apply the 1muNp selection
                if(!kIcarus202507Contained1muNp(sr->slc.at(i))) continue;
                int idx = kIcarus202507MuonIdx(sr->slc.at(i));
                Trk.push_back(sr->slc.at(i).reco.pfp.at(idx).trk.end.x);
            }
        }
        return Trk;
    });


    // return the end positions on the x-axis of nuMuCC tracks after the stitching
    // NOTE: for the non-stitched tracks take the segment closest to the vertex
    const SpillMultiVar kXendA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(sr->hdr.ismc) {
            // MC: loop on the slices
            if(MCmuons.empty()) return Trk;
            for(size_t i=0; i<MCmuons.size(); i++) {
                // take all segments associated to the same muon
                std::vector<PFP> pfp;
                for(int j=0; j<sr->nslc; j++) {
                    if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                    for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                        if(MCmuons.at(i).G4ID!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                        if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                        PFP P;
                        P.vertex.x = sr->slc.at(j).vertex.x;
                        P.vertex.y = sr->slc.at(j).vertex.y;
                        P.vertex.z = sr->slc.at(j).vertex.z;
                        P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                        P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                        P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                        P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                        P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                        P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                        P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                        pfp.push_back(P);
                    }
                }
                if(pfp.empty() || pfp.size()>2) continue;
                // check if the muon has been stitched
                bool stitch = false;
                for(size_t j=0; j<sth.size(); j++) {
                    if(MCmuons.at(i).G4ID==sth.at(j).G4ID && sth.at(j).nuMuCC==true) {
                        Trk.push_back(sth.at(j).end2.x);
                        stitch = true;
                        break;
                    }
                }
                if(stitch) continue;
                if(pfp.size() == 1) {
                    /*if(std::abs(pfp.at(0).end.x) > 208 && std::abs(pfp.at(0).end.x) < 212) {
                        std::cout << " " << std::endl;
                        std::cout << "pfp reco: start (" << pfp.at(0).start.x << ", " << pfp.at(0).start.y << ", " << pfp.at(0).start.z <<
                                ") - end (" << pfp.at(0).end.x << ", " << pfp.at(0).end.y << ", " << pfp.at(0).end.z << ")" << std::endl;
                        for(int j=0; j<sr->nslc; j++) {
                            if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                            for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                                if(MCmuons.at(i).G4ID!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                                if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                                PFP P;
                                P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.truth.p.start.x;
                                P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.truth.p.start.y;
                                P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.truth.p.start.z;
                                P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.truth.p.end.x;
                                P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.truth.p.end.y;
                                P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.truth.p.end.z;
                                std::cout << "pfp truth: start (" << P.start.x << ", " << P.start.y << ", " << P.start.z <<
                                                    ") - end (" << P.end.x << ", " << P.end.y << ", " << P.end.z << ")" << std::endl;
                            }
                        }
                    }*/
                    Trk.push_back(pfp.at(0).end.x);
                    continue;
                }
                // sort the segments along the z-axis according to the distance from the vertex
                for(size_t j=0; j<pfp.size(); j++) {
                    caf::SRVector3D vtx;
                    vtx.x = pfp.at(j).vertex.x;
                    vtx.y = pfp.at(j).vertex.y;
                    vtx.z = pfp.at(j).vertex.z;
                    std::vector<float> ordz;
                    ordz.push_back(vtx.z);
                    for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                    std::sort(ordz.begin(), ordz.end());
                    // if the vertex is in the middle try again with the vertices associated to other segments
                    if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                    // sort start and end points of the segments according to the distance from the vertex
                    for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                    // fill the histogram
                    for(size_t k=0; k<pfp.size(); k++) {
                        if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                            (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                                Trk.push_back(pfp.at(k).end.x);
                                break;
                        }
                    }
                    break;
                }
            }
        }
        else {
            // DATA: loop on the slices
            for(int i=0; i<sr->nslc; i++) {
                // apply the 1muNp selection
                if(!kIcarus202507Contained1muNp(sr->slc.at(i))) continue;
                int idx = kIcarus202507MuonIdx(sr->slc.at(i));
                // check if the muon has been stitched
                bool push = true;
                for(size_t k=0; k<sth.size(); k++) {
                    if(sr->slc.at(i).reco.pfp.at(idx).id==sth.at(k).id || sr->slc.at(i).reco.pfp.at(idx).id==sth.at(k).id2) {
                        Trk.push_back(sth.at(k).end2.x);
                        push = false;
                        break;
                    }
                }
                if(push) Trk.push_back(sr->slc.at(i).reco.pfp.at(idx).trk.end.x);
            }
        }
        return Trk;
    });


    // return the end positions on the z-axis of nuMuCC tracks before the stitching
    // NOTE: for the split tracks take the segment closest to the vertex
    const SpillMultiVar kZendB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(sr->hdr.ismc) {
            // MC: loop on the slices
            // select mu primary daugthers of numuCC interactions, contained
            // in a single cryostat and with a minimal length of the track
            std::vector<int> nuMuCC;
            for(size_t i=0; i<sr->mc.nnu; i++) {
                for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                    if(sr->mc.nu.at(i).prim.at(j).pdg != std::abs(13)) continue;
                    if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                    if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                    nuMuCC.push_back(sr->mc.nu.at(i).prim.at(j).G4ID);
                }
            }
            if(nuMuCC.empty()) return Trk;
            for(size_t i=0; i<nuMuCC.size(); i++) {
                // take all segments associated to the same muon
                std::vector<PFP> pfp;
                for(int j=0; j<sr->nslc; j++) {
                    if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                    for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                        if(nuMuCC.at(i)!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                        if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                        PFP P;
                        P.vertex.x = sr->slc.at(j).vertex.x;
                        P.vertex.y = sr->slc.at(j).vertex.y;
                        P.vertex.z = sr->slc.at(j).vertex.z;
                        P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                        P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                        P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                        P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                        P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                        P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                        P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                        pfp.push_back(P);
                    }
                }
                if(pfp.empty() || pfp.size()>2) continue;
                if(pfp.size() == 1) {
                    Trk.push_back(pfp.at(0).end.z);
                    continue;
                }
                // sort the segments along the z-axis according to the distance from the vertex
                for(size_t j=0; j<pfp.size(); j++) {
                    caf::SRVector3D vtx;
                    vtx.x = pfp.at(j).vertex.x;
                    vtx.y = pfp.at(j).vertex.y;
                    vtx.z = pfp.at(j).vertex.z;
                    std::vector<float> ordz;
                    ordz.push_back(vtx.z);
                    for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                    std::sort(ordz.begin(), ordz.end());
                    // if the vertex is in the middle try again with the vertices associated to other segments
                    if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                    // sort start and end points of the segments according to the distance from the vertex
                    for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                    // fill the histogram
                    for(size_t k=0; k<pfp.size(); k++) {
                        if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                            (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                                Trk.push_back(pfp.at(k).end.z);
                                break;
                        }
                    }
                    break;
                }
            }
        }
        else {
            // DATA: loop on the slices
            for(int i=0; i<sr->nslc; i++) {
                // apply the 1muNp selection
                if(!kIcarus202507Contained1muNp(sr->slc.at(i))) continue;
                int idx = kIcarus202507MuonIdx(sr->slc.at(i));
                Trk.push_back(sr->slc.at(i).reco.pfp.at(idx).trk.end.z);
            }
        }
        return Trk;
    });


    // return the end positions on the z-axis of nuMuCC tracks after the stitching
    // NOTE: for the non-stitched tracks take the segment closest to the vertex
    const SpillMultiVar kZendA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(sr->hdr.ismc) {
            if(MCmuons.empty()) return Trk;
            for(size_t i=0; i<MCmuons.size(); i++) {
                // take all segments associated to the same muon
                std::vector<PFP> pfp;
                for(int j=0; j<sr->nslc; j++) {
                    if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                    for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                        if(MCmuons.at(i).G4ID!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                        if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                        PFP P;
                        P.vertex.x = sr->slc.at(j).vertex.x;
                        P.vertex.y = sr->slc.at(j).vertex.y;
                        P.vertex.z = sr->slc.at(j).vertex.z;
                        P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                        P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                        P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                        P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                        P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                        P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                        P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                        pfp.push_back(P);
                    }
                }
                if(pfp.empty() || pfp.size()>2) continue;
                // check if the muon has been stitched
                bool stitch = false;
                for(size_t j=0; j<sth.size(); j++) {
                    if(MCmuons.at(i).G4ID==sth.at(j).G4ID && sth.at(j).nuMuCC==true) {
                        Trk.push_back(sth.at(j).end2.z);
                        stitch = true;
                        break;
                    }
                }
                if(stitch) continue;
                if(pfp.size() == 1) {
                    Trk.push_back(pfp.at(0).end.z);
                    continue;
                }
                // sort the segments along the z-axis according to the distance from the vertex
                for(size_t j=0; j<pfp.size(); j++) {
                    caf::SRVector3D vtx;
                    vtx.x = pfp.at(j).vertex.x;
                    vtx.y = pfp.at(j).vertex.y;
                    vtx.z = pfp.at(j).vertex.z;
                    std::vector<float> ordz;
                    ordz.push_back(vtx.z);
                    for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                    std::sort(ordz.begin(), ordz.end());
                    // if the vertex is in the middle try again with the vertices associated to other segments
                    if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                    // sort start and end points of the segments according to the distance from the vertex
                    for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                    // fill the histogram
                    for(size_t k=0; k<pfp.size(); k++) {
                        if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                            (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                                Trk.push_back(pfp.at(k).end.z);
                                break;
                        }
                    }
                    break;
                }
            }
        }
        else {
            // DATA: loop on the slices
            for(int i=0; i<sr->nslc; i++) {
                // apply the 1muNp selection
                if(!kIcarus202507Contained1muNp(sr->slc.at(i))) continue;
                int idx = kIcarus202507MuonIdx(sr->slc.at(i));
                // check if the muon has been stitched
                bool push = true;
                for(size_t k=0; k<sth.size(); k++) {
                    if(sr->slc.at(i).reco.pfp.at(idx).id==sth.at(k).id || sr->slc.at(i).reco.pfp.at(idx).id==sth.at(k).id2) {
                        Trk.push_back(sth.at(k).end2.z);
                        push = false;
                        break;
                    }
                }
                if(push) Trk.push_back(sr->slc.at(i).reco.pfp.at(idx).trk.end.z);
            }
        }
        return Trk;
    });


    // return the track length before the stitching
    const SpillMultiVar kSthLengthB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> length;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(sr->hdr.ismc) {
            // MC
            for(size_t i=0; i<sth.size(); i++) {
                if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
                length.push_back(sth.at(i).len);
            }
        }
        else {
            // DATA
            for(size_t i=0; i<sth.size(); i++) length.push_back(sth.at(i).len);
        }
        return length;
    });


    // return the track length after the stitching
    const SpillMultiVar kSthLengthA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> length;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(sr->hdr.ismc) {
            // MC
            for(size_t i=0; i<sth.size(); i++) {
                if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
                length.push_back(sth.at(i).Len);
            }
        }
        else {
            // DATA
            for(size_t i=0; i<sth.size(); i++) length.push_back(sth.at(i).Len);
        }
        return length;
    });


    // return the difference of the muon track length after - before the stitching
    const SpillMultiVar kSthLengthAminB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> length;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(sr->hdr.ismc) {
            // MC
            for(size_t i=0; i<sth.size(); i++) {
                if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
                length.push_back(sth.at(i).Len - sth.at(i).len);
            }
        }
        else {
            // DATA
            for(size_t i=0; i<sth.size(); i++) length.push_back(sth.at(i).Len - sth.at(i).len);
        }
        return length;
    });


    // return the track energy before the stitching
    const SpillMultiVar kSthEnergyB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> energy;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(sr->hdr.ismc) {
            // MC
            for(size_t i=0; i<sth.size(); i++) {
                if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
                energy.push_back( std::sqrt(sth.at(i).p_muon*sth.at(i).p_muon + mmu*mmu)*1e3 );
            }
        }
        else {
            // DATA
            for(size_t i=0; i<sth.size(); i++) energy.push_back( std::sqrt(sth.at(i).p_muon*sth.at(i).p_muon + mmu*mmu)*1e3 );
        }
        return energy;
    });


    // return the track energy after the stitching
    const SpillMultiVar kSthEnergyA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> energy;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(sr->hdr.ismc) {
            // MC
            for(size_t i=0; i<sth.size(); i++) {
                if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
                energy.push_back( std::sqrt( sth.at(i).P_muon*sth.at(i).P_muon + mmu*mmu)*1e3 );
            }
        }
        else {
            // DATA
            for(size_t i=0; i<sth.size(); i++) energy.push_back( std::sqrt( sth.at(i).P_muon*sth.at(i).P_muon + mmu*mmu)*1e3 );
        }
        return energy;
    });


    // MC only: return the end positions on the x-axis of split nuMuCC tracks before the stitching
    // NOTE: for the split tracks take the segment closest to the vertex
    const SpillMultiVar kSpTrkXendB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(!sr->hdr.ismc) return Trk;
        // select mu primary daugthers of numuCC interactions, contained
        // in a single cryostat and with a minimal length of the track
        std::vector<int> nuMuCC;
        for(size_t i=0; i<sr->mc.nnu; i++) {
            for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                if(sr->mc.nu.at(i).prim.at(j).pdg != std::abs(13)) continue;
                if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                nuMuCC.push_back(sr->mc.nu.at(i).prim.at(j).G4ID);
            }
        }
        if(nuMuCC.empty()) return Trk;
        for(size_t i=0; i<nuMuCC.size(); i++) {
            // take all segments associated to the same muon and select the split tracks
            std::vector<PFP> pfp;
            for(int j=0; j<sr->nslc; j++) {
                if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                    if(nuMuCC.at(i)!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                    if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                    PFP P;
                    P.vertex.x = sr->slc.at(j).vertex.x;
                    P.vertex.y = sr->slc.at(j).vertex.y;
                    P.vertex.z = sr->slc.at(j).vertex.z;
                    P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                    P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                    P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                    P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                    P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                    P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                    P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                    pfp.push_back(P);
                }
            }
            if(pfp.size()!=2) continue; // consider only tracks broken into 2 segments
            // sort the segments along the z-axis according to the distance from the vertex
            for(size_t j=0; j<pfp.size(); j++) {
                caf::SRVector3D vtx;
                vtx.x = pfp.at(j).vertex.x;
                vtx.y = pfp.at(j).vertex.y;
                vtx.z = pfp.at(j).vertex.z;
                std::vector<float> ordz;
                ordz.push_back(vtx.z);
                for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                std::sort(ordz.begin(), ordz.end());
                // if the vertex is in the middle try again with the vertices associated to other segments
                if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                // sort start and end points of the segments according to the distance from the vertex
                for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                // fill the histogram
                for(size_t k=0; k<pfp.size(); k++) {
                    if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                        (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                            Trk.push_back(pfp.at(k).end.x);
                            break;
                    }
                }
                break;
            }
        }
        return Trk;
    });


    // MC only: return the end positions on the x-axis of split nuMuCC tracks after the stitching
    // NOTE: for the non-stitched tracks take the segment closest to the vertex
    const SpillMultiVar kSpTrkXendA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(!sr->hdr.ismc) return Trk;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(MCmuons.empty()) return Trk;
        for(size_t i=0; i<MCmuons.size(); i++) {
            // take all segments associated to the same muon
            std::vector<PFP> pfp;
            for(int j=0; j<sr->nslc; j++) {
                if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                    if(MCmuons.at(i).G4ID!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                    if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                    PFP P;
                    P.vertex.x = sr->slc.at(j).vertex.x;
                    P.vertex.y = sr->slc.at(j).vertex.y;
                    P.vertex.z = sr->slc.at(j).vertex.z;
                    P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                    P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                    P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                    P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                    P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                    P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                    P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                    pfp.push_back(P);
                }
            }
            if(pfp.size()!=2) continue; // consider only tracks broken into 2 segments
            // check if the muon has been stitched
            bool stitch = false;
            for(size_t j=0; j<sth.size(); j++) {
                if(MCmuons.at(i).G4ID==sth.at(j).G4ID && sth.at(j).nuMuCC==true) {
                    Trk.push_back(sth.at(j).end2.x);
                    stitch = true;
                    break;
                }
            }
            if(stitch) continue;
            // sort the segments along the z-axis according to the distance from the vertex
            for(size_t j=0; j<pfp.size(); j++) {
                caf::SRVector3D vtx;
                vtx.x = pfp.at(j).vertex.x;
                vtx.y = pfp.at(j).vertex.y;
                vtx.z = pfp.at(j).vertex.z;
                std::vector<float> ordz;
                ordz.push_back(vtx.z);
                for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                std::sort(ordz.begin(), ordz.end());
                // if the vertex is in the middle try again with the vertices associated to other segments
                if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                // sort start and end points of the segments according to the distance from the vertex
                for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                // fill the histogram
                for(size_t k=0; k<pfp.size(); k++) {
                    if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                        (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                            Trk.push_back(pfp.at(k).end.x);
                            break;
                    }
                }
                break;
            }
        }
        return Trk;
    });


    // MC only: return the end positions on the z-axis of split nuMuCC tracks before the stitching
    // NOTE: for the split tracks take the segment closest to the vertex
    const SpillMultiVar kSpTrkZendB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(!sr->hdr.ismc) return Trk;
        // select mu primary daugthers of numuCC interactions, contained
        // in a single cryostat and with a minimal length of the track
        std::vector<int> nuMuCC;
        for(size_t i=0; i<sr->mc.nnu; i++) {
            for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                if(sr->mc.nu.at(i).prim.at(j).pdg != std::abs(13)) continue;
                if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                nuMuCC.push_back(sr->mc.nu.at(i).prim.at(j).G4ID);
            }
        }
        if(nuMuCC.empty()) return Trk;
        for(size_t i=0; i<nuMuCC.size(); i++) {
            // take all segments associated to the same muon and select the split tracks
            std::vector<PFP> pfp;
            for(int j=0; j<sr->nslc; j++) {
                if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                    if(nuMuCC.at(i)!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                    if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                    PFP P;
                    P.vertex.x = sr->slc.at(j).vertex.x;
                    P.vertex.y = sr->slc.at(j).vertex.y;
                    P.vertex.z = sr->slc.at(j).vertex.z;
                    P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                    P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                    P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                    P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                    P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                    P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                    P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                    pfp.push_back(P);
                }
            }
            if(pfp.size()!=2) continue; // consider only tracks broken into 2 segments
            // sort the segments along the z-axis according to the distance from the vertex
            for(size_t j=0; j<pfp.size(); j++) {
                caf::SRVector3D vtx;
                vtx.x = pfp.at(j).vertex.x;
                vtx.y = pfp.at(j).vertex.y;
                vtx.z = pfp.at(j).vertex.z;
                std::vector<float> ordz;
                ordz.push_back(vtx.z);
                for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                std::sort(ordz.begin(), ordz.end());
                // if the vertex is in the middle try again with the vertices associated to other segments
                if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                // sort start and end points of the segments according to the distance from the vertex
                for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                // fill the histogram
                for(size_t k=0; k<pfp.size(); k++) {
                    if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                        (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                            Trk.push_back(pfp.at(k).end.z);
                            break;
                    }
                }
                break;
            }
        }
        return Trk;
    });


    // MC only: return the end positions on the z-axis of split nuMuCC tracks after the stitching
    // NOTE: for the non-stitched tracks take the segment closest to the vertex
    const SpillMultiVar kSpTrkZendA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(!sr->hdr.ismc) return Trk;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        if(MCmuons.empty()) return Trk;
        for(size_t i=0; i<MCmuons.size(); i++) {
            // take all segments associated to the same muon
            std::vector<PFP> pfp;
            for(int j=0; j<sr->nslc; j++) {
                if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                    if(MCmuons.at(i).G4ID!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                    if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                    PFP P;
                    P.vertex.x = sr->slc.at(j).vertex.x;
                    P.vertex.y = sr->slc.at(j).vertex.y;
                    P.vertex.z = sr->slc.at(j).vertex.z;
                    P.start.x = sr->slc.at(j).reco.pfp.at(k).trk.start.x;
                    P.start.y = sr->slc.at(j).reco.pfp.at(k).trk.start.y;
                    P.start.z = sr->slc.at(j).reco.pfp.at(k).trk.start.z;
                    P.end.x = sr->slc.at(j).reco.pfp.at(k).trk.end.x;
                    P.end.y = sr->slc.at(j).reco.pfp.at(k).trk.end.y;
                    P.end.z = sr->slc.at(j).reco.pfp.at(k).trk.end.z;
                    P.slcID = sr->slc.at(j).reco.pfp.at(k).slcID;
                    pfp.push_back(P);
                }
            }
            if(pfp.size()!=2) continue; // consider only tracks broken into 2 segments
            // check if the muon has been stitched
            bool stitch = false;
            for(size_t j=0; j<sth.size(); j++) {
                if(MCmuons.at(i).G4ID==sth.at(j).G4ID && sth.at(j).nuMuCC==true) {
                    Trk.push_back(sth.at(j).end2.z);
                    stitch = true;
                    break;
                }
            }
            if(stitch) continue;
            // sort the segments along the z-axis according to the distance from the vertex
            for(size_t j=0; j<pfp.size(); j++) {
                caf::SRVector3D vtx;
                vtx.x = pfp.at(j).vertex.x;
                vtx.y = pfp.at(j).vertex.y;
                vtx.z = pfp.at(j).vertex.z;
                std::vector<float> ordz;
                ordz.push_back(vtx.z);
                for(size_t k=0; k<pfp.size(); k++) ordz.push_back( 0.5*(pfp.at(k).start.z+pfp.at(k).end.z) );
                std::sort(ordz.begin(), ordz.end());
                // if the vertex is in the middle try again with the vertices associated to other segments
                if(ordz.at(0)!=vtx.z && ordz.back()!=vtx.z) continue;
                // sort start and end points of the segments according to the distance from the vertex
                for(size_t k=0; k<pfp.size(); k++) Icarus202503SortPFP(pfp.at(k), vtx);
                // fill the histogram
                for(size_t k=0; k<pfp.size(); k++) {
                    if( (ordz.at(0)==vtx.z && ordz.at(1)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ||
                        (ordz.back()==vtx.z && ordz.at(ordz.size()-2)==0.5*(pfp.at(k).start.z+pfp.at(k).end.z)) ) {
                            Trk.push_back(pfp.at(k).end.z);
                            break;
                    }
                }
                break;
            }
        }
        return Trk;
    });


    // MC only: return (Lreco-Ltrue)/Ltrue before the stitching
    const SpillMultiVar kSthLengthRatioB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> length;
        if(!sr->hdr.ismc) return length;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        for(size_t i=0; i<sth.size(); i++) {
            if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
            length.push_back( std::abs(sth.at(i).len-sth.at(i).MClen) / sth.at(i).MClen );
        }
        return length;
    });


    // MC only: return (Lreco-Ltrue)/Ltrue after the stitching
    const SpillMultiVar kSthLengthRatioA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> length;
        if(!sr->hdr.ismc) return length;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        for(size_t i=0; i<sth.size(); i++) {
            if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
            length.push_back( std::abs(sth.at(i).Len-sth.at(i).MClen) / sth.at(i).MClen );
        }
        return length;
    });


    // MC only: return (Ereco-Etrue)/Etrue before the stitching
    const SpillMultiVar kSthEnergyRatioB([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> energy;
        if(!sr->hdr.ismc) return energy;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        for(size_t i=0; i<sth.size(); i++) {
            if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
            float MCEn = std::sqrt(sth.at(i).MCp_muon*sth.at(i).MCp_muon + mmu*mmu)*1e3;
            float E = std::sqrt(sth.at(i).p_muon*sth.at(i).p_muon + mmu*mmu)*1e3;
            energy.push_back( std::abs(E-MCEn) / MCEn );
        }
        return energy;
    });


    // MC only: return (Ereco-Etrue)/Etrue after the stitching
    const SpillMultiVar kSthEnergyRatioA([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> energy;
        if(!sr->hdr.ismc) return energy;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        calibration calib;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        for(size_t i=0; i<sth.size(); i++) {
            if(sth.at(i).nuMuCC==false) continue;   // skip non-nuMuCC tracks
            float MCEn = std::sqrt(sth.at(i).MCp_muon*sth.at(i).MCp_muon + mmu*mmu)*1e3;
            float E = std::sqrt(sth.at(i).P_muon*sth.at(i).P_muon + mmu*mmu)*1e3;
            energy.push_back( std::abs(E-MCEn) / MCEn );
        }
        return energy;
    });


    // MC only: return the PDG of the split tracks
    const SpillMultiVar kSpTrkPDG([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(!sr->hdr.ismc) return Trk;
        std::vector<int> nuMuCC;
        // select mu primary daugthers of numuCC interactions, contained in a single cryostat
        // and with a minimal length of the track
        for(size_t i=0; i<sr->mc.nnu; i++) {
            for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                if(sr->mc.nu.at(i).prim.at(j).pdg != std::abs(13)) continue;
                if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                nuMuCC.push_back(sr->mc.nu.at(i).prim.at(j).G4ID);
            }
        }
        if(nuMuCC.empty()) return Trk;
        // select the split tracks
        for(size_t i=0; i<nuMuCC.size(); i++) {
            size_t counter = 0;
            for(int j=0; j<sr->nslc; j++) {
                if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                    if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                    if(nuMuCC.at(i)==sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) counter++;
                }
            }
            if(counter<2) continue; // skip non-split tracks
            // fill the histogram
            for(int j=0; j<sr->nslc; j++) {
                if(!kIcarus202412RecoFiducial(sr->slc.at(j))) continue;
                for(size_t k=0; k<sr->slc.at(j).reco.npfp; k++) {
                    if(nuMuCC.at(i)!=sr->slc.at(j).reco.pfp.at(k).trk.truth.p.G4ID) continue;
                    if(!kIcarus202401MuonTrack(sr->slc.at(j).reco.pfp.at(k), param::trkScore, param::trkLmin)) continue;
                    if(std::abs(sr->slc.at(j).reco.pfp.at(k).trk.truth.p.pdg)==11) Trk.push_back(0);         // e+,e-
                    else if(std::abs(sr->slc.at(j).reco.pfp.at(k).trk.truth.p.pdg)==13) Trk.push_back(1);    // mu+,mu-
                    else if(std::abs(sr->slc.at(j).reco.pfp.at(k).trk.truth.p.pdg)==211) Trk.push_back(2);   // pi+,pi-
                    else if(std::abs(sr->slc.at(j).reco.pfp.at(k).trk.truth.p.pdg)==2212) Trk.push_back(3);  // p
                    else if(std::abs(sr->slc.at(j).reco.pfp.at(k).trk.truth.p.pdg)==22) Trk.push_back(4);    // gamma
                    else Trk.push_back(5);                                                                   // other
                }
            }
        }
        return Trk;
    });


    // MC only: return the PDG of the stitched tracks
    const SpillMultiVar kSthTrkPDG([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> Trk;
        if(!sr->hdr.ismc) return Trk;
        std::vector<PFP> sth;
        std::vector<muon> MCmuons;
        std::vector<int> G4ID;
        calibration calib;
        bool found1 = false;
        bool found2 = false;
        sth = stitch(sr, MCmuons, calib, false, sr->hdr.ismc);
        for(size_t i=0; i<sth.size(); i++) {
            //if(sth.at(i).nuMuCC==false) continue; // skip non-nuMuCC interactions
            // skip G4IDs already considered
            for(size_t j=0; j<G4ID.size(); j++) {
                if(G4ID.at(j)==sth.at(i).G4ID) found1 = true;
                else if(G4ID.at(j)==sth.at(i).G4ID2) found2 = true;
            }
            if(!found1) {
                G4ID.push_back(sth.at(i).G4ID);
                if(std::abs(sth.at(i).pdg)==11) Trk.push_back(0);           // e+,e-
                else if(std::abs(sth.at(i).pdg)==13) Trk.push_back(1);      // mu+,mu-
                else if(std::abs(sth.at(i).pdg)==211) Trk.push_back(2);     // pi+,pi-
                else if(std::abs(sth.at(i).pdg)==2212) Trk.push_back(3);    // p
                else if(std::abs(sth.at(i).pdg)==22) Trk.push_back(4);      // gamma
                else Trk.push_back(5);                                      // other
            }
            if(!found2) {
                G4ID.push_back(sth.at(i).G4ID2);
                if(std::abs(sth.at(i).pdg2)==11) Trk.push_back(0);          // e+,e-
                else if(std::abs(sth.at(i).pdg2)==13) Trk.push_back(1);     // mu+,mu-
                else if(std::abs(sth.at(i).pdg2)==211) Trk.push_back(2);    // pi+,pi-
                else if(std::abs(sth.at(i).pdg2)==2212) Trk.push_back(3);   // p
                else if(std::abs(sth.at(i).pdg2)==22) Trk.push_back(4);     // gamma
                else Trk.push_back(5);                                      // other
            }
        }
        return Trk;
    });

} // namespace ana