////////////////////////////////////////////////////////////////////
// Definition of the stitching function between slices
//
// Author: Alessandro Maria Ricci
////////////////////////////////////////////////////////////////////

#ifndef STITCH2_H
#define STITCH2_H

#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202412.h"

namespace ana {

    // stitching function in the event
    std::vector<double> stitch2(const caf::SRSpillProxy* sr, std::vector<double> &hist, std::vector<int> &sthG4ID, bool mc) {
        std::vector<double> stitch;
        std::vector<muon> MCmuons;
        std::vector<PFP> pfp;
        size_t nStitch = 0;
        if(mc) {
            // select mu primary daugthers of numuCC interactions, contained in a single cryostat
            // and with a minimal length of the track
            for(size_t i=0; i<sr->mc.nnu; i++) {
                for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                    if(sr->mc.nu.at(i).prim.at(j).pdg != 13 && sr->mc.nu.at(i).prim.at(j).pdg != -13) continue;
                    if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                    if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                    muon mu = {sr->mc.nu.at(i).prim.at(j).G4ID, 0, false, false, false};
                    MCmuons.push_back(mu);
                }
            }
            if(MCmuons.empty()) return stitch;
        }
        // loop on the slices to import the pfp
        for(int i=0; i<sr->nslc; i++) {
            // check if the slice is in the fiducial volume
            if(!kIcarus202412RecoFiducial(sr->slc.at(i))) continue;
            for(size_t j=0; j<sr->slc.at(i).reco.npfp; j++) {
                // check if the pfp is a muon
                if(!kIcarus202401MuonTrack(sr->slc.at(i).reco.pfp.at(j), param::trkScore, param::trkLmin)) continue;
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
                P.slcID = sr->slc.at(i).reco.pfp.at(j).slcID;
                P.id = sr->slc.at(i).reco.pfp.at(j).id;
                // MC
                if(mc) {
                    P.G4ID = sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.G4ID;
                    P.pdg = sr->slc.at(i).reco.pfp.at(j).trk.truth.p.pdg;
                    //P.energy_comp = sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.energy_completeness;
                    //P.hit_comp = sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.hit_completeness;
                    // quality cut on hit_comp and E_comp
                    /*if(param::EHcomp) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(MCmuons.at(k).G4ID == P.G4ID && P.energy_comp >= 0.1 && P.hit_comp >= 0.1) {
                                MCmuons.at(k).counter++;
                                if(MCmuons.at(k).counter >= 2) {
                                    MCmuons.at(k).split = true;
                                }
                            }
                        }
                    }*/
                }
                pfp.push_back(P);
            }
        }
        // quality cut on hit_comp and E_comp
        /*if(mc && param::EHcomp) {
            // delete all elements that satisfy the condition split = false
            MCmuons.erase(std::remove_if(MCmuons.begin(), MCmuons.end(), [](const muon& mu) {
                return mu.split == false;
            }), MCmuons.end());
        }*/
        // order the pfp according to the length
        std::sort(pfp.begin(), pfp.end(), [](const PFP a, const PFP b) {return a.len > b.len;});
        // stitching
        for(size_t i=0; i<pfp.size(); i++) {
            for(size_t j=i+1; j<pfp.size(); j++) {
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
                    if(mc) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(pfp.at(i).G4ID==MCmuons.at(k).G4ID && pfp.at(j).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                if(pfp.at(i).slcID==pfp.at(j).slcID) {
                                    //hist.push_back(1);
                                    //hist.push_back(8);
                                }
                                else {
                                    //hist.push_back(5);
                                    //hist.push_back(13);
                                }
                                MCmuons.at(k).found = true;
                            }
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
                        if(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= param::slc::Tss) {
                            if(mc) {
                                for(size_t k=0; k<MCmuons.size(); k++) {
                                    if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                        //hist.push_back(1);
                                        //hist.push_back(9);
                                        //MCmuons.at(k).found = true;
                                    }
                                    //else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                    //         (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) && MCmuons.at(k).found == false)
                                    //            hist.push_back(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0]);
                                }
                            }
                            //continue;
                        }
                    }
                    else {
                        if(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= param::evt::Tss) {
                            if(mc) {
                                for(size_t k=0; k<MCmuons.size(); k++) {
                                    if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                        //hist.push_back(5);
                                        //hist.push_back(14);
                                        //MCmuons.at(k).found = true;
                                    }
                                }
                            }
                            //continue;
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
                    //if(phi*180/3.141592 >= param::slc::theta) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    //hist.push_back(1);
                                    //hist.push_back(10);
                                    //MCmuons.at(k).found = true;
                                }
                                //else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                //         (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) && MCmuons.at(k).found == false)
                                //            hist.push_back(phi*180/3.141592);
                            }
                        }
                        //continue;
                    //}
                }
                else {
                    if(phi*180/3.141592 >= param::evt::theta) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    //hist.push_back(5);
                                    //hist.push_back(15);
                                    //MCmuons.at(k).found = true;
                                }
                            }
                        }
                        //continue;
                    }
                }
                // check distance between tracks
                std::vector<float> length;
                length.push_back(Ra);
                length.push_back(Rb);
                std::sort(length.begin(), length.end());
                float dist = std::hypot(trk.at(1).start.x-trk.at(0).end.x, trk.at(1).start.y-trk.at(0).end.y, trk.at(1).start.z-trk.at(0).end.z);
                if(trk.at(0).slcID==trk.at(1).slcID) {
                    //if(dist/length[0] >= param::slc::Dss) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    //hist.push_back(1);
                                    //hist.push_back(11);
                                    //MCmuons.at(k).found = true;
                                }
                                else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                         (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) && MCmuons.at(k).found == false) {
                                            hist.push_back(dist/length[0]);
                                            MCmuons.at(k).found = true;
                                }
                            }
                        }
                        //continue;
                    //}
                }
                else {
                    //if(dist/length[0] >= param::evt::Dss) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    //hist.push_back(5);
                                    //hist.push_back(16);
                                    //hist.push_back(dist/length[0]);
                                    //MCmuons.at(k).found = true;
                                }
                                else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) ||
                                         (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) && MCmuons.at(k).found == false) {
                                            //hist.push_back(dist/length[0]);
                                            //MCmuons.at(k).found = true;
                                         }
                            }
                        }
                        //continue;
                    //}
                }
                // fill last bins
                if(mc) {
                    for(size_t k=0; k<MCmuons.size(); k++) {
                        if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found==false) {
                            //if(trk.at(0).slcID==trk.at(1).slcID) hist.push_back(0);
                            //else hist.push_back(4);
                            //MCmuons.at(k).found = true;
                        }
                        else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) || (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) && MCmuons.at(k).wrong==false) {
                            //if(trk.at(0).slcID==trk.at(1).slcID) hist.push_back(2);
                            //else hist.push_back(6);
                            //if(ordz[1] != vtx.z) continue;
                            //hist.push_back( (trk.at(0).len+trk.at(1).len)/std::max(trk.at(0).len, trk.at(1).len) );
                            //if(trk.at(0).len<50 && trk.at(0).len+trk.at(1).len<50) hist.push_back(0);
                            //else if(trk.at(0).len>50 && trk.at(0).len+trk.at(1).len>50) hist.push_back(1);
                            //else if(trk.at(0).len<50 && trk.at(0).len+trk.at(1).len>50) hist.push_back(2);
                            //else hist.push_back(3);
                            MCmuons.at(k).wrong = true;
                        }
                    }
                }
                nStitch++;
                stitch.push_back(nStitch);                           // number of stitching
                stitch.push_back(trk.at(0).slcID);                   // slcID 1
                stitch.push_back(trk.at(1).slcID);                   // slcID 2
                stitch.push_back(trk.at(0).id);                      // pfpID 1
                stitch.push_back(trk.at(1).id);                      // pfpID 2

                stitch.push_back(vtx.x);                             // vertex
                stitch.push_back(trk.at(0).start.x);                 // start 1
                stitch.push_back(trk.at(0).end.x);                   // end 1
                stitch.push_back(trk.at(1).start.x);                 // start 2
                stitch.push_back(trk.at(1).end.x);                   // end 2

                stitch.push_back(vtx.y);                             // vertex
                stitch.push_back(trk.at(0).start.y);                 // start 1
                stitch.push_back(trk.at(0).end.y);                   // end 1
                stitch.push_back(trk.at(1).start.y);                 // start 2
                stitch.push_back(trk.at(1).end.y);                   // end 2

                stitch.push_back(vtx.z);                             // vertex
                stitch.push_back(trk.at(0).start.z);                 // start 1
                stitch.push_back(trk.at(0).end.z);                   // end 1
                stitch.push_back(trk.at(1).start.z);                 // start 2
                stitch.push_back(trk.at(1).end.z);                   // end 2

                stitch.push_back(trk.at(0).len);                     // length before stitching
                stitch.push_back(trk.at(0).len + trk.at(1).len);     // length after stitching
                if(mc) for(size_t k=0; k<MCmuons.size(); k++) {
                    sthG4ID.push_back(MCmuons.at(k).G4ID);           // mu G4ID
                }
            }
        }
        //for(size_t i=0; i<MCmuons.size(); i++) if(MCmuons.at(i).found==true) hist.push_back(18);
        return stitch;
    }


    // stitching function
    const SpillMultiVar kStitch2([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> sth;
        std::vector<double> hist;
        std::vector<int> sthG4ID;
        bool mc = false;
        if(sr->hdr.ismc) mc = true;
        sth = stitch2(sr, hist, sthG4ID, mc);
        return hist;
    });

} // namespace ana

#endif // STITCH2_H