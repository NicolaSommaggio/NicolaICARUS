////////////////////////////////////////////////////////////////////
// Definition of the stitching function between slices
//
// Author: Alessandro Maria Ricci
////////////////////////////////////////////////////////////////////

#ifndef STITCH_H
#define STITCH_H

#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202412.h"

namespace ana {

    namespace param {
        bool print = false;      // print on terminal
        bool EHcomp = false;     // enable quality cut on hit_comp and E_comp
        bool checkVtx = false;   // check if the vertex is between the two tracks
        float mcLmin = 50.;      // minimal length of the track on mc
        namespace slc {
            float Tss = 0.7;     // overlap of the tracks on the z-axis
            float theta = 35.;   // angle between the directions of the tracks
            float Dss = 0.7;     // 3D distance between the tracks (ratio)
        }
        namespace evt {
            float Tss = 0.7;     // overlap of the tracks on the z-axis
            float theta = 30.;   // angle between the directions of the tracks
            float Dss = 0.9;     // 3D distance between the tracks (ratio)
        }
    } // namespace param

    struct muon {
        int G4ID;
        size_t counter;
        bool found;
        bool wrong;
        bool split;
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

    // stitching function in the event
    std::vector<double> stitch(const caf::SRSpillProxy* sr, std::vector<double> &hist, std::vector<int> &sthG4ID, bool mc) {
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
            if(!kIcarus202401RecoFiducial(sr->slc.at(i))) continue;
            for(size_t j=0; j<sr->slc.at(i).reco.npfp; j++) {
                // check if the pfp is a muon
                if(!kIcarus202401MuonTrack(sr->slc.at(i).reco.pfp.at(j))) continue;
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
                    P.energy_comp = sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.energy_completeness;
                    P.hit_comp = sr->slc.at(i).reco.pfp.at(j).trk.truth.bestmatch.hit_completeness;
                    // quality cut on hit_comp and E_comp
                    if(param::EHcomp) {
                        for(size_t k=0; k<MCmuons.size(); k++) {
                            if(MCmuons.at(k).G4ID == P.G4ID && P.energy_comp >= 0.1 && P.hit_comp >= 0.1) {
                                MCmuons.at(k).counter++;
                                if(MCmuons.at(k).counter >= 2) {
                                    MCmuons.at(k).split = true;
                                }
                            }
                        }
                    }
                }
                pfp.push_back(P);
            }
        }
        // quality cut on hit_comp and E_comp
        if(mc && param::EHcomp) {
            // delete all elements that satisfy the condition split = false
            MCmuons.erase(std::remove_if(MCmuons.begin(), MCmuons.end(), [](const muon& mu) {
                return mu.split == false;
            }), MCmuons.end());
        }
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
                                    hist.push_back(1);
                                    hist.push_back(8);
                                }
                                else {
                                    hist.push_back(5);
                                    hist.push_back(13);
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
                    /*if(trk.at(0).len<30 || trk.at(1).len<30) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(pfp.at(i).G4ID==MCmuons.at(k).G4ID && pfp.at(j).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    if(pfp.at(i).slcID==pfp.at(j).slcID) {
                                        hist.push_back(1);
                                        hist.push_back(8);
                                    }
                                    else {
                                        hist.push_back(5);
                                        hist.push_back(13);
                                    }
                                    MCmuons.at(k).found = true;
                                }
                            }
                        }
                        continue;
                    }*/
                    // sort start and end end point of each track according to the 3D distance and z-axis
                    Icarus202503SortPFP(trk.at(0), trk.at(1));
                }
                else {
                    trk.push_back(pfp.at(j));
                    trk.push_back(pfp.at(i));
                    /*if(trk.at(0).len<30 || trk.at(1).len<30) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(pfp.at(i).G4ID==MCmuons.at(k).G4ID && pfp.at(j).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    if(pfp.at(i).slcID==pfp.at(j).slcID) {
                                        hist.push_back(1);
                                        hist.push_back(8);
                                    }
                                    else {
                                        hist.push_back(5);
                                        hist.push_back(13);
                                    }
                                    MCmuons.at(k).found = true;
                                }
                            }
                        }
                        continue;
                    }*/
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
                                        hist.push_back(1);
                                        hist.push_back(9);
                                        MCmuons.at(k).found = true;
                                    }
                                }
                            }
                            continue;
                        }
                    }
                    else {
                        if(std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= param::evt::Tss) {
                            if(mc) {
                                for(size_t k=0; k<MCmuons.size(); k++) {
                                    if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                        hist.push_back(5);
                                        hist.push_back(14);
                                        MCmuons.at(k).found = true;
                                    }
                                }
                            }
                            continue;
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
                    if(phi*180/3.141592 >= param::slc::theta) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    hist.push_back(1);
                                    hist.push_back(10);
                                    MCmuons.at(k).found = true;
                                }
                            }
                        }
                        continue;
                    }
                }
                else {
                    if(phi*180/3.141592 >= param::evt::theta) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    hist.push_back(5);
                                    hist.push_back(15);
                                    MCmuons.at(k).found = true;
                                }
                            }
                        }
                        continue;
                    }
                }
                // check distance between tracks
                std::vector<float> length;
                length.push_back(Ra);
                length.push_back(Rb);
                std::sort(length.begin(), length.end());
                float dist = std::hypot(trk.at(1).start.x-trk.at(0).end.x, trk.at(1).start.y-trk.at(0).end.y, trk.at(1).start.z-trk.at(0).end.z);
                if(trk.at(0).slcID==trk.at(1).slcID) {
                    if(dist/length[0] >= param::slc::Dss) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    hist.push_back(1);
                                    hist.push_back(11);
                                    MCmuons.at(k).found = true;
                                }
                            }
                        }
                        continue;
                    }
                }
                else {
                    if(dist/length[0] >= param::evt::Dss) {
                        if(mc) {
                            for(size_t k=0; k<MCmuons.size(); k++) {
                                if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found == false) {
                                    hist.push_back(5);
                                    hist.push_back(16);
                                    MCmuons.at(k).found = true;
                                }
                            }
                        }
                        continue;
                    }
                }
                // fill last bins
                if(mc) {
                    for(size_t k=0; k<MCmuons.size(); k++) {
                        if(trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID && MCmuons.at(k).found==false) {
                            if(trk.at(0).slcID==trk.at(1).slcID) hist.push_back(0);
                            else hist.push_back(4);
                            MCmuons.at(k).found = true;
                        }
                        else if( (trk.at(0).G4ID==MCmuons.at(k).G4ID && trk.at(1).G4ID!=MCmuons.at(k).G4ID) || (trk.at(0).G4ID!=MCmuons.at(k).G4ID && trk.at(1).G4ID==MCmuons.at(k).G4ID) && MCmuons.at(k).wrong==false) {
                            if(trk.at(0).slcID==trk.at(1).slcID) hist.push_back(2);
                            else hist.push_back(6);
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
        for(size_t i=0; i<MCmuons.size(); i++) if(MCmuons.at(i).found==true) hist.push_back(18);
        return stitch;
    }


    // print on terminal the output of stitch
    void PrintStitch(std::vector<double> &stitch) {
        for(size_t i=0; i<stitch.size(); i=i+22) {
            std::cout << " " << std::endl;
            std::cout << "===================" << std::endl;
            std::cout << "Stitching No: " << stitch.at(0+i) << std::endl;
            std::cout << " " << std::endl;
            std::cout << "slcID 1 = " << stitch.at(1+i) << std::endl;
            std::cout << "slcID 2 = " << stitch.at(2+i) << std::endl;
            std::cout << "pfpID 1 = " << stitch.at(3+i) << std::endl;
            std::cout << "pfpID 2 = " << stitch.at(4+i) << std::endl;
            std::cout << " " << std::endl;
            std::cout << "geom = {v,s1,e1,s2,e2}" << std::endl;
            for(size_t j=0; j<5; j++) {
                std::cout << "x[" << j << "] = " << stitch.at(5+j+i) <<
                    "\t y[" << j << "] = " << stitch.at(10+j+i) <<
                    "\t z[" << j << "] = " << stitch.at(15+j+i) << std::endl;
            }
            std::cout << " " << std::endl;
            std::cout << "L1 = " << stitch.at(20+i) << std::endl;
            std::cout << "L1+L2 = " << stitch.at(21+i) << std::endl;
            std::cout << " " << std::endl;
        }
    }


    // stitching function
    const SpillMultiVar kStitch([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> sth;
        std::vector<double> hist;
        std::vector<int> sthG4ID;
        bool mc = false;
        if(sr->hdr.ismc) mc = true;
        sth = stitch(sr, hist, sthG4ID, mc);
        if(!sth.empty() && param::print) {
            std::cout << " " << std::endl;
            std::cout << " " << std::endl;
            std::cout << "==============" << std::endl;
            std::cout << "Run: " << sr->hdr.run << std::endl;
            std::cout << "Subrun: " << sr->hdr.subrun << std::endl;
            std::cout << "Event: " << sr->hdr.evt << std::endl;
            std::cout << "==============" << std::endl;
            PrintStitch(sth);
        }
        return hist;
    });



    // adapted from kIcarus202401MuonIdx in NumuVarsIcarus202401.cxx
    bool kIcarus202503MuonTrack(const caf::SRPFPProxy &pfp, size_t &counter) {
        // The (dis)qualification of a slice is based upon the track level information.
        bool muTrack = false;
        if(pfp.trackScore < 0.4) return muTrack;
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
        if(trk.len >= 20 && counter<4) counter=4;
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len >= 20 );
        if(MaybeMuonContained) muTrack = true;
        return muTrack;
    }


    // return the number of mu, primary daugher of nuMuCC interactions, contained in a single TPC volume
    const SpillMultiVar kMu([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> mu;
        for(size_t i=0; i<sr->mc.nnu; i++) {
            if(!sr->mc.nu.at(i).iscc) continue;
            for(int j=0; j<sr->mc.nu.at(i).nprim; j++) {
                if( !(sr->mc.nu.at(i).prim.at(j).pdg == 13 || sr->mc.nu.at(i).prim.at(j).pdg == -13)) continue;
                if(!sr->mc.nu.at(i).prim.at(j).contained) continue;
                if(sr->mc.nu.at(i).prim.at(j).length < param::mcLmin) continue;
                int MCG4ID = sr->mc.nu.at(i).prim.at(j).G4ID;
                size_t counter = 0;
                // loop on the slices
                for(int k=0; k<sr->nslc; k++) {
                    // check if the slice is in the fiducial volume
                    if(!kIcarus202401RecoFiducial(sr->slc.at(k))) continue;
                    // loop on the pfp particles
                    for(size_t m=0; m<sr->slc.at(k).reco.npfp; m++) {
                        int pfpG4ID = sr->slc.at(k).reco.pfp.at(m).trk.truth.bestmatch.G4ID;
                        if(MCG4ID != pfpG4ID) continue;
                        if(counter<1) counter=1; // muons that passed the RecoFiducial cut
                        // check if the pfp is a muon
                        if(!kIcarus202503MuonTrack(sr->slc.at(k).reco.pfp.at(m), counter)) continue;
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