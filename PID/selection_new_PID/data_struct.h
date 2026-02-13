#ifndef DATASTRUCT_H
#define DATASTRUCT_H

#include <vector>
#include "TObject.h"

// =====================
// reco_track
// =====================
class reco_track : public TObject {
public:
    std::vector<double> dedx;
    std::vector<double> rr;
    std::vector<double> pitch;
    std::vector<double> hitx;
    std::vector<double> hity;
    std::vector<double> hitz;

    double length = 0;
    double trackscore = 0;
    double start_distance_from_reco_vertex = 0;
    double end_distance_from_reco_vertex = 0;

    bool is_start_contained = false;
    bool is_end_contained = false;
    bool all_in_1_tpc = false;
    bool is_primary = false;

    double energy_deposited_as_pro = 0;
    double energy_deposited_as_pi  = 0;

    std::vector<double> start;
    std::vector<double> end;
    std::vector<double> likelihood_ratios;

    double deposited_energy = 0;
    int ipfp = -1;

    reco_track() {}
    virtual ~reco_track() {}

    ClassDef(reco_track,1);
};

// =====================
// true_track
// =====================
class true_track : public TObject {
public:
    int ipfp = -1;
    int true_class = -1;
    int pdg = -1;
    
    true_track() {}
    virtual ~true_track() {}

    ClassDef(true_track,1);
};

// =====================
// track
// =====================
class track : public TObject {
public:
    reco_track reco;
    true_track truth;

    track() {}
    virtual ~track() {}

    ClassDef(track,1);
};

// =====================
// RecoSlice
// =====================
class RecoSlice : public TObject {
public:
    int slice_number = -1;
    int run = -1;
    int subrun = -1;
    int evt = -1;

    std::string true_slice_classifications;

    bool is_clear_cosmic = false;
    double true_neutrino_energy = 0;
    bool isinFV = false;
    double reco_true_vertex_distance = 0;
    bool all_tracks_contained = false;
    double light_charge_baricenter = 0;
    std::vector<double> reco_vertex;
    std::vector<double> true_vertex;

    double truth_matching_efficiency;
    double truth_matching_purity;

    std::vector<track> tracks;

    RecoSlice() {}
    virtual ~RecoSlice() {}

    ClassDef(RecoSlice,1);
};

#endif
