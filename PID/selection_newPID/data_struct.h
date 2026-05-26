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

    std::vector<double> dedx_ind1; //<-- DONE
    std::vector<double> rr_ind1; //<-- DONE
    std::vector<double> pitch_ind1; //<-- DONE
    std::vector<double> dedx_ind2; //<-- DONE
    std::vector<double> rr_ind2; //<-- DONE
    std::vector<double> pitch_ind2; //<-- DONE

    int bestplane = -1; //<-- DONE

    //int predicted_class = -1; //<-- DONE
    //std::vector<double> predicted_probabilities; //<-- DONE

    int id_ipfp_reco = -1; //<-- DONE
    
    //std::vector<double> hitx;
    //std::vector<double> hity;
    //std::vector<double> hitz;


    std::vector<double> dir; //<-- DONE

    double ek_coll = -1; //<-- DONE
    double ek_ind1 = -1; //<-- DONE
    double ek_ind2 = -1; //<-- DONE

    double ek_bestplane = -1; //<-- DONe

    double chi2_as_mu = -1;
    double chi2_as_pro = -1;

    double length = -1;
    double trackscore = -1;
    double start_distance_from_reco_vertex = -1;
    double end_distance_from_reco_vertex = -1;

    bool is_start_contained = false;
    bool is_end_contained = false;
    bool all_in_1_tpc = false;
    bool is_primary = false;

    double energy_deposited_as_pro = -1;
    double energy_deposited_as_pi  = -1;
    double energy_deposited_as_mu = -1;

    double shower_energy = -1;

    std::vector<double> start;
    std::vector<double> end;
    std::vector<double> likelihood_ratios;
    //std::vector<double> likelihood_ratios_mediumbinning;

    double daughter_depE = -1;
    double daughter_angle_end = -1;

    //int nhits_coll;

    double deposited_energy = -1;
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
    //std::string reco_slice_classification; //<-- DONE

    bool is_true_1muNp = false; //<-- DONE
    //bool is_reco_1muNp = false; //<-- DONE
    bool is_true_1mu0p0pi = false; //<-- DONE

    bool is_charged_current = false;
    bool is_neutral_current = false;
    int genie_mode=-1;
    bool is_in_FV = false;
    int neutrino_pdg;
    
    bool is_clear_cosmic = false;
    double true_neutrino_energy = 0;
    bool isinFV = false;
    double reco_true_vertex_distance = 0;
    bool all_tracks_contained = false;
    double barycenter=0; //This is z!!!
	double barycenter_x=0;
	double bar_charge=0;    
    double delta=0;
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
