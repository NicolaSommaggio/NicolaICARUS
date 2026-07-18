#ifndef SLICESTRUCT_H
#define SLICESTRUCT_H

#include <vector>
#include "TObject.h"

class _pfp : public TObject {
public:
    std::vector<double> _lr;
    double _depE = -999;
    std::vector<double> _KE;
    double _length = -999;
    std::vector<double> _daughter_vars;
    std::vector<double> _dedx;
    std::vector<double> _rr;
    double _theta_xw = -999;

    _pfp() {}
    virtual ~_pfp() {}

    ClassDef(_pfp,1);
};

class _slice : public TObject {
public:

    int _run = -999;
    int _evt = -999;
    int _slice_counter = -999;
    std::vector<_pfp> _protons;
    _pfp _mu;
    std::string _reco_class = "N/D";
    std::string _true_class = "N/D";
    double _nuE = -999;
    double _mu_pro_angle = -999;
    double _nu_score = -999;
  
    double _crlongtrkdiry = -999;
    double _bar_flash = -999;
    double _bar_flash_x = -999;
    
    _slice() {}
    virtual ~_slice() {}

    ClassDef(_slice,1);
};

#endif
