#ifndef DAUGHTERSINFO_H
#define DAUGHTERSINFO_H

#include <vector>
#include "Rtypes.h"  // Per ClassDef

struct daughtersInfo {
    std::vector<double> dEdx;
    std::vector<double> rr;
    int pdg;

    daughtersInfo() : pdg(0) {}  // costruttore default
    virtual ~daughtersInfo();     // distruttore virtuale

    ClassDef(daughtersInfo,1)    // ROOT macro per la dictionary
};

#endif
