#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "testBDT.h"



using namespace ana;

void macro_testBDT(){

std::vector<double> features_prova = {5.47148e-01,  5.86870e-01,  5.68587e-01,  5.13775e-01,  7.09001e-01,
  9.97885e-02,  5.22486e-02, -7.35480e-02,  4.56859e-01, -4.81556e-02,
 -1.70308e-01,  3.95098e-01, -1.24626e-01,  4.25900e-01,  4.96458e-01,
  1.13418e+02, -1.00000e+00, -1.00000e+00};
std::vector<double> p_prova = model.predict_proba(features_prova);
for(const auto &p : p_prova){cout << p << " ";}
cout << endl;
  
//1/2 MC OVERLAYS GRAY
const std::string fdata = "/exp/icarus/data/users/nsommagg/SELECTION_FOLDER/*.root";

SpectrumLoader loader(fdata);         

const Binning kBinz = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, testBDT ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

}