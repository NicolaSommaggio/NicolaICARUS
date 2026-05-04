///////////////////////////////////////////////////////////////////////
// Function to compute the track momentum
// adapted from LArSoft/larreco/RecoAlg/TrackMomentumCalculator.cxx
//
// Author: Alessandro Maria Ricci
///////////////////////////////////////////////////////////////////////

#ifndef TRACK_MOMENTUM_CALCULATOR_H
#define TRACK_MOMENTUM_CALCULATOR_H

// C++
#include <array>
#include <iostream>

// ROOT
#include "TGraph.h"
#include "TSpline.h"

#include "Selection.h"

namespace ana {

  constexpr double LAr_density{1.396};

  constexpr auto range_gramper_cm() {
    std::array<double, 73> Range_grampercm{
      {0.9833,   1.36,      1.786,    2.507,   3.321,   4.859,   6.598,   8.512,   10.58,   12.78,
       15.1,     17.52,     20.04,    25.31,   30.84,   36.59,   42.5,    54.73,   67.32,   86.66,
       106.3,    139.4,     172.5,    205.6,   238.5,   271.1,   303.5,   335.7,   367.7,   431.0,
       493.4,    555.2,     616.3,    736.8,   855.2,   1030.0,  1202.0,  1482.0,  1758.0,  2029.0,
       2297.0,   2562.0,    2825.0,   3085.0,  3343.0,  3854.0,  4359.0,  4859.0,  5354.0,  6333.0,
       7298.0,   8726.0,    10130.0,  12430.0, 14690.0, 16920.0, 19100.0, 21260.0, 23380.0, 25480.0,
       27550.0,  31610.0,   35580.0,  39460.0, 43260.0, 50620.0, 57680.0, 67780.0, 77340.0, 92220.0,
       1.06e+05, 1.188e+05, 1.307e+05}};
    for (double& value : Range_grampercm) {
      value /= LAr_density; // convert to cm
    }
    return Range_grampercm;
  }

  constexpr auto Range_grampercm = range_gramper_cm();
  constexpr std::array<double, 73> KE_MeV{
    {10.0,    12.0,    14.0,    17.0,    20.0,    25.0,    30.0,    35.0,    40.0,    45.0,
     50.0,    55.0,    60.0,    70.0,    80.0,    90.0,    100.0,   120.0,   140.0,   170.0,
     200.0,   250.0,   300.0,   350.0,   400.0,   450.0,   500.0,   550.0,   600.0,   700.0,
     800.0,   900.0,   1000.0,  1200.0,  1400.0,  1700.0,  2000.0,  2500.0,  3000.0,  3500.0,
     4000.0,  4500.0,  5000.0,  5500.0,  6000.0,  7000.0,  8000.0,  9000.0,  10000.0, 12000.0,
     14000.0, 17000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, 45000.0, 50000.0, 55000.0,
     60000.0, 70000.0, 80000.0, 90000.0, 1e+05,   1.2e+05, 1.4e+05, 1.7e+05, 2e+05,   2.5e+05,
     3e+05,   3.5e+05, 4e+05}};
  TGraph const KEvsR{73, Range_grampercm.data(), KE_MeV.data()};
  TSpline3 const KEvsR_spline3{"KEvsRS", &KEvsR};

  constexpr std::array<double, 73> KE_MeV_Pi{
    {11.1928, 13.4567, 15.6946,  19.0487,  22.3945,  27.9408,  33.4579,  38.9494,  44.4198,
     49.8612, 55.2838, 60.6757,  66.0632,  76.7631,  87.3839,  97.9617,  108.461,  129.341,
     150.023, 180.793, 211.275,  261.704,  311.585,  361.284,  410.705,  459.809,  508.805,
     557.727, 606.592, 704.011,  801.047,  898.111,  994.983,  1188.51,  1381.61,  1671.44,
     1961.46, 2442.32, 2925.25,  3406.69,  3888.88,  4370.85,  4853.71,  5335.02,  5816.14,
     6778.21, 7739.41, 8699.89,  9658.32,  11572.9,  13480.8,  16335.1,  19171.0,  23866.9,
     28529.9, 33169.1, 37734.7,  42284.0,  46770.4,  51233.2,  55648.5,  64349.8,  72904.4,
     81303.4, 89561.4, 1.06e+05, 1.21e+05, 1.43e+05, 1.65e+05, 1.98e+05, 2.28e+05, 2.55e+05,
     2.77e+05}};
  TGraph const KEvsRPi{73, Range_grampercm.data(), KE_MeV_Pi.data()};
  TSpline3 const KEvsRPi_spline3{"KEvsRS_pion", &KEvsRPi};

  double GetTrackMomentum(double trkrange, int pdg) {
    /* Muon range-momentum tables from CSDA (Argon density = 1.4 g/cm^3)
       website:
       http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf

       CSDA table values:
       double Range_grampercm[30] = {9.833E-1, 1.786E0, 3.321E0,
       6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
       2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3,
       4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
       4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5}; double KE_MeV[30] = {10, 14,
       20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000,
       4000, 8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000,
       200000, 300000, 400000};

       Functions below are obtained by fitting polynomial fits to KE_MeV vs
       Range (cm) graph. A better fit was obtained by splitting the graph into
       two: Below range<=200cm,a polynomial of power 4 was a good fit; above
       200cm, a polynomial of power 6 was a good fit

       Fit errors for future purposes:
       Below 200cm, Forpoly4 fit: p0 err=1.38533;p1 err=0.209626; p2
       err=0.00650077; p3 err=6.42207E-5; p4 err=1.94893E-7; Above 200cm,
       Forpoly6 fit: p0 err=5.24743;p1 err=0.0176229; p2 err=1.6263E-5; p3
       err=5.9155E-9; p4 err=9.71709E-13; p5 err=7.22381E-17;p6
       err=1.9709E-21;*/

    ///////////////////////////////////////////////////////////////////////////
    //*********For muon, the calculations are valid up to 1.91E4 cm range
    //corresponding to a Muon KE of 40 GeV**********//
    ///////////////////////////////////////////////////////////////////////////

    /*Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
      website: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html

      CSDA values:
      double KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300,
      350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
      1500, 2000, 2500, 3000, 4000, 5000};

      double Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296,
      2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1,
      9.184E1, 1.138E2, 1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2,
      2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
      7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3};

      Functions below are obtained by fitting power and polynomial fits to
      KE_MeV vs Range (cm) graph. A better fit was obtained by splitting the
      graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
      polynomial of power 6 was a good fit

      Fit errors for future purposes:
      For power function fit: a=0.388873; and b=0.00347075
      Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3
      err=4.45542E-7; p4 err=4.16428E-10; p5 err=1.81679E-13;p6
      err=2.96958E-17;*/

    ///////////////////////////////////////////////////////////////////////////
    //*********For proton, the calculations are valid up to 3.022E3 cm range
    //corresponding to a Muon KE of 5 GeV**********//
    ///////////////////////////////////////////////////////////////////////////

    /* Pion range-momentum extracted from Bethe-Bloch equation.  The
       computation are all one using
       https://github.com/sungbinoh/PDSPTreeAnalyzer.git (last commit 4a3d1d2)

       The computation is done with the spline of the Range_grampercm and
       KE_MeV_Pion.*/

    ///////////////////////////////////////////////////////////////////////////
    //*********For pion, the calculations are valid only when there is no
    //inelastic scattering present. So this method should be used
    //carefully.**********//
    ///////////////////////////////////////////////////////////////////////////

    if (trkrange < 0 || std::isnan(trkrange)) {
      std::cout << "TrackMomentumCalculator:" << "Invalid track range " << trkrange << " return -1" << std::endl;
      return -1.;
    }

    double KE, Momentum, M;
    constexpr double Muon_M = mmu * 1e3, Proton_M = 938.272, Pion_M = 139.57039;

    if (abs(pdg) == 13) {
      M = Muon_M;
      KE = KEvsR_spline3.Eval(trkrange);
    }
    else if (abs(pdg) == 2212) {
      M = Proton_M;
      if (trkrange > 0 && trkrange <= 80)
        KE = 29.9317 * std::pow(trkrange, 0.586304);
      else if (trkrange > 80 && trkrange <= 3.022E3)
        KE = 149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
             (4.34587E-6 * trkrange * trkrange * trkrange) +
             (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
             (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
             (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange * trkrange);
      else
        KE = -999;
    }
    else if (abs(pdg) == 211) {
      M = Pion_M;
      KE = KEvsRPi_spline3.Eval(trkrange);
    }
    else
      KE = -999;

    if (KE < 0)
      Momentum = -999;
    else
      Momentum = std::sqrt((KE * KE) + (2 * M * KE));

    Momentum = Momentum / 1000;

    return Momentum;
  }

} // namespace ana

#endif // TRACK_MOMENTUM_CALCULATOR_H