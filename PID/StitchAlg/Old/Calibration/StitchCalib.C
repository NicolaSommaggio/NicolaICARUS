////////////////////////////////////////////////////////////////////
// Macro to stitch the split tracks of muons
//
// Author: Alessandro Maria Ricci
////////////////////////////////////////////////////////////////////

// C++
#include <string>

// ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

// sbnana
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "StitchCalib.h"
#include "StitchCalib2.h"

using namespace ana;

// MACRO
void StitchCalib() {

    // Target files
    //const std::string TargetFiles = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_CV_2ndV/mc/reconstructed/icaruscode_v09_89_01_01/flatcaf/*/*/*.root";
    //const std::string TargetFiles = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_IndGap1WireFil_2ndV/mc/reconstructed/icaruscode_v09_89_01_02/flatcaf/*/*/*.root";
    //const std::string TargetFiles = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_CathBen_2ndV/mc/reconstructed/icaruscode_v09_89_02_00p01/flatcaf/*/*/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/Complete_MC_final/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/Complete_MC_final/mc-v09_84_00_01-202403-cnaf-corrsce_run3*.concatflat.caf.root";
    const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/CAFs_9435_prod_CNAF/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/cafs_Prescaled_Run2_v0984/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/Offbeam_Unblind_DATA_bnbmaj/*.root";

    // Note: the first letter in the object name identifies the type

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum
    SpectrumLoader NuLoader(TargetFiles);
    std::cout << "Target files loaded." << std::endl;

    // Create the binning schemes for the Vars we wish to plot
    const Binning bMu = Binning::Simple(6, 0, 6);
    //const Binning bStitch = Binning::Simple(19, 0, 19);
    const Binning bStitch = Binning::Simple(40, 0, 4);
    //const Binning bStitch = Binning::Simple(36, 0, 180);
    //const Binning bStitch = Binning::Simple(10, 1, 2);
    //const Binning bStitch = Binning::Simple(4, 0, 4);

    // Create a Spectrum showing the number of mu, primary daughters of nuMuCC interactions, contained in a single TPC volume
    Spectrum sMu( "Nmu",
                  bMu,
                  NuLoader,
                  kMu,
                  kNoSpillCut );

    // Create a Spectrum showing how the algorithm stitches
    Spectrum sStitch( "sStitch",                    // A label for the Spectrum
		              bStitch,                      // Define the binning
		              NuLoader,                     // Associate this Spectrum with the Loader object (and its target CAF)
		              kStitch,                      // The Var to plot
		              kNoSpillCut );                // The SpillCut to use (none in this case)

    // Create a Spectrum showing how the algorithm stitches
    Spectrum sStitch2( "sStitch2",                  // A label for the Spectrum
                       bStitch,                     // Define the binning
                       NuLoader,                    // Associate this Spectrum with the Loader object (and its target CAF)
                       kStitch2,                    // The Var to plot
                       kNoSpillCut );               // The SpillCut to use (none in this case)

    // Now that each Spectrum is defined, use the Go() method to populate the Spectrum objects
    NuLoader.Go();
    std::cout << "Spectra created and populated." << std::endl;

    // ---- DRAW -----
    // Write the Spectrum objects to a file as a TH1 object.
    TFile fStitch("Stitch.root", "recreate");

    TH1* hMu = sMu.ToTH1(sMu.Livetime(), kLivetime);
    hMu->SetStats(true);
    hMu->Write("Nmu");

    TH1* hStitch = sStitch.ToTH1(sStitch.Livetime(), kLivetime);
    hStitch->SetStats(true);
    hStitch->SetLineColor(kRed);
    hStitch->Scale(1./hStitch->Integral());
    hStitch->Write("Stitch");

    TH1* hStitch2 = sStitch2.ToTH1(sStitch2.Livetime(), kLivetime);
    hStitch2->SetStats(true);
    hStitch2->SetLineColor(kBlue);
    hStitch2->Scale(1./hStitch2->Integral());
    hStitch2->Write("Stitch2");

    TCanvas* cOverlap = new TCanvas();
    hStitch->Draw("hist");
    hStitch2->Draw("SameHist");
    cOverlap->Write("Overlap");

    fStitch.Close();
    std::cout << "Spectra saved." << std::endl;
}
