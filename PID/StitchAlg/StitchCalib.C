///////////////////////////////////////////////////////////////////////////
// Macro to calibrate the algorithm to stitch the split tracks of muons
//
// Author: Alessandro Maria Ricci
///////////////////////////////////////////////////////////////////////////

// C++
#include <string>

// ROOT
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"

// sbnana
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "Spectra.h"

using namespace ana;

// MACRO
void StitchCalib() {

    // Target files

    // MC with 1D deconvolution
    //const std::string TargetFiles = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_CV_2ndV/mc/reconstructed/icaruscode_v09_89_01_01/flatcaf/*/*/*.root";
    //const std::string TargetFiles = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_IndGap1WireFil_2ndV/mc/reconstructed/icaruscode_v09_89_01_02/flatcaf/*/*/*.root";
    //const std::string TargetFiles = "/pnfs/sbn/data/sbn_fd/poms_production/2024A_ICARUS_MC_Sys_NuCos/2024A_MC_Sys_NuCos_CathBen_2ndV/mc/reconstructed/icaruscode_v09_89_02_00p01/flatcaf/*/*/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/Complete_MC_final/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/Complete_MC_final/mc-v09_84_00_01-202403-cnaf-corrsce_run3*.concatflat.caf.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/CAFs_9435_prod_CNAF/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/cafs_Prescaled_Run2_v0984/*.root";
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/Offbeam_Unblind_DATA_bnbmaj/*.root";
    
    // Overlays with 2D deconvolution
    const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026_[2,3]/*flat.caf.root";

    // Real data with 2D deconvolution
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/Data_prescaled_Run2_calibrated/*.root";

    // CNAF
    // 30% del Run 2
    //const std::string TargetFiles = "/storage/gpfs_data/icarus/plain/data/prod/run2-v09_84_00_01-202403-cnaf/concat-caf/control-sample-30/*.root";
    // Run 9435 with 2D deconvolution
    //const std::string TargetFiles = "/storage/gpfs_data/icarus/local/users/cfarnese/9435_spring_production_2025/9435_2025_FNAL/cafs_calibrated/*.root";
    // Run 9435 with 1D deconvolution
    //const std::string TargetFiles = "/storage/gpfs_data/icarus/local/users/cfarnese/9435_1muNp_icarusonly/9435_2025_FNAL_2/*.root";

    // Note: the first letter in the object name identifies the type

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum
    SpectrumLoader NuLoader(TargetFiles);
    std::cout << "Target files loaded." << std::endl;

    // Create the binning schemes for the Vars we wish to plot
    const Binning bTssIS = Binning::Simple(50, 0, 5);
    const Binning bTssES = Binning::Simple(25, 0, 5);
    const Binning bThetaIS = Binning::Simple(36, 0, 180);
    const Binning bThetaES = Binning::Simple(36, 0, 180);
    const Binning bDssIS = Binning::Simple(50, 0, 5);
    const Binning bDssES = Binning::Simple(50, 0, 5);

    //===================================================================================================
    //===================================================================================================

    // Create a Spectrum showing the values of the parameter Tss for intra-slice muon-muon stitching in MC sample
    Spectrum sTssMuMuIS( "TssMuMuIS",
		                 bTssIS,
		                 NuLoader,
		                 kTssMuMuIS,
		                 kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Tss for intra-slice muon-other stitching in MC sample
    Spectrum sTssMuOtherIS( "TssMuOtherIS",
		                    bTssIS,
		                    NuLoader,
		                    kTssMuOtherIS,
		                    kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Tss for intra-slice stitching in the real data
    Spectrum sTssDataIS( "TssDataIS",
                         bTssIS,
                         NuLoader,
                         kTssDataIS,
                         kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Tss for extra-slice muon-muon stitching in MC sample
    Spectrum sTssMuMuES( "TssMuMuES",
		                 bTssES,
		                 NuLoader,
		                 kTssMuMuES,
		                 kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Tss for extra-slice muon-other stitching in MC sample
    Spectrum sTssMuOtherES( "TssMuOtherES",
		                    bTssES,
		                    NuLoader,
		                    kTssMuOtherES,
		                    kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Tss for extra-slice stitching in the real data
    Spectrum sTssDataES( "TssDataES",
                         bTssES,
                         NuLoader,
                         kTssDataES,
                         kNoSpillCut );

    // Create a Spectrum showing the values of the parameter theta for intra-slice muon-muon stitching in MC sample
    Spectrum sThetaMuMuIS( "ThetaMuMuIS",
		                   bThetaIS,
		                   NuLoader,
		                   kThetaMuMuIS,
		                   kNoSpillCut );

    // Create a Spectrum showing the values of the parameter theta for intra-slice muon-other stitching in MC sample
    Spectrum sThetaMuOtherIS( "ThetaMuOtherIS",
		                      bThetaIS,
		                      NuLoader,
		                      kThetaMuOtherIS,
		                      kNoSpillCut );

    // Create a Spectrum showing the values of the parameter theta for intra-slice stitching in the real data
    Spectrum sThetaDataIS( "ThetaDataIS",
                           bThetaIS,
                           NuLoader,
                           kThetaDataIS,
                           kNoSpillCut );

    // Create a Spectrum showing the values of the parameter theta for extra-slice muon-mu stitching in MC sample
    Spectrum sThetaMuMuES( "ThetaMuMuES",
		                   bThetaES,
		                   NuLoader,
		                   kThetaMuMuES,
		                   kNoSpillCut );

    // Create a Spectrum showing the values of the parameter theta for extra-slice muon-other stitching in MC sample
    Spectrum sThetaMuOtherES( "ThetaMuOtherES",
		                      bThetaES,
		                      NuLoader,
		                      kThetaMuOtherES,
		                      kNoSpillCut );

    // Create a Spectrum showing the values of the parameter theta for extra-slice stitching in the real data
    Spectrum sThetaDataES( "ThetaDataES",
                           bThetaES,
                           NuLoader,
                           kThetaDataES,
                           kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Dss for intra-slice muon-muon stitching in MC sample
    Spectrum sDssMuMuIS( "DssMuMuIS",
		                 bDssIS,
		                 NuLoader,
		                 kDssMuMuIS,
		                 kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Dss for intra-slice muon-other stitching in MC sample
    Spectrum sDssMuOtherIS( "DssMuOtherIS",
                            bDssIS,
                            NuLoader,
                            kDssMuOtherIS,
                            kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Dss for intra-slice stitching in the real data
    Spectrum sDssDataIS( "DssDataIS",
                         bDssIS,
                         NuLoader,
                         kDssDataIS,
                         kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Dss for extra-slice muon-muon stitching in MC sample
    Spectrum sDssMuMuES( "DssMuMuES",
                         bDssES,
                         NuLoader,
                         kDssMuMuES,
                         kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Dss for extra-slice muon-other stitching in MC sample
    Spectrum sDssMuOtherES( "DssMuOtherES",
                            bDssES,
                            NuLoader,
                            kDssMuOtherES,
                            kNoSpillCut );

    // Create a Spectrum showing the values of the parameter Dss for extra-slice stitching in the real data
    Spectrum sDssDataES( "DssDataES",
                         bDssES,
                         NuLoader,
                         kDssDataES,
                         kNoSpillCut );

    //===================================================================================================
    //===================================================================================================

    // Now that each Spectrum is defined, use the Go() method to populate the Spectrum objects
    NuLoader.Go();
    std::cout << "Spectra created and populated." << std::endl;

    // ---- DRAW -----
    // Write the Spectrum objects to a file as a TH1 object.
    TFile fStitch("StitchCalib_newPID.root", "recreate");

    //==================================================================================

    TH1* hTssMuMuIS = sTssMuMuIS.ToTH1(sTssMuMuIS.POT(), kPOT);
    hTssMuMuIS->SetLineColor(kRed);
    hTssMuMuIS->Scale(1./hTssMuMuIS->Integral());
    hTssMuMuIS->SetTitle("MC Distribution of Tss for intra-slice stitching");
    hTssMuMuIS->GetXaxis()->SetTitle("Tss");
    hTssMuMuIS->GetXaxis()->SetTitleSize(0.045);
    hTssMuMuIS->GetYaxis()->SetTitle("Distribution [%]");
    hTssMuMuIS->GetYaxis()->SetTitleSize(0.045);

    TH1* hTssMuOtherIS = sTssMuOtherIS.ToTH1(sTssMuOtherIS.POT(), kPOT);
    hTssMuOtherIS->SetLineColor(kBlue);
    hTssMuOtherIS->Scale(1./hTssMuOtherIS->Integral());

    TPaveStats* stTssMuMuIS = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stTssMuMuIS->AddText(Form("Entries = %.0f", hTssMuMuIS->GetEntries()));
    stTssMuMuIS->AddText(Form("Mean = %.3f", hTssMuMuIS->GetMean()));
    stTssMuMuIS->AddText(Form("StdDev = %.3f", hTssMuMuIS->GetStdDev()));
    stTssMuMuIS->SetTextColor(kRed);
    hTssMuMuIS->GetListOfFunctions()->Add(stTssMuMuIS);

    TPaveStats* stTssMuOtherIS = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stTssMuOtherIS->AddText(Form("Entries = %.0f", hTssMuOtherIS->GetEntries()));
    stTssMuOtherIS->AddText(Form("Mean = %.3f", hTssMuOtherIS->GetMean()));
    stTssMuOtherIS->AddText(Form("StdDev = %.3f", hTssMuOtherIS->GetStdDev()));
    stTssMuOtherIS->SetTextColor(kBlue);
    hTssMuOtherIS->GetListOfFunctions()->Add(stTssMuOtherIS);

    TCanvas* cTssISCalib = new TCanvas();
    hTssMuMuIS->Draw("hist");
    hTssMuOtherIS->Draw("SameHist");
    cTssISCalib->Write("MC calib: Tss intra-slice overlap");

    //==================================================================================

    TH1* hTssMuMuES = sTssMuMuES.ToTH1(sTssMuMuES.POT(), kPOT);
    hTssMuMuES->SetLineColor(kRed);
    hTssMuMuES->Scale(1./hTssMuMuES->Integral());
    hTssMuMuES->SetTitle("MC distribution of Tss for extra-slice stitching");
    hTssMuMuES->GetXaxis()->SetTitle("Tss");
    hTssMuMuES->GetXaxis()->SetTitleSize(0.045);
    hTssMuMuES->GetYaxis()->SetTitle("Distribution [%]");
    hTssMuMuES->GetYaxis()->SetTitleSize(0.045);

    TH1* hTssMuOtherES = sTssMuOtherES.ToTH1(sTssMuOtherES.POT(), kPOT);
    hTssMuOtherES->SetLineColor(kBlue);
    hTssMuOtherES->Scale(1./hTssMuOtherES->Integral());

    TPaveStats* stTssMuMuES = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stTssMuMuES->AddText(Form("Entries = %.0f", hTssMuMuES->GetEntries()));
    stTssMuMuES->AddText(Form("Mean = %.3f", hTssMuMuES->GetMean()));
    stTssMuMuES->AddText(Form("StdDev = %.3f", hTssMuMuES->GetStdDev()));
    stTssMuMuES->SetTextColor(kRed);
    hTssMuMuES->GetListOfFunctions()->Add(stTssMuMuES);

    TPaveStats* stTssMuOtherES = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stTssMuOtherES->AddText(Form("Entries = %.0f", hTssMuOtherES->GetEntries()));
    stTssMuOtherES->AddText(Form("Mean = %.3f", hTssMuOtherES->GetMean()));
    stTssMuOtherES->AddText(Form("StdDev = %.3f", hTssMuOtherES->GetStdDev()));
    stTssMuOtherES->SetTextColor(kBlue);
    hTssMuOtherES->GetListOfFunctions()->Add(stTssMuOtherES);

    TCanvas* cTssESCalib = new TCanvas();
    hTssMuMuES->Draw("hist");
    hTssMuOtherES->Draw("SameHist");
    cTssESCalib->Write("MC calib: Tss extra-slice overlap");

    //==================================================================================

    TH1* hThetaMuMuIS = sThetaMuMuIS.ToTH1(sThetaMuMuIS.POT(), kPOT);
    hThetaMuMuIS->SetLineColor(kRed);
    hThetaMuMuIS->Scale(1./hThetaMuMuIS->Integral());
    hThetaMuMuIS->SetTitle("MC distribution of #theta for intra-slice stitching");
    hThetaMuMuIS->GetXaxis()->SetTitle("#theta [#circ]");
    hThetaMuMuIS->GetXaxis()->SetTitleSize(0.045);
    hThetaMuMuIS->GetYaxis()->SetTitle("Distribution [%]");
    hThetaMuMuIS->GetYaxis()->SetTitleSize(0.045);

    TH1* hThetaMuOtherIS = sThetaMuOtherIS.ToTH1(sThetaMuOtherIS.POT(), kPOT);
    hThetaMuOtherIS->SetLineColor(kBlue);
    hThetaMuOtherIS->Scale(1./hThetaMuOtherIS->Integral());

    TPaveStats* stThetaMuMuIS = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stThetaMuMuIS->AddText(Form("Entries = %.0f", hThetaMuMuIS->GetEntries()));
    stThetaMuMuIS->AddText(Form("Mean = %.3f", hThetaMuMuIS->GetMean()));
    stThetaMuMuIS->AddText(Form("StdDev = %.3f", hThetaMuMuIS->GetStdDev()));
    stThetaMuMuIS->SetTextColor(kRed);
    hThetaMuMuIS->GetListOfFunctions()->Add(stThetaMuMuIS);

    TPaveStats* stThetaMuOtherIS = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stThetaMuOtherIS->AddText(Form("Entries = %.0f", hThetaMuOtherIS->GetEntries()));
    stThetaMuOtherIS->AddText(Form("Mean = %.3f", hThetaMuOtherIS->GetMean()));
    stThetaMuOtherIS->AddText(Form("StdDev = %.3f", hThetaMuOtherIS->GetStdDev()));
    stThetaMuOtherIS->SetTextColor(kBlue);
    hThetaMuOtherIS->GetListOfFunctions()->Add(stThetaMuOtherIS);

    TCanvas* cThetaISCalib = new TCanvas();
    hThetaMuMuIS->Draw("hist");
    hThetaMuOtherIS->Draw("SameHist");
    cThetaISCalib->Write("MC calib: theta intra-slice overlap");

    //==================================================================================

    TH1* hThetaMuMuES = sThetaMuMuES.ToTH1(sThetaMuMuES.POT(), kPOT);
    hThetaMuMuES->SetLineColor(kRed);
    hThetaMuMuES->Scale(1./hThetaMuMuES->Integral());
    hThetaMuMuES->SetTitle("MC distribution of #theta for extra-slice stitching");
    hThetaMuMuES->GetXaxis()->SetTitle("#theta [#circ]");
    hThetaMuMuES->GetXaxis()->SetTitleSize(0.045);
    hThetaMuMuES->GetYaxis()->SetTitle("Distribution [%]");
    hThetaMuMuES->GetYaxis()->SetTitleSize(0.045);

    TH1* hThetaMuOtherES = sThetaMuOtherES.ToTH1(sThetaMuOtherES.POT(), kPOT);
    hThetaMuOtherES->SetLineColor(kBlue);
    hThetaMuOtherES->Scale(1./hThetaMuOtherES->Integral());

    TPaveStats* stThetaMuMuES = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stThetaMuMuES->AddText(Form("Entries = %.0f", hThetaMuMuES->GetEntries()));
    stThetaMuMuES->AddText(Form("Mean = %.3f", hThetaMuMuES->GetMean()));
    stThetaMuMuES->AddText(Form("StdDev = %.3f", hThetaMuMuES->GetStdDev()));
    stThetaMuMuES->SetTextColor(kRed);
    hThetaMuMuES->GetListOfFunctions()->Add(stThetaMuMuES);

    TPaveStats* stThetaMuOtherES = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stThetaMuOtherES->AddText(Form("Entries = %.0f", hThetaMuOtherES->GetEntries()));
    stThetaMuOtherES->AddText(Form("Mean = %.3f", hThetaMuOtherES->GetMean()));
    stThetaMuOtherES->AddText(Form("StdDev = %.3f", hThetaMuOtherES->GetStdDev()));
    stThetaMuOtherES->SetTextColor(kBlue);
    hThetaMuOtherES->GetListOfFunctions()->Add(stThetaMuOtherES);

    TCanvas* cThetaESCalib = new TCanvas();
    hThetaMuMuES->Draw("hist");
    hThetaMuOtherES->Draw("SameHist");
    cThetaESCalib->Write("MC calib: theta extra-slice overlap");

    //==================================================================================

    TH1* hDssMuMuIS = sDssMuMuIS.ToTH1(sDssMuMuIS.POT(), kPOT);
    hDssMuMuIS->SetLineColor(kRed);
    hDssMuMuIS->Scale(1./hDssMuMuIS->Integral());
    hDssMuMuIS->SetTitle("MC Distribution of Dss for intra-slice stitching");
    hDssMuMuIS->GetXaxis()->SetTitle("Dss");
    hDssMuMuIS->GetXaxis()->SetTitleSize(0.045);
    hDssMuMuIS->GetYaxis()->SetTitle("Distribution [%]");
    hDssMuMuIS->GetYaxis()->SetTitleSize(0.045);

    TH1* hDssMuOtherIS = sDssMuOtherIS.ToTH1(sDssMuOtherIS.POT(), kPOT);
    hDssMuOtherIS->SetLineColor(kBlue);
    hDssMuOtherIS->Scale(1./hDssMuOtherIS->Integral());

    TPaveStats* stDssMuMuIS = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stDssMuMuIS->AddText(Form("Entries = %.0f", hDssMuMuIS->GetEntries()));
    stDssMuMuIS->AddText(Form("Mean = %.3f", hDssMuMuIS->GetMean()));
    stDssMuMuIS->AddText(Form("StdDev = %.3f", hDssMuMuIS->GetStdDev()));
    stDssMuMuIS->SetTextColor(kRed);
    hDssMuMuIS->GetListOfFunctions()->Add(stDssMuMuIS);

    TPaveStats* stDssMuOtherIS = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stDssMuOtherIS->AddText(Form("Entries = %.0f", hDssMuOtherIS->GetEntries()));
    stDssMuOtherIS->AddText(Form("Mean = %.3f", hDssMuOtherIS->GetMean()));
    stDssMuOtherIS->AddText(Form("StdDev = %.3f", hDssMuOtherIS->GetStdDev()));
    stDssMuOtherIS->SetTextColor(kBlue);
    hDssMuOtherIS->GetListOfFunctions()->Add(stDssMuOtherIS);

    TCanvas* cDssISCalib = new TCanvas();
    hDssMuMuIS->Draw("hist");
    hDssMuOtherIS->Draw("SameHist");
    cDssISCalib->Write("MC calib: Dss intra-slice overlap");

    //==================================================================================

    TH1* hDssMuMuES = sDssMuMuES.ToTH1(sDssMuMuES.POT(), kPOT);
    hDssMuMuES->SetLineColor(kRed);
    hDssMuMuES->Scale(1./hDssMuMuES->Integral());
    hDssMuMuES->SetTitle("MC distribution of Dss for extra-slice stitching");
    hDssMuMuES->GetXaxis()->SetTitle("Dss");
    hDssMuMuES->GetXaxis()->SetTitleSize(0.045);
    hDssMuMuES->GetYaxis()->SetTitle("Distribution [%]");
    hDssMuMuES->GetYaxis()->SetTitleSize(0.045);

    TH1* hDssMuOtherES = sDssMuOtherES.ToTH1(sDssMuOtherES.POT(), kPOT);
    hDssMuOtherES->SetLineColor(kBlue);
    hDssMuOtherES->Scale(1./hDssMuOtherES->Integral());

    TPaveStats* stDssMuMuES = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stDssMuMuES->AddText(Form("Entries = %.0f", hDssMuMuES->GetEntries()));
    stDssMuMuES->AddText(Form("Mean = %.3f", hDssMuMuES->GetMean()));
    stDssMuMuES->AddText(Form("StdDev = %.3f", hDssMuMuES->GetStdDev()));
    stDssMuMuES->SetTextColor(kRed);
    hDssMuMuES->GetListOfFunctions()->Add(stDssMuMuES);

    TPaveStats* stDssMuOtherES = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stDssMuOtherES->AddText(Form("Entries = %.0f", hDssMuOtherES->GetEntries()));
    stDssMuOtherES->AddText(Form("Mean = %.3f", hDssMuOtherES->GetMean()));
    stDssMuOtherES->AddText(Form("StdDev = %.3f", hDssMuOtherES->GetStdDev()));
    stDssMuOtherES->SetTextColor(kBlue);
    hDssMuOtherES->GetListOfFunctions()->Add(stDssMuOtherES);

    TCanvas* cDssESCalib = new TCanvas();
    hDssMuMuES->Draw("hist");
    hDssMuOtherES->Draw("SameHist");
    cDssESCalib->Write("MC calib: Dss extra-slice overlap");

    //==================================================================================

    TH1* hTssDataIS = sTssDataIS.ToTH1(sTssDataIS.POT(), kPOT);
    hTssDataIS->SetStats(true);
    hTssDataIS->Scale(1./hTssDataIS->Integral());
    hTssDataIS->SetTitle("Distribution of Tss on data for intra-slice stitching");
    hTssDataIS->GetXaxis()->SetTitle("Tss");
    hTssDataIS->GetXaxis()->SetTitleSize(0.045);
    hTssDataIS->GetYaxis()->SetTitle("Distribution [%]");
    hTssDataIS->GetYaxis()->SetTitleSize(0.045);
    hTssDataIS->Write("Data calib: Tss intra-slice");

    TH1* hTssDataES = sTssDataES.ToTH1(sTssDataES.POT(), kPOT);
    hTssDataES->SetStats(true);
    hTssDataES->Scale(1./hTssDataES->Integral());
    hTssDataES->SetTitle("Distribution of Tss on data for extra-slice stitching");
    hTssDataES->GetXaxis()->SetTitle("Tss");
    hTssDataES->GetXaxis()->SetTitleSize(0.045);
    hTssDataES->GetYaxis()->SetTitle("Distribution [%]");
    hTssDataES->GetYaxis()->SetTitleSize(0.045);
    hTssDataES->Write("Data calib: Tss extra-slice");

    TH1* hThetaDataIS = sThetaDataIS.ToTH1(sThetaDataIS.POT(), kPOT);
    hThetaDataIS->SetStats(true);
    hThetaDataIS->Scale(1./hThetaDataIS->Integral());
    hThetaDataIS->SetTitle("Distribution of #theta on data for intra-slice stitching");
    hThetaDataIS->GetXaxis()->SetTitle("#theta [#circ]");
    hThetaDataIS->GetXaxis()->SetTitleSize(0.045);
    hThetaDataIS->GetYaxis()->SetTitle("Distribution [%]");
    hThetaDataIS->GetYaxis()->SetTitleSize(0.045);
    hThetaDataIS->Write("Data calib: theta intra-slice");

    TH1* hThetaDataES = sThetaDataES.ToTH1(sThetaDataES.POT(), kPOT);
    hThetaDataES->SetStats(true);
    hThetaDataES->Scale(1./hThetaDataES->Integral());
    hThetaDataES->SetTitle("Distribution of #theta on data for extra-slice stitching");
    hThetaDataES->GetXaxis()->SetTitle("#theta [#circ]");
    hThetaDataES->GetXaxis()->SetTitleSize(0.045);
    hThetaDataES->GetYaxis()->SetTitle("Distribution [%]");
    hThetaDataES->GetYaxis()->SetTitleSize(0.045);
    hThetaDataES->Write("Data calib: theta extra-slice");

    TH1* hDssDataIS = sDssDataIS.ToTH1(sDssDataIS.POT(), kPOT);
    hDssDataIS->SetStats(true);
    hDssDataIS->Scale(1./hDssDataIS->Integral());
    hDssDataIS->SetTitle("Distribution of Dss on data for intra-slice stitching");
    hDssDataIS->GetXaxis()->SetTitle("Dss");
    hDssDataIS->GetXaxis()->SetTitleSize(0.045);
    hDssDataIS->GetYaxis()->SetTitle("Distribution [%]");
    hDssDataIS->GetYaxis()->SetTitleSize(0.045);
    hDssDataIS->Write("Data calib: Dss intra-slice");

    TH1* hDssDataES = sDssDataES.ToTH1(sDssDataES.POT(), kPOT);
    hDssDataES->SetStats(true);
    hDssDataES->Scale(1./hDssDataES->Integral());
    hDssDataES->SetTitle("Distribution of Dss on data for extra-slice stitching");
    hDssDataES->GetXaxis()->SetTitle("Dss");
    hDssDataES->GetXaxis()->SetTitleSize(0.045);
    hDssDataES->GetYaxis()->SetTitle("Distribution [%]");
    hDssDataES->GetYaxis()->SetTitleSize(0.045);
    hDssDataES->Write("Data calib: Dss extra-slice");

    //===================================================================================================
    //===================================================================================================

    fStitch.Close();
    std::cout << "Spectra saved." << std::endl;
}
