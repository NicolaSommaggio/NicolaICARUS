////////////////////////////////////////////////////////////////////
// Macro to stitch the split tracks of muons
//
// Author: Alessandro Maria Ricci
////////////////////////////////////////////////////////////////////

// C++
#include <string>
#include <sstream>

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
void Stitch() {

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
    //1 file overlays
    //const std::string TargetFiles = "/pnfs/icarus/persistent/users/cfarnese/MC_overlays_2026/group_6664.flat.caf.root";

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

    //ifstream evts("evtNS.txt");
    //std::string linea;
    //while(std::getline(evts,linea))
    //{
    //    std::stringstream ss(linea);
    //    int run, evt, slice_index;
    //    ss >> run >> evt >> slice_index;
    //    std::array<int,3> temp_evt;
    //    temp_evt[0] = run;
    //    temp_evt[1] = evt;
    //    temp_evt[2] = slice_index;
    //    cout << run << " " << evt << " " << slice_index << endl;
    //    NS_evts.push_back(temp_evt);
    //}


    SpectrumLoader NuLoader(TargetFiles);
    std::cout << "Target files loaded." << std::endl;

    // Create the binning schemes for the Vars we wish to plot
    const Binning bMu = Binning::Simple(6, 0, 6);
    const Binning bStitch = Binning::Simple(33, 0, 33);
    const Binning bX = Binning::Simple(200, -400, 400);
    const Binning bY = Binning::Simple(32, -185, 135);
    const Binning bZ = Binning::Simple(125, -1000, 1000);
    const Binning bLen = Binning::Simple(100, 0, 1000);
    const Binning bLen2 = Binning::Simple(40, 0, 400);
    const Binning bE = Binning::Simple(300, 0, 3000);
    const Binning bRatio = Binning::Simple(50, 0, 1);

    //===================================================================================================
    //===================================================================================================

    // Create a Spectrum showing the number of mu, primary daughters of nuMuCC interactions, contained in a single TPC volume
    /*
    Spectrum sMu( "Nmu",                            // A label for the Spectrum
                  bMu,                              // Define the binning
                  NuLoader,                         // Associate this Spectrum with the Loader object (and its target CAF)
                  kMu,                              // The Var to plot
                  kNoSpillCut );                    // The SpillCut to use (none in this case)

    // MC only: create a Spectrum showing how the algorithm stitches
    */
    Spectrum sStitch( "Stitch",
		              bStitch,
		              NuLoader,
		              kStitch,
		              kNoSpillCut );

    // Return the end positions on the x-axis of all tracks before the stitching
    /*
    Spectrum sXendB( "XendB",
                     bX,
		             NuLoader,
                     kXendB,
		             kNoSpillCut );

    // Return the end positions on the x-axis of all tracks after the stitching
    Spectrum sXendA( "XendA",
                      bX,
		              NuLoader,
                      kXendA,
		              kNoSpillCut );

    // Return the end positions on the z-axis of all tracks before the stitching
    Spectrum sZendB( "ZendB",
                      bZ,
		              NuLoader,
                      kZendB,
		              kNoSpillCut );

    // Return the end positions on the z-axis of all tracks after the stitching
    Spectrum sZendA( "ZendA",
                      bZ,
		              NuLoader,
                      kZendA,
		              kNoSpillCut );

    // Return the difference of the muon track length after - before the stitching
    Spectrum sSthLengthAminB( "SthLengthAminB",
                              bLen2,
                              NuLoader,
                              kSthLengthAminB,
                              kNoSpillCut );

    // Return the track length before and after the stitching
    Spectrum sSthLength( "SthLength",
                         NuLoader,
                         bLen,
                         kSthLengthB,
                         bLen,
                         kSthLengthA,
                         kNoSpillCut );

    // Return the track energy before and after the stitching
    Spectrum sSthEnergy( "SthEnergy",
                         NuLoader,
                         bE,
                         kSthEnergyB,
                         bE,
                         kSthEnergyA,
                         kNoSpillCut );

    // MC only: return the end positions on the x-axis of the split nuMuCC tracks before the stitching
    Spectrum sSpTrkXendB( "SpTrkXendB",
                          bX,
		                  NuLoader,
                          kSpTrkXendB,
		                  kNoSpillCut );

    // MC only: return the end positions on the x-axis of the split nuMuCC tracks after the stitching
    Spectrum sSpTrkXendA( "SpTrkXendA",
                          bX,
		                  NuLoader,
                          kSpTrkXendA,
		                  kNoSpillCut );

    // MC only: return the end positions on the z-axis of the split nuMuCC tracks before the stitching
    Spectrum sSpTrkZendB( "SpTrkZendB",
                          bZ,
		                  NuLoader,
                          kSpTrkZendB,
		                  kNoSpillCut );

    // MC only: return the end positions on the z-axis of the split nuMuCC tracks after the stitching
    Spectrum sSpTrkZendA( "SpTrkZendA",
                          bZ,
		                  NuLoader,
                          kSpTrkZendA,
		                  kNoSpillCut );

    // MC only: return |Lreco-Ltrue|/Ltrue before the stitching
    Spectrum sSthLengthRatioB( "SthLengthRatioB",
                               bRatio,
                               NuLoader,
                               kSthLengthRatioB,
                               kNoSpillCut );

    // MC only: return |Lreco-Ltrue|/Ltrue after the stitching
    Spectrum sSthLengthRatioA( "SthLengthRatioA",
                               bRatio,
                               NuLoader,
                               kSthLengthRatioA,
                               kNoSpillCut );

    // MC only: return |Ereco-Etrue|/Etrue before the stitching
    Spectrum sSthEnergyRatioB( "SthEnergyRatioB",
                               bRatio,
                               NuLoader,
                               kSthEnergyRatioB,
                               kNoSpillCut );

    // MC only: return |Ereco-Etrue|/Etrue after the stitching
    Spectrum sSthEnergyRatioA( "SthEnergyRatioA",
                               bRatio,
                               NuLoader,
                               kSthEnergyRatioA ,
                               kNoSpillCut );

    // MC only: return |Ereco-Etrue|/Etrue before and after the stitching
    Spectrum sSthEnergyRatio( "SthEnergyRatio",
                              NuLoader,
                              bRatio,
                              kSthEnergyRatioB,
                              bRatio,
                              kSthEnergyRatioA,
                              kNoSpillCut );

    // MC only: return the PDG of the split tracks
    Spectrum sSpTrkPDG( "SpTrkPDG",
                        bMu,
		                NuLoader,
                        kSpTrkPDG,
		                kNoSpillCut );

    // MC only: return the PDG of the stitched tracks
    Spectrum sSthTrkPDG( "SthTrkPDG",
                         bMu,
		                 NuLoader,
                         kSthTrkPDG,
		                 kNoSpillCut );

    */
    //===================================================================================================
    //===================================================================================================

    // Now that each Spectrum is defined, use the Go() method to populate the Spectrum objects
    NuLoader.Go();
    std::cout << "Spectra created and populated." << std::endl;

    // ---- DRAW -----
    // Write the Spectrum objects to a file as a TH1 object.
    TFile fStitch("Stitch_chi2_ts4.root", "recreate");

    //==================================================================================

    /*
    TH1* hMu = sMu.ToTH1(sMu.POT(), kPOT);
    hMu->SetStats(true);
    hMu->Write("MC only: nuMuCC passing pre-selection");
    */

    TH1* hStitch = sStitch.ToTH1(sStitch.POT(), kPOT);
    hStitch->SetStats(true);
    hStitch->Write("MC only: stitching performance");

    //==================================================================================

    /*
    TH1* hXendB = sXendB.ToTH1(sXendB.POT(), kPOT);
    //auto hXendBInt = hXendB->Integral();
    //hXendB->Scale(1./hXendBInt);
    hXendB->SetMarkerStyle(kFullCircle);
    hXendB->SetTitle("Distribution of the muon x-end points");
    hXendB->GetXaxis()->SetTitle("Muon x-end [cm]");
    hXendB->GetXaxis()->SetTitleSize(0.045);
    hXendB->GetYaxis()->SetTitleSize(0.045);

    TH1* hXendA = sXendA.ToTH1(sXendA.POT(), kPOT);
    //hXendA->Scale(1./hXendBInt);
    hXendA->SetMarkerStyle(kFullCircle);
    hXendA->SetMarkerColor(kRed);
    hXendA->SetLineColor(kRed);

    TPaveStats* stXendB = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stXendB->AddText(Form("Entries = %.0f", hXendB->GetEntries()));
    stXendB->AddText(Form("Mean = %.3f", hXendB->GetMean()));
    stXendB->AddText(Form("StdDev = %.3f", hXendB->GetStdDev()));
    hXendB->GetListOfFunctions()->Add(stXendB);

    TPaveStats* stXendA = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stXendA->AddText(Form("Entries = %.0f", hXendA->GetEntries()));
    stXendA->AddText(Form("Mean = %.3f", hXendA->GetMean()));
    stXendA->AddText(Form("StdDev = %.3f", hXendA->GetStdDev()));
    stXendA->SetTextColor(kRed);
    hXendA->GetListOfFunctions()->Add(stXendA);

    TCanvas* cXend = new TCanvas();
    hXendB->Draw();
    hXendA->Draw("Same");
    cXend->Write("Muon x-end AS vs BS");

    //==================================================================================

    TH1* hZendB = sZendB.ToTH1(sZendB.POT(), kPOT);
    //auto hZendBInt = hZendB->Integral();
    //hZendB->Scale(1./hZendBInt);
    hZendB->SetMarkerStyle(kFullCircle);
    hZendB->SetTitle("Distribution of the muon z-end points");
    hZendB->GetXaxis()->SetTitle("Muon z-end [cm]");
    hZendB->GetXaxis()->SetTitleSize(0.045);
    hZendB->GetYaxis()->SetTitleSize(0.045);

    TH1* hZendA = sZendA.ToTH1(sZendA.POT(), kPOT);
    //hZendA->Scale(1./hZendBInt);
    hZendA->SetMarkerStyle(kFullCircle);
    hZendA->SetMarkerColor(kRed);
    hZendA->SetLineColor(kRed);

    TPaveStats* stZendB = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stZendB->AddText(Form("Entries = %.0f", hZendB->GetEntries()));
    stZendB->AddText(Form("Mean = %.3f", hZendB->GetMean()));
    stZendB->AddText(Form("StdDev = %.3f", hZendB->GetStdDev()));
    hZendB->GetListOfFunctions()->Add(stZendB);

    TPaveStats* stZendA = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stZendA->AddText(Form("Entries = %.0f", hZendA->GetEntries()));
    stZendA->AddText(Form("Mean = %.3f", hZendA->GetMean()));
    stZendA->AddText(Form("StdDev = %.3f", hZendA->GetStdDev()));
    stZendA->SetTextColor(kRed);
    hZendA->GetListOfFunctions()->Add(stZendA);

    TCanvas* cZend = new TCanvas();
    hZendB->Draw();
    hZendA->Draw("Same");
    cZend->Write("Muon z-end AS vs BS");

    //==================================================================================

    TH1* hSthLengthAminB = sSthLengthAminB.ToTH1(sSthLengthAminB.POT(), kPOT);
    hSthLengthAminB->SetStats(true);
    hSthLengthAminB->SetMarkerStyle(kFullCircle);
    hSthLengthAminB->SetTitle("Distribution of the difference of the muon track length AS-BS");
    hSthLengthAminB->GetXaxis()->SetTitle("Muon track length difference AS-BS [cm]");
    hSthLengthAminB->GetXaxis()->SetTitleSize(0.045);
    hSthLengthAminB->GetYaxis()->SetTitleSize(0.045);
    hSthLengthAminB->Write("Muon track length diff. AS-BS");

    TH2* hSthLength = sSthLength.ToTH2(sSthLength.POT(), kPOT);
    hSthLength->SetStats(true);
    hSthLength->SetMarkerStyle(kFullCircle);
    hSthLength->SetTitle("Distribution of the muon track length");
    hSthLength->GetXaxis()->SetTitle("Muon length BS [cm]");
    hSthLength->GetXaxis()->SetTitleSize(0.045);
    hSthLength->GetYaxis()->SetTitle("Muon length AS [cm]");
    hSthLength->GetYaxis()->SetTitleSize(0.045);

    TLine* line = new TLine(0, 0, 1000, 1000);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);

    TCanvas* cSthLength = new TCanvas();
    hSthLength->Draw("SCAT");
    line->Draw("same");
    cSthLength->Write("Muon track length AS vs BS");

    //==================================================================================

    TH2* hSthEnergy = sSthEnergy.ToTH2(sSthEnergy.POT(), kPOT);
    hSthEnergy->SetStats(true);
    hSthEnergy->SetMarkerStyle(kFullCircle);
    hSthEnergy->SetTitle("Distribution of the muon track energy");
    hSthEnergy->GetXaxis()->SetTitle("Muon energy BS [MeV]");
    hSthEnergy->GetXaxis()->SetTitleSize(0.045);
    hSthEnergy->GetYaxis()->SetTitle("Muon energy AS [MeV]");
    hSthEnergy->GetYaxis()->SetTitleSize(0.045);

    TLine* line2 = new TLine(0, 0, 3000, 3000);
    line2->SetLineColor(kRed);
    line2->SetLineWidth(2);

    TCanvas* cSthEnergy = new TCanvas();
    hSthEnergy->Draw("SCAT");
    line2->Draw("same");
    cSthEnergy->Write("Muon track energy AS vs BS");

    //==================================================================================

    TH1* hSpTrkXendB = sSpTrkXendB.ToTH1(sSpTrkXendB.POT(), kPOT);
    //auto hSpTrkXendBInt = hSpTrkXendB->Integral();
    //hSpTrkXendB->Scale(1./hSpTrkXendBInt);
    hSpTrkXendB->SetMarkerStyle(kFullCircle);
    hSpTrkXendB->SetTitle("MC distribution of the split muon x-end points");
    hSpTrkXendB->GetXaxis()->SetTitle("Muon x-end [cm]");
    hSpTrkXendB->GetXaxis()->SetTitleSize(0.045);
    hSpTrkXendB->GetYaxis()->SetTitleSize(0.045);

    TH1* hSpTrkXendA = sSpTrkXendA.ToTH1(sSpTrkXendA.POT(), kPOT);
    //hSpTrkXendA->Scale(1./hSpTrkXendBInt);
    hSpTrkXendA->SetMarkerStyle(kFullCircle);
    hSpTrkXendA->SetMarkerColor(kRed);
    hSpTrkXendA->SetLineColor(kRed);

    TPaveStats* stSpTrkXendB = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stSpTrkXendB->AddText(Form("Entries = %.0f", hSpTrkXendB->GetEntries()));
    stSpTrkXendB->AddText(Form("Mean = %.3f", hSpTrkXendB->GetMean()));
    stSpTrkXendB->AddText(Form("StdDev = %.3f", hSpTrkXendB->GetStdDev()));
    hSpTrkXendB->GetListOfFunctions()->Add(stSpTrkXendB);

    TPaveStats* stSpTrkXendA = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stSpTrkXendA->AddText(Form("Entries = %.0f", hSpTrkXendA->GetEntries()));
    stSpTrkXendA->AddText(Form("Mean = %.3f", hSpTrkXendA->GetMean()));
    stSpTrkXendA->AddText(Form("StdDev = %.3f", hSpTrkXendA->GetStdDev()));
    stSpTrkXendA->SetTextColor(kRed);
    hSpTrkXendA->GetListOfFunctions()->Add(stSpTrkXendA);

    TCanvas* cSpTrkXend = new TCanvas();
    hSpTrkXendB->Draw();
    hSpTrkXendA->Draw("Same");
    cSpTrkXend->Write("MC only: split muon x-end AS vs BS");

    //==================================================================================

    TH1* hSpTrkZendB = sSpTrkZendB.ToTH1(sSpTrkZendB.POT(), kPOT);
    //auto hSpTrkZendBInt = hSpTrkZendB->Integral();
    //hSpTrkZendB->Scale(1./hSpTrkZendBInt);
    hSpTrkZendB->SetMarkerStyle(kFullCircle);
    hSpTrkZendB->SetTitle("MC Distribution of the split muon z-end points");
    hSpTrkZendB->GetXaxis()->SetTitle("Muon z-end [cm]");
    hSpTrkZendB->GetXaxis()->SetTitleSize(0.045);
    hSpTrkZendB->GetYaxis()->SetTitleSize(0.045);

    TH1* hSpTrkZendA = sSpTrkZendA.ToTH1(sSpTrkZendA.POT(), kPOT);
    //hSpTrkZendA->Scale(1./hSpTrkZendBInt);
    hSpTrkZendA->SetMarkerStyle(kFullCircle);
    hSpTrkZendA->SetMarkerColor(kRed);
    hSpTrkZendA->SetLineColor(kRed);

    TPaveStats* stSpTrkZendB = new TPaveStats(0.78, 0.78, 0.93, 0.90, "brNDC");
    stSpTrkZendB->AddText(Form("Entries = %.0f", hSpTrkZendB->GetEntries()));
    stSpTrkZendB->AddText(Form("Mean = %.3f", hSpTrkZendB->GetMean()));
    stSpTrkZendB->AddText(Form("StdDev = %.3f", hSpTrkZendB->GetStdDev()));
    hSpTrkZendB->GetListOfFunctions()->Add(stSpTrkZendB);

    TPaveStats* stSpTrkZendA = new TPaveStats(0.78, 0.65, 0.93, 0.77, "brNDC");
    stSpTrkZendA->AddText(Form("Entries = %.0f", hSpTrkZendA->GetEntries()));
    stSpTrkZendA->AddText(Form("Mean = %.3f", hSpTrkZendA->GetMean()));
    stSpTrkZendA->AddText(Form("StdDev = %.3f", hSpTrkZendA->GetStdDev()));
    stSpTrkZendA->SetTextColor(kRed);
    hSpTrkZendA->GetListOfFunctions()->Add(stSpTrkZendA);

    TCanvas* cSpTrkZend = new TCanvas();
    hSpTrkZendB->Draw();
    hSpTrkZendA->Draw("Same");
    cSpTrkZend->Write("MC only: split muon z-end AS vs BS");

    //==================================================================================
    
    TH1* hSthLengthRatioB = sSthLengthRatioB.ToTH1(sSthLengthRatioB.POT(), kPOT);
    hSthLengthRatioB->SetStats(true);
    hSthLengthRatioB->SetMarkerStyle(kFullCircle);
    hSthLengthRatioB->SetTitle("Distribution of |Lreco-Ltrue|/Ltrue before the stitching");
    hSthLengthRatioB->GetXaxis()->SetTitle("Muon |Lreco-Ltrue|/Ltrue");
    hSthLengthRatioB->GetXaxis()->SetTitleSize(0.045);
    hSthLengthRatioB->GetYaxis()->SetTitleSize(0.045);
    hSthLengthRatioB->Write("MC only: muon track length ratio BS");

    TH1* hSthLengthRatioA = sSthLengthRatioA.ToTH1(sSthLengthRatioA.POT(), kPOT);
    hSthLengthRatioA->SetStats(true);
    hSthLengthRatioA->SetMarkerStyle(kFullCircle);
    hSthLengthRatioA->SetTitle("Distribution of |Lreco-Ltrue|/Ltrue after the stitching");
    hSthLengthRatioA->GetXaxis()->SetTitle("Muon |Lreco-Ltrue|/Ltrue");
    hSthLengthRatioA->GetXaxis()->SetTitleSize(0.045);
    hSthLengthRatioA->GetYaxis()->SetTitleSize(0.045);
    hSthLengthRatioA->Write("MC only: muon track length ratio AS");

    //==================================================================================

    TH1* hSthEnergyRatioB = sSthEnergyRatioB.ToTH1(sSthEnergyRatioB.POT(), kPOT);
    hSthEnergyRatioB->SetStats(true);
    hSthEnergyRatioB->SetMarkerStyle(kFullCircle);
    hSthEnergyRatioB->SetTitle("Distribution of |Ereco-Etrue|/Etrue before the stitching");
    hSthEnergyRatioB->GetXaxis()->SetTitle("Muon |Ereco-Etrue|/Etrue");
    hSthEnergyRatioB->GetXaxis()->SetTitleSize(0.045);
    hSthEnergyRatioB->GetYaxis()->SetTitleSize(0.045);
    hSthEnergyRatioB->Write("MC only: muon track energy ratio BS");

    TH1* hSthEnergyRatioA = sSthEnergyRatioA.ToTH1(sSthEnergyRatioA.POT(), kPOT);
    hSthEnergyRatioA->SetStats(true);
    hSthEnergyRatioA->SetMarkerStyle(kFullCircle);
    hSthEnergyRatioA->SetTitle("Distribution of |Ereco-Etrue|/Etrue after the stitching");
    hSthEnergyRatioA->GetXaxis()->SetTitle("Muon |Ereco-Etrue|/Etrue");
    hSthEnergyRatioA->GetXaxis()->SetTitleSize(0.045);
    hSthEnergyRatioA->GetYaxis()->SetTitleSize(0.045);
    hSthEnergyRatioA->Write("MC only: muon track energy ratio AS");

    TH2* hSthEnergyRatio = sSthEnergyRatio.ToTH2(sSthEnergyRatio.POT(), kPOT);
    hSthEnergyRatio->SetStats(true);
    hSthEnergyRatio->SetMarkerStyle(kFullCircle);
    hSthEnergyRatio->SetTitle("Distribution of |Ereco-Etrue|/Etrue for muon tracks");
    hSthEnergyRatio->GetXaxis()->SetTitle("Ratio before the stitching");
    hSthEnergyRatio->GetXaxis()->SetTitleSize(0.045);
    hSthEnergyRatio->GetYaxis()->SetTitle("Ratio after the stitching");
    hSthEnergyRatio->GetYaxis()->SetTitleSize(0.045);
    hSthEnergyRatio->Write("MC only: muon track energy ratio AS vs BS");

    TCanvas* cSthEnergyRatio = new TCanvas();
    hSthEnergyRatio->Draw("SCAT");
    cSthEnergyRatio->Write("Muon track energy ratio AS vs BS");

    //==================================================================================

    TH1* hSpTrkPDG = sSpTrkPDG.ToTH1(sSpTrkPDG.POT(), kPOT);
    hSpTrkPDG->SetStats(true);
    hSpTrkPDG->SetTitle("MC distribution of the split track PDGs");
    hSpTrkPDG->GetXaxis()->SetTitle("PDGs");
    hSpTrkPDG->GetXaxis()->SetTitleSize(0.045);
    hSpTrkPDG->GetYaxis()->SetTitleSize(0.045);
    hSpTrkPDG->Write("MC only: split track PDGs");

    TH1* hSthTrkPDG = sSthTrkPDG.ToTH1(sSthTrkPDG.POT(), kPOT);
    hSthTrkPDG->SetStats(true);
    hSthTrkPDG->GetXaxis()->SetTitle("PDGs");
    hSthTrkPDG->GetXaxis()->SetTitleSize(0.045);
    hSthTrkPDG->GetYaxis()->SetTitleSize(0.045);
    hSthTrkPDG->Write("MC only: stitched track PDGs");
    */

    fStitch.Close();
    std::cout << "Spectra saved." << std::endl;
}
