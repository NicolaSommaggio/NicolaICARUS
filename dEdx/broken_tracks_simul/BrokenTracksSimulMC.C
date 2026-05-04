#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "helper_stitch_simulz0.h"
#include "helper_broken_xz.h"
#include <TLegend.h>

using namespace ana;

double mediana(std::vector<double> dummy)
{
    std::sort(dummy.begin(), dummy.end());
    double mediana=0;
    int idx=0;
    if(int(dummy.size())%2==1)
    {
        idx=(dummy.size()-1)/2;
        mediana=dummy.at(idx);
    }
    if(dummy.size()%2==0)
    {
        idx =dummy.size()/2-1;
        mediana=( dummy.at(dummy.size()/2) + dummy.at(idx) )/2.;
    }
    return mediana;
}

void BrokenTracksSimulMC()
{

//path MC senza cosmici 
//const std::string fdata = "/storage/gpfs_data/icarus/plain/user/cfarnese/test_genie_v0984_largestat_nuonly/concat_singleneu_*.flat.caf.root";

//1 file MCsenza cosmici
//const std::string fdata= "/storage/gpfs_data/icarus/plain/user/cfarnese/test_genie_v0984_largestat_nuonly/concat_singleneu_2.flat.caf.root";

//path MC + cosmici 
//const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/new_studies_Large_Prod_CNAF/Complete_MC_final/*run*.concatflat.caf.root";

//1 file MC +COSMICI
//const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/new_studies_Large_Prod_CNAF/Complete_MC_final/mc-v09_84_00_01-202403-cnaf-corrsce_run10.concatflat.caf.root";


//YZ Variations 1 file (nu + cosmics)
//const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/New_CAFS_Variations_Nov2024/YZ/CV_yz_FNAL_10new_0.flat.caf.root";

//YZ variations 
const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/New_CAFS_Variations_Nov2024/YZ/CV_yz_FNAL_10new_*.flat.caf.root";

SpectrumLoader loader(fdata);

const Binning kBinz = Binning::Simple(1000,0.,0.);

Spectrum s1("", kBinz, loader, simul_broken_xz ,kCRTPMTNeutrino );

loader.Go();

TH1D* hendx = new TH1D("end x", " ", 100, -400., 400.);
TH1D* hendz = new TH1D("end z", " ", 200, -900, 900);

TH1D *h_median = new TH1D("h_median", "", 200, 0, 20 );

TH2D *dEdx_range = new TH2D("dEdx_range", " ", 500, 0, 25, 300, 0, 30);

cout << "endx size " << endx.size() << " endz size " << endz.size() << " dEdx size " << HITdeMU.size() << endl;

for(int track=0; track<int(endx.size()); track++)
{
    hendx->Fill(endx[track]);
}

for(int track=0; track<int(endz.size()); track++)
{
    hendz->Fill(endz[track]);
}

for(int track=0; track<int(HITrrMU.size()); track++)
{
    std::vector<double> dummy;
    for(int hit=0; hit<int(HITrrMU[track].size())-2; hit++)
    {
        dEdx_range->Fill(HITrrMU[track][hit], HITdeMU[track][hit]);
        if(HITrrMU[track][hit]<5.){dummy.push_back(HITdeMU[track][hit]);}
    }
    if(dummy.size()>0){h_median->Fill(mediana(dummy));}
}


TFile *file = TFile::Open("ConfrontoDatiMC_copy.root", "READ");
TDirectory *d1 = (TDirectory*)file->Get("muon");
TDirectory *d2 = (TDirectory*)d1->Get("dati");
TDirectory *d3 = (TDirectory*)d1->Get("mc");

TH1D *median_dati_originale = (TH1D*)d2->Get("h_median");
TH1D *median_dati = (TH1D*)median_dati_originale->Clone("h_median_dati");

TH1D *median_mc_originale = (TH1D*)d3->Get("h_median");
TH1D *median_mc = (TH1D*)median_mc_originale->Clone("h_median_mc");

TFile *filee = TFile::Open("End.root", "READ");

TH1D *hendx_dati_originale = (TH1D*)filee->Get("end_x_dati");
TH1D *hendx_dati = (TH1D*)hendx_dati_originale->Clone("end_x_dati");

TH1D *hendz_dati_originale = (TH1D*)filee->Get("end_z_dati");
TH1D *hendz_dati = (TH1D*)hendz_dati_originale->Clone("end_z_dati");

TH2D *dEdx_range_mc_originale = (TH2D*)filee->Get("dEdx_range_mc");
TH2D *dEdx_range_mc_noC = (TH2D*)dEdx_range_mc_originale->Clone("dEdx_range_mc_noC");

TFile *fileBroken = TFile::Open("BrokenTracksSimulMC.root", "READ");
TH1D *hmedian_mc_simul_originale = (TH1D*)fileBroken->Get("h_median"); //questi vengono importati già scalati
TH1D* hmedian_mc_simul = (TH1D*)hmedian_mc_simul_originale->Clone("hmedian_mc_simul");


TFile *fileBroken_withCos = TFile::Open("BrokenTracksSimulMC_withCOS.root", "READ");
TH1D *hmedian_simul_cos_originale = (TH1D*)fileBroken_withCos->Get("h_median");
TH1D *hmedian_simul_cos = (TH1D*)hmedian_simul_cos_originale->Clone("hmedian_simul_cos");

TFile *f = new TFile("BrokenTracksSimulMC_YZvariations.root", "RECREATE");
f->cd();

hmedian_simul_cos->Write();
hmedian_mc_simul->Write();

hendx->Write("hendx_not_scaled");
hendx->Scale(1./hendx->Integral());
hendx->Write();
hendz->Write("hendz_not_scaled");
hendz->Scale(1./hendz->Integral());
hendz->Write();

h_median->Scale(1./h_median->Integral());
h_median->Write();
median_dati->Write("median_dati");
median_mc->Write("median mc");

hendx_dati->Scale(1./hendx_dati->Integral());
hendx_dati->Write();
hendz_dati->Scale(1./hendz_dati->Integral());
hendz_dati->Write();
dEdx_range_mc_noC->Write();
dEdx_range->Write();


TCanvas *canvas_endz = new TCanvas();
hendz->Draw("");
hendz->SetLineWidth(2);
hendz_dati->Draw("same");
hendz_dati->SetLineWidth(2);
hendz_dati->SetLineColor(kRed);
TLegend *ll = new TLegend();
ll->AddEntry(hendz, "MC post breaking sim.");
ll->AddEntry(hendz_dati, "DATI");
ll->Draw("same");
canvas_endz->Write("canvas_endz");


TCanvas *c = new TCanvas();
h_median->Draw("hist");
h_median->SetLineWidth(2);
h_median->GetXaxis()->SetTitle("median [MeV/cm]");
h_median->GetYaxis()->SetTitle("counts (area normalized)");
median_dati->Draw("hist same");
median_dati->SetLineColor(kRed);
median_dati->SetLineWidth(2);
hmedian_simul_cos->Draw("hist same");
hmedian_simul_cos->SetLineColor(kMagenta);
hmedian_simul_cos->SetLineWidth(2);

//median_mc->Draw("hist same");
//median_mc->SetLineColor(kOrange);
//median_mc->SetLineWidth(2);
//hmedian_mc_simul->Draw("hist same");
//hmedian_mc_simul->SetLineWidth(2);
//hmedian_mc_simul->SetLineColor(kGreen);
TLegend *l = new TLegend();
l->AddEntry(h_median, "MC + COSMICS YZ variation post break. sim.");
l->AddEntry(median_dati, "DATI");
l->AddEntry(hmedian_simul_cos, "MC + COSMICS post break. sim.");
//l->AddEntry(median_mc, "MC");
//l->AddEntry(hmedian_mc_simul, "MC #nu only post breaking sim.");
l->Draw("same");
c->Write("mediane");
 


}




   



