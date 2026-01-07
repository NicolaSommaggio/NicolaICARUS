#include <TFile.h>
#include <TDirectory.h>

void fileMaker(){

TFile *GeneralFile = new TFile("GeneralFile.root", "NEW"); 

TDirectory *wire_gap_and_cathode_crossing = GeneralFile->mkdir("wire_gap_and_cathode_crossing");

TDirectory *Eint_MC = GeneralFile->mkdir("Eint_MC");

TDirectory *Eint_MC_not_contained = GeneralFile->mkdir("Eint_MC_not_contained");

TDirectory *Eint_DATA = GeneralFile->mkdir("Eint_DATA");

TDirectory *Eint_DATA_not_contained = GeneralFile->mkdir("Eint_DATA_not_contained");

TDirectory *dedx_MC = GeneralFile->mkdir("dedx_MC");

TDirectory *dedx_MC_with_cosm = GeneralFile->mkdir("dedx_MC_with_cosm");

TDirectory *dedx_DATA = GeneralFile->mkdir("dedx_DATA");

TDirectory *TPCdedx = GeneralFile->mkdir("TPCdedx");

GeneralFile->Close();

}

void CHIfileMaker()
{

    TFile *CHIout = new TFile("CHIout.root", "NEW");

    TDirectory *DATI = CHIout->mkdir("DATI");

    TDirectory *MC = CHIout->mkdir("MC");

    CHIout->Close();

}


void GeneralFileMaker(const char *filename)
{
    TFile *GRAFICI = new TFile(filename, "NEW");
    TDirectory *MUONS = GRAFICI->mkdir("muon");
    TDirectory *PROTONS = GRAFICI->mkdir("proton");

    TDirectory *DATI_P = PROTONS->mkdir("dati");
    TDirectory *MC_P = PROTONS->mkdir("mc");
    TDirectory *DATI_M = MUONS->mkdir("dati");
    TDirectory *MC_M = MUONS->mkdir("mc");
    TDirectory *NC_M = MUONS->mkdir("nc");
    GRAFICI->Close();
}

void FitFileMaker(const char *filename)
{
    TFile *fitResults = new TFile(filename, "NEW");
    TDirectory *mc = (TDirectory*)fitResults->mkdir("mc");
    TDirectory *dati = (TDirectory*)fitResults->mkdir("dati");
    TDirectory *nc = (TDirectory*)fitResults->mkdir("nc");

    TDirectory *mcdEdx = (TDirectory*)mc->mkdir("dEdx");
    TDirectory *mcEint = (TDirectory*)mc->mkdir("Eint");
    TDirectory *mcprotondedx = (TDirectory*)mcdEdx->mkdir("proton");
    TDirectory *mcprotonEint = (TDirectory*)mcEint->mkdir("proton");
    TDirectory *mcmuondedx = (TDirectory*)mcdEdx->mkdir("muon");
    TDirectory *mcmuonEint = (TDirectory*)mcEint->mkdir("muon");

    TDirectory *datidEdx = (TDirectory*)dati->mkdir("dEdx");
    TDirectory *datiEint = (TDirectory*)dati->mkdir("Eint");
    TDirectory *datiprotondedx = (TDirectory*)datidEdx->mkdir("proton");
    TDirectory *datiprotonEint = (TDirectory*)datiEint->mkdir("proton");
    TDirectory *datimuondedx = (TDirectory*)datidEdx->mkdir("muon");
    TDirectory *datimuonEint = (TDirectory*)datiEint->mkdir("muon");

    TDirectory *ncdEdx = (TDirectory*)nc->mkdir("dEdx");
    TDirectory *ncEint = (TDirectory*)nc->mkdir("Eint");

    fitResults->Close();
    delete fitResults;

}

void TPCfile()
{
    TFile *f = new TFile("outTPC.root", "NEW" );
    TDirectory * d = (TDirectory*)f->mkdir("MC");
    TDirectory * d1 = (TDirectory*)f->mkdir("DATI");
}

void likelihood_file()
{
    TFile *f = new TFile("OUTPUT_likelihood.root", "NEW");

    TDirectory *d = (TDirectory*)f->mkdir("muon");
    TDirectory *d1 = (TDirectory*)f->mkdir("proton");
    TDirectory *d2 = (TDirectory*)f->mkdir("muon_hypothesis");
    TDirectory *d3 = (TDirectory*)f->mkdir("proton_hypothesis");

    f->Close();
    delete f;
}


void confronto_file(std::string filename)
{
    TFile *f = new TFile(filename.c_str(), "NEW");

    TDirectory *d = (TDirectory*)f->mkdir("muon");
    TDirectory *d1 = (TDirectory*)f->mkdir("proton");
    TDirectory *d2 = (TDirectory*)d->mkdir("dati");
    TDirectory *d4 = (TDirectory*)d->mkdir("mc");
    TDirectory *d3 = (TDirectory*)d1->mkdir("mc");
    TDirectory *d5 = (TDirectory*)d1->mkdir("dati");

    f->Close();
    delete f;
}


void plotMPV_file(std::string filename)
{
    TFile *f = new TFile(filename.c_str(), "NEW");
    TDirectory *d = (TDirectory*)f->mkdir("muon");
    TDirectory *d1 = (TDirectory*)f->mkdir("proton");
    f->Close();
    delete f;
}


void hReference_file()
{
    TFile *f = new TFile("hReference.root", "NEW");
    f->Close();
    delete f;
}


void out_TPC_2d_file()
{
    TFile * f = new TFile("outTPC_dati_1d_vs_2d.root", "NEW");

    TDirectory *ddati = (TDirectory*)f->mkdir("DATI");
    TDirectory *ddati_10 = (TDirectory*)f->mkdir("DATI_10");
    TDirectory *ddati2d = (TDirectory*)f->mkdir("DATI_2D");
    TDirectory *ddati2d_10 = (TDirectory*)f->mkdir("DATI_2D_10");
    //TDirectory *dmc = (TDirectory*)f->mkdir("MC");
    //TDirectory *dmc2d = (TDirectory*)f->mkdir("MC_2D");

    TDirectory *dmuon = (TDirectory*)f->mkdir("muon");
    TDirectory *dmuon_dati = (TDirectory*)dmuon->mkdir("dati");
    TDirectory *dmuon_dati_2d_10 = (TDirectory*)dmuon->mkdir("dati2d10");
    TDirectory *dmuon_dati_2d = (TDirectory*)dmuon->mkdir("dati2d");

    f->Close();

}

void update_out_TPC_2d_file()
{
    TFile * f = new TFile("outTPC_dati_1d_vs_2d.root", "UPDATE");

    TDirectory *dmuon= (TDirectory*)f->mkdir("proton");
    TDirectory *dmuon_dati = (TDirectory*)dmuon->mkdir("dati");
    TDirectory *dmuon_dati_2d_10 = (TDirectory*)dmuon->mkdir("dati2d10");
    TDirectory *dmuon_dati_2d = (TDirectory*)dmuon->mkdir("dati2d");

    f->Close();

}

void confronto1d2d()
{
    TFile * f = new TFile("Confronto1d2d.root","NEW");
    TDirectory *dmuon = (TDirectory*)f->mkdir("muon");
    TDirectory *dproton = (TDirectory*)f->mkdir("proton");
    TDirectory *dpion = (TDirectory*)f->mkdir("pion");
    
    std::array<TDirectory*,3> directories = {dmuon,dproton,dpion};
    for(auto dir : directories)
    {
        TDirectory * dir1d = (TDirectory*)dir->mkdir("mc1D");
        TDirectory * dir2d = (TDirectory*)dir->mkdir("mc2D");
        TDirectory * dir2d_DNN_PT1 = (TDirectory*)dir->mkdir("mc2D_DNN_PT1");
        TDirectory * dir2d_DNN_PT10 = (TDirectory*)dir->mkdir("mc2D_DNN_PT10");
    }

}

