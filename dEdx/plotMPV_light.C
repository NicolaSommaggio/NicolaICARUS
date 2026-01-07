#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "ReadTree.C"


//funzione che legge dati da file txt (con la possibiltà di leggere solo alcune colonne)
std::vector<std::vector<double>> LeggiColonneDaFile(const std::string& filename, const std::vector<size_t>& colonne_da_leggere,const std::string& rowFilter = "all")
{ 
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
        return {};
    }

    std::vector<std::vector<double>> dati(colonne_da_leggere.size());
    std::string line;
    size_t riga_totale = 0;  // tutte le righe del file

    while (std::getline(file, line)) {
        // Salta righe vuote o commenti
        if (line.empty() || line[0] == '#') {
            ++riga_totale;
            continue;
        }

        // Filtro riga
        if ((rowFilter == "even" && riga_totale % 2 != 0) ||
            (rowFilter == "odd"  && riga_totale % 2 != 1)) {
            ++riga_totale;
            continue;
        }

        std::istringstream iss(line);
        std::vector<double> valori;
        double val;

        while (iss >> val) {
            valori.push_back(val);
        }

        if (*std::max_element(colonne_da_leggere.begin(), colonne_da_leggere.end()) >= valori.size()) {
            std::cerr << "Riga " << riga_totale << " ignorata: troppe poche colonne.\n";
            ++riga_totale;
            continue;
        }

        for (size_t i = 0; i < colonne_da_leggere.size(); ++i) {
            dati[i].push_back(valori[colonne_da_leggere[i]]);
        }

        ++riga_totale;
    }

    return dati;
}


double propagateError(double x, double dx, double y, double dy) {
    double z = (2 * (x - y)) / (x + y);

    // Derivate parziali
    double dz_dx = (4 * y) / std::pow(x + y, 2);
    double dz_dy = (-4 * x) / std::pow(x + y, 2);

    // Propagazione dell'errore
    double dz = std::sqrt(std::pow(dz_dx * dx, 2) + std::pow(dz_dy * dy, 2));

    return dz;
}


std::pair<double,double> mediaPesata(const std::vector<double>& valori, const std::vector<double>& errori) {

    double sommaPesata = 0.0;
    double sommaPesi = 0.0;

    for (size_t i = 0; i < valori.size(); ++i) 
    {
        double peso = 1.0 / (errori[i] * errori[i]);
        sommaPesata += valori[i] * peso;
        sommaPesi += peso;
    }

    double media = sommaPesata / sommaPesi;
    double erroreMedia = std::sqrt(1.0 / sommaPesi);

    std::pair<double,double> mypair;
    mypair.first=media;
    mypair.second=erroreMedia;

    return mypair;
}


void plotMPV(){

gStyle->SetOptStat(0);


/////////////////////////////////// importazione risultati del fit e creazione dei TGraph sovrapposti ///////////////////////////////
std::vector<size_t> colonne = {0, 1, 2};
auto dati = LeggiColonneDaFile("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/dEdx/dump_mpv_1d.txt", colonne, "all");

std::vector<double>& x    = dati[0];
std::vector<double>& y    = dati[1];
std::vector<double>& yerr = dati[2];
std::vector<double> xerr ;

std::vector<size_t> colonne_ref={0,1,2};
auto ref = LeggiColonneDaFile("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/dEdx/dump_mpv_2d.txt", colonne_ref, "all");

std::vector<double>& x_ref    = ref[0];
std::vector<double>& y_ref    = ref[1];
std::vector<double>& yerr_ref = ref[2];
std::vector<double> xerr_ref ;


for(int i=0; i<x.size(); i++)
{
    xerr.push_back(0.15);
    xerr_ref.push_back(0.15);
}

TGraphErrors *g = new TGraphErrors(x.size(), x.data(), y.data(), xerr.data(), yerr.data());

TGraphErrors *g_ref = new TGraphErrors(x_ref.size(), x_ref.data(), y_ref.data(), xerr_ref.data(), yerr_ref.data());


TFile *f= TFile::Open("plotMPV_1d2d.root", "UPDATE");
TDirectory *d = (TDirectory*)f->Get("muon");
d->cd();


g_ref->Write("dati_mpvs_2d_muon",TObject::kOverwrite);
g->Write("dati_mpvs_1d_muon",TObject::kOverwrite);


////////////////////////////////////// differenza DATI MC MPV /////////////////////////////////////////////////////////////
std::vector<double> differenza;
std::vector<double> err_differenza;


TGraphErrors *g_diff = new TGraphErrors();

TGraphErrors *gdiff_dEdx = new TGraphErrors();

std::vector<double> dummy_valori;
std::vector<double> dummy_pesi;
std::vector<double> dummy_valori_dEdx;
std::vector<double> dummy_pesi_dEdx;
std::vector<double> dummy_rr;
std::vector<double> dummy_media;
std::vector<double> dummy_err_media;


for(int i=2; i<x.size(); i++)//si skippano rr=0.75 e rr=1.05
{

    double diff = 2*(y[i]-y_ref[i])/(y[i]+y_ref[i]); 
    double err_diff=propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] );

    int n = g_diff->GetN();
    g_diff->SetPoint(n, x[i], diff );
    g_diff->SetPointError(n,xerr[i],err_diff );

    int n_dEdx = gdiff_dEdx->GetN();
    gdiff_dEdx->SetPoint(n_dEdx, (y[i]+y_ref[i])/2, diff );
    gdiff_dEdx->SetPointError(n_dEdx, 0, err_diff );
}

gdiff_dEdx->Write("gdiff_dEdx"); //anche questo diviso la media

//differenza media in funzione del dEDX
TCanvas *diff_dEdx = new TCanvas("diff_dEdx");
gdiff_dEdx->SetMarkerStyle(7);
gdiff_dEdx->GetXaxis()->SetTitle("(dati+mc)/2 [MeV/cm]");
gdiff_dEdx->GetYaxis()->SetTitle("(dati-mc)/media");
gdiff_dEdx->Draw("AP");
diff_dEdx->Write(0,TObject::kOverwrite);

//differenza media in funzione del rr
TCanvas *diffDatiMC = new TCanvas();
g_diff->Draw("A p");
g_diff->SetMarkerStyle(7);
diffDatiMC->Write("differenza / media",TObject::kOverwrite);

g_diff->Write("differenza_su_media_graph",TObject::kOverwrite);
//differenza/media in funzione del rr rebinnato


f->Close();

}


