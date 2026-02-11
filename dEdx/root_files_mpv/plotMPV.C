#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "ReadTree.C"
#include "chiFunction.C"


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
auto dati = LeggiColonneDaFile("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/confronto_muon_dati.txt", colonne, "all");

std::vector<double>& x    = dati[0];
std::vector<double>& y    = dati[1];
std::vector<double>& yerr = dati[2];
std::vector<double> xerr ;

std::vector<size_t> colonne_ref={0,1,2};
auto ref = LeggiColonneDaFile("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/confronto_muon_mc.txt", colonne_ref, "all");

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


//TFile *f= TFile::Open("mpvPlot.root", "UPDATE");
TFile *f= TFile::Open("new_plotMPV.root", "UPDATE");
TDirectory *d = (TDirectory*)f->Get("muon");
d->cd();

TCanvas *canvas = new TCanvas();


g_ref->GetXaxis()->SetTitle("Residual Range [cm]");
g_ref->GetYaxis()->SetTitle("MPV dE/dx [MeV/cm]");
g_ref->Draw("P");
g_ref->SetMarkerStyle(7);
g_ref->SetMarkerColor(kRed);
g_ref->SetLineColor(kRed);
g->Draw("P same");
g->SetMarkerStyle(7);
g->SetMarkerColor(kBlack);
g->SetLineColor(kBlack);
TLegend* l = new TLegend(0.65, 0.75, 0.85, 0.89);
l->AddEntry(g, "DATI");
l->AddEntry(g_ref, "MC");
l->Draw("same");

canvas->Write();

g_ref->Write("mc_mpvs_muon");
g->Write("dati_mpvs_muon");


////////////////////////////////////// differenza DATI MC MPV /////////////////////////////////////////////////////////////
std::vector<double> differenza;
std::vector<double> err_differenza;


TGraphErrors *g_diff = new TGraphErrors();

TGraphErrors *trend = new TGraphErrors();

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
    //condizione sul rr -> prendo in considerazione solo gli ultimi 25 cm
    //if(x[i]>25.5)continue;

    double diff = 2*(y[i]-y_ref[i])/(y[i]+y_ref[i]); //dati - mc
    double err_diff=propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] );

    int n = g_diff->GetN();
    g_diff->SetPoint(n, x[i], diff );
    g_diff->SetPointError(n,xerr[i],err_diff );

    int n_dEdx = gdiff_dEdx->GetN();
    gdiff_dEdx->SetPoint(n_dEdx, (y[i]+y_ref[i])/2, diff );
    gdiff_dEdx->SetPointError(n_dEdx, 0, err_diff );
    

    //RAGGRUPPARE I PUNTI IN CLUSTER
    if(i>=2)
    {
        dummy_valori.push_back(2.*(y[i]-y_ref[i])/(y[i]+y_ref[i]));
        dummy_pesi.push_back(propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] ));
        dummy_valori_dEdx.push_back(2.*(y[i]-y_ref[i])/(y[i]+y_ref[i]));
        dummy_pesi_dEdx.push_back(propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] ));
        dummy_rr.push_back(x[i]);
        dummy_media.push_back((y[i]+y_ref[i])/2.);
        dummy_err_media.push_back(0.5*std::sqrt(yerr[i]*yerr[i] + yerr_ref[i]*yerr_ref[i]));
    }

    
    if(int(dummy_valori.size())==3)
    {
        int nn = trend->GetN();
        std::pair<double,double> media = mediaPesata(dummy_valori, dummy_pesi);
        trend->SetPoint(nn, dummy_rr[1], media.first);
        trend->SetPointError(nn, 0.,  media.second);

        dummy_valori.clear();
        dummy_pesi.clear();
        dummy_rr.clear();

    }
    
}

gdiff_dEdx->Write("gdiff_dEdx");
trend->Write("gtrend");

//differenza media in funzione del dEDX
TCanvas *diff_dEdx = new TCanvas("diff_dEdx");
gdiff_dEdx->SetMarkerStyle(7);
gdiff_dEdx->GetXaxis()->SetTitle("(dati+mc)/2 [MeV/cm]");
gdiff_dEdx->GetYaxis()->SetTitle("(dati-mc)/media");
gdiff_dEdx->Draw("AP");
diff_dEdx->Write();

//differenza media in funzione del rr
TCanvas *diffDatiMC = new TCanvas();
g_diff->Draw("A p");
g_diff->SetMarkerStyle(7);
diffDatiMC->Write("differenza / media");

g_diff->Write("differenza_su_media_graph");
//differenza/media in funzione del rr rebinnato
TCanvas *diffSpaced = new TCanvas();
trend->Draw("A P");
trend->SetMarkerStyle(7);
diffSpaced->Write("trend");


//////////////////////////////////////////// FIT della differenza/media in funzione del RR ////////////////////////////////////////////////

TF1* myFunc = new TF1("myFunc", "[0]*exp(-(x - [1])/[2]) / (1 + exp(-(x - [3])/[4])) + [5]", 1.1, 150);
myFunc->SetParameters(-0.04 , -0.5 , 42 , 3.3 , 1 , 0.035);

TF1* myFunc_25only = new TF1("myFunc_25only", "[0]*exp(-(x - [1])/[2]) / (1 + exp(-(x - [3])/[4])) + [5]", 1., 50.);
myFunc_25only->SetParameters(-0.04 , -0.5 , 42 , 3.3 , 1 , 0.035);



g_diff->Fit(myFunc, "RE");

g_diff->Fit(myFunc_25only, "RE");



//fit differenza/media in funzione del rr per tutto il rr
TCanvas* c = new TCanvas("Fit", "Fit", 800, 600);
trend->SetMarkerStyle(7);
trend->Draw("A P same");
myFunc->SetLineColor(kRed);
myFunc->Draw("same");
c->Write();

myFunc->Write("fit_in_funzione_delRR_fullrange");

//fit differenza/media in funzione del rr solo per gli ultimi 25 cm
TCanvas* cc = new TCanvas("Fit_25_only", "Fit_25_only", 800, 600);
trend->SetMarkerStyle(7);
trend->Draw(" A P same");
myFunc_25only->SetLineColor(kMagenta);
myFunc_25only->Draw("same");
cc->Write();
myFunc_25only->Write("fit_in_funzione_delRR_25only");

TGraph *fitfull = new TGraph();
TGraph *fit25only = new TGraph();
for(int i =0; i<25000; i++ )
{
    int point = fit25only->GetN();
    fit25only->SetPoint(point, i/1000., myFunc_25only->Eval(i/1000.));
}
for(int i =0; i<150000; i++ )
{
    int point = fitfull->GetN();
    fitfull->SetPoint(point,i/1000.,myFunc->Eval(i/1000.) );
}
fitfull->Write("fitfull");
fit25only->Write("fit25only");

f->Close();

}


void plotMPVpro(){

gStyle->SetOptStat(0);


/////////////////////////////////// importazione risultati del fit e creazione dei TGraph sovrapposti ///////////////////////////////
std::vector<size_t> colonne = {0, 1, 2};
auto dati = LeggiColonneDaFile("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/confronto_proton_dati.txt", colonne, "all");

std::vector<double>& x    = dati[0];
std::vector<double>& y    = dati[1];
std::vector<double>& yerr = dati[2];
std::vector<double> xerr ;

std::vector<size_t> colonne_ref={0,1,2};
auto ref = LeggiColonneDaFile("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/confronto_proton_mc.txt", colonne_ref, "all");

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


//TFile *f= TFile::Open("mpvPlot.root", "UPDATE");
TFile *f= TFile::Open("new_plotMPV.root", "UPDATE");
TDirectory *d = (TDirectory*)f->Get("proton");
d->cd();

TCanvas *canvas = new TCanvas();


g_ref->GetXaxis()->SetTitle("Residual Range [cm]");
g_ref->GetYaxis()->SetTitle("MPV dE/dx [MeV/cm]");
g_ref->Draw("P");
g_ref->SetMarkerStyle(7);
g_ref->SetMarkerColor(kRed);
g_ref->SetLineColor(kRed);
g->Draw("P same");
g->SetMarkerStyle(7);
g->SetMarkerColor(kBlack);
g->SetLineColor(kBlack);
TLegend* l = new TLegend(0.65, 0.75, 0.85, 0.89);
l->AddEntry(g, "DATI");
l->AddEntry(g_ref, "MC");
l->Draw("same");

canvas->Write();

g_ref->Write("mc_mpvs_proton");
g->Write("dati_mpvs_proton");


////////////////////////////////////// differenza DATI MC MPV /////////////////////////////////////////////////////////////
std::vector<double> differenza;
std::vector<double> err_differenza;


TGraphErrors *g_diff = new TGraphErrors();

TGraphErrors *trend = new TGraphErrors();

TGraphErrors *gdiff_dEdx = new TGraphErrors();


std::vector<double> dummy_valori;
std::vector<double> dummy_pesi;
std::vector<double> dummy_rr;



for(int i=2; i<x.size(); i++)//si skippano rr=0.75 e rr=1.05 (solo per la differenza, le reference iniziano a rr_centroide = 0.75 cm)
{
    //condizione sul rr -> prendo in considerazione solo gli ultimi 25 cm
    //if(x[i]>25.5)continue;

    double diff = 2*(y[i]-y_ref[i])/(y[i]+y_ref[i]); //dati - mc
    double err_diff=propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] );

    int n = g_diff->GetN();
    g_diff->SetPoint(n, x[i], diff );
    g_diff->SetPointError(n,xerr[i],err_diff );

    int n_dEdx = gdiff_dEdx->GetN();
    gdiff_dEdx->SetPoint(n_dEdx, (y[i]+y_ref[i])/2, diff );
    gdiff_dEdx->SetPointError(n_dEdx, 0, err_diff );
    

    //RAGGRUPPARE PUNTI A CLUSTER DI 3
    if(i>=2)
    {
        dummy_valori.push_back(2.*(y[i]-y_ref[i])/(y[i]+y_ref[i]));
        dummy_pesi.push_back(propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] ));
        dummy_rr.push_back(x[i]);
    }

    if(int(dummy_valori.size())==3)
    {
        int nn = trend->GetN();
        std::pair<double,double> media = mediaPesata(dummy_valori, dummy_pesi);
        trend->SetPoint(nn, dummy_rr[1], media.first);
        trend->SetPointError(nn, 0.,  media.second);
       
        dummy_valori.clear();
        dummy_pesi.clear();
        dummy_rr.clear();
    }
    
}

gdiff_dEdx->Write("gdiff_dEdx");
trend->Write("gtrend");

//differenza/media in funzione del dEdx
TCanvas *diff_dEdx = new TCanvas("diff_dEdx_pro");
gdiff_dEdx->SetMarkerStyle(7);
gdiff_dEdx->GetXaxis()->SetTitle("(dati+mc)/2 [MeV/cm]");
gdiff_dEdx->GetYaxis()->SetTitle("(dati-mc)/media");
gdiff_dEdx->Draw("AP");
diff_dEdx->Write();

//differenza/media in funzione del rr 
TCanvas *diffDatiMC = new TCanvas();
g_diff->Draw("A p");
g_diff->SetMarkerStyle(7);
diffDatiMC->Write("differenza / media");

g_diff->Write("differenza_su_media_graph");

// differenza/media in funzione del RR rebinnata
TCanvas *diffSpaced = new TCanvas();
trend->Draw("A P");
trend->SetMarkerStyle(7);
diffSpaced->Write("trend");


f->Close();

}

// funzione che calcola le mediane e l'Eint fino a 12.5 cm, con possibilià di correggere il dEdx MC
std::vector<TH1D*> confrontoEintDatiMC(std::string dati_o_mc, std::string particle, TFile *file, TF1* correction_function, bool correcting, double median_threshold)
{
    
    EventsData dat = load_data(dati_o_mc, particle, "Np", "full");
    cout << dat.tree->GetEntries() << " MC total muon tracks" << endl;

    std::string corr;
    if(correcting==true){corr="corrected";}
    else{corr="";}

    TH1D *Eint = new TH1D(Form("Eint_%s_%s_%s", particle.c_str(), dati_o_mc.c_str(), corr.c_str() ), "", 120., 0., 120. );
    TH1D *Mediana_Eint = new TH1D(Form("Mediana_Eint_%s_%s_%s", particle.c_str(), dati_o_mc.c_str(), corr.c_str() ), "", 120, 0., 120. );
    TH1D *hmediana = new TH1D(Form("hmediana_%s_%s_%s", particle.c_str(), dati_o_mc.c_str(), corr.c_str()), "", 200, 0, 20);

    ofstream controllo(Form("controllo_%s_%s.txt", particle.c_str(), dati_o_mc.c_str()));
    ofstream controlloEint(Form("controllo_Eint_%s_%s.txt", particle.c_str(), dati_o_mc.c_str()));

    std::vector<double> vec_rr;
    std::vector<double> vec_Eint;
    std::vector<double> vec_dE;
    std::vector<double> vec_pitch;

  double factor=1;
  int counter=0;
  int counter2=0;


    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);
        std::vector<double> trackdEdx;

        //calcolo della mediana
        for(int hit=0; hit<int(dat.track.rr->size())-2; hit++) //removing the last 2 hits (lower rr) to calculate the median
        { 
            if(dat.track.rr->at(hit)<5.)
            {
                if(correcting==true){factor = (2+correction_function->Eval(dat.track.dE->at(hit)))/(2-correction_function->Eval(dat.track.dE->at(hit)));}
                else{factor=1.;}

                trackdEdx.push_back(dat.track.dE->at(hit)*factor);
            } 
            
            vec_rr.push_back(dat.track.rr->at(hit));
            vec_pitch.push_back(dat.track.pitch->at(hit));

            if(correcting==true){factor = (2+correction_function->Eval(dat.track.dE->at(hit)))/(2-correction_function->Eval(dat.track.dE->at(hit)));}
            else{factor=1.;}
            vec_dE.push_back(dat.track.dE->at(hit)*factor);
            
            if(hit==int(dat.track.rr->size())-2)
              {
                vec_rr.push_back(dat.track.rr->at(hit+1));
                vec_rr.push_back(dat.track.rr->at(hit+2));
                vec_pitch.push_back(dat.track.pitch->at(hit+1));
                vec_pitch.push_back(dat.track.pitch->at(hit+2));

                if(correcting==true){factor = (2+correction_function->Eval(dat.track.dE->at(hit)))/(2-correction_function->Eval(dat.track.dE->at(hit)));}
                else{factor=1.;}
                vec_dE.push_back(dat.track.dE->at(hit+1)*factor);
                if(correcting==true){factor = (2+correction_function->Eval(dat.track.dE->at(hit)))/(2-correction_function->Eval(dat.track.dE->at(hit)));}
                else{factor=1.;}
                vec_dE.push_back(dat.track.dE->at(hit+2)*factor);
              }
              
        }
        if(int(trackdEdx.size())==0)
        {
            counter2+=1;

            if(dati_o_mc=="mc")
            {
                for(int i=0; i<3; i++){controllo << dat.vertex.reco_vertex->at(i) << " ";}
                for(int i=0; i<3; i++){controllo << dat.vertex.true_vertex->at(i) << " ";}
                for(int i=0; i<3; i++){controllo << dat.track.start_reco->at(i) << " ";}
                for(int i=0; i<3; i++){controllo << dat.track.start_true->at(i) << " ";}
                for(int i=0; i<3; i++){controllo << dat.track.end_reco->at(i) << " ";}
                for(int i=0; i<3; i++){controllo << dat.track.end_true->at(i) << " ";}
                controllo << dat.track.len_reco << " ";
                controllo << dat.track.len_true << endl;
                
            }
            else if(dati_o_mc=="dati")
            {
                for(int i=0; i<3; i++){controllo << dat.vertex.reco_vertex->at(i) << " ";}
                for(int i=0; i<3; i++){controllo << dat.track.start_reco->at(i) << " ";} 
                for(int i=0; i<3; i++){controllo << dat.track.end_reco->at(i) << " ";} 
                controllo << dat.track.len_reco << endl;
            }
             
        }
        
        std::vector<double> vec_rr_inverted=vec_rr;
        std::vector<double> vec_dE_inverted=vec_dE;
        std::vector<double> vec_pitch_inverted=vec_pitch;
        std::reverse(vec_rr_inverted.begin(), vec_rr_inverted.end());
        std::reverse(vec_dE_inverted.begin(), vec_dE_inverted.end());
        std::reverse(vec_pitch_inverted.begin(), vec_pitch_inverted.end());
        
        //calcolo Eint
        double Eint_val=0;
        for(int i=0; i<int(vec_rr_inverted.size()); i++)
        {
          Eint_val= Eint_val + vec_dE_inverted[i]*vec_pitch_inverted[i];
          vec_Eint.push_back(Eint_val);
        }
        
        
        //Eint con tutte le hit e senza condizione sulla mediana 
        int found_hit=-10;
        double best_distance=10e6;
        for(int hit=0; hit<int(vec_rr_inverted.size()); hit++)
            {
              if(vec_rr_inverted[hit]<12.5)
                {
                    double distance=std::abs(vec_rr_inverted[hit]-12.5);
                    if(distance<=best_distance){best_distance=distance; found_hit=hit;}
                }   
            }
            if(found_hit > 0 && found_hit < int(vec_Eint.size())-1)
            {
                double nextEint = (vec_dE_inverted[found_hit+1]*vec_pitch_inverted[found_hit+1])/(vec_rr_inverted[found_hit+1]-vec_rr_inverted[found_hit])*best_distance;
                Eint->Fill(vec_Eint[found_hit]+nextEint); 
                controlloEint << found_hit << " " << vec_rr_inverted[found_hit] << " " << vec_rr_inverted[found_hit+1] << " "; 
                for(int idx=0; idx<found_hit; idx++)
                {
                    controlloEint << vec_Eint[idx] << " ";
                }
                controlloEint << nextEint << endl;
            }

        double med =0;
        if(trackdEdx.size()>0) 
        {
            med= mediana(trackdEdx);
        }
        hmediana->Fill(med);
        

        if(med>median_threshold)
        {
            counter+=1;
            int found_hit=-10;
            double best_distance=10e6;
            for(int hit=0; hit<int(vec_rr_inverted.size()); hit++)
            {
                if(vec_rr_inverted[hit]<12.5)
                {
                    double distance=std::abs(vec_rr_inverted[hit]-12.5);
                    if(distance<=best_distance){best_distance=distance; found_hit=hit;}
                } 
            }
            if(found_hit > 0 && found_hit < int(vec_Eint.size()))
            {
                double nextEint = (vec_dE_inverted[found_hit+1]*vec_pitch_inverted[found_hit+1])/(vec_rr_inverted[found_hit+1]-vec_rr_inverted[found_hit])*best_distance;
                Mediana_Eint->Fill(vec_Eint[found_hit]+nextEint); 
            }
        }
        
       
       vec_rr.clear();
       vec_dE.clear();
       vec_pitch.clear();
       vec_Eint.clear();
       vec_rr_inverted.clear();
       vec_dE_inverted.clear();
       vec_pitch_inverted.clear(); 
        
        
        
    }

    cout << counter << " tracks with median > " << median_threshold << endl;
    cout << counter2 << " tracks with no hit with rr < 5 cm" <<endl;

    /*
    TDirectory *d = (TDirectory*)file->Get(particle.c_str());
    if(correcting==true){file->cd();}
    else{d->cd();}
    Eint->Scale(1./Eint->Integral());
    Eint->Write(0,TObject::kOverwrite);
    Mediana_Eint->Scale(1./Mediana_Eint->Integral());
    Mediana_Eint->Write(0, TObject::kOverwrite);
    hmediana->Scale(1./hmediana->Integral());
    hmediana->Write(0, TObject::kOverwrite);
    */

    std::vector<TH1D*> retunredHisto;
    retunredHisto.push_back(Eint);
    retunredHisto.push_back(Mediana_Eint);
    retunredHisto.push_back(hmediana);

    return retunredHisto;
    
}

Double_t correction_function(Double_t *x, Double_t *par)
{
    if(x[0]<=1.64000) return 0.0354575;
    else if(x[0]>=15.0000) return -0.0164076;
    else return exp(par[0]+par[1]*x[0])+par[2]+par[3]*exp(-0.5*pow((x[0]-par[4])/par[5],2));
}

Double_t correction_function_25only(Double_t *x, Double_t *par)
{
    if(x[0]<=2) return 0.0127159;
    else if (x[0]>=15) return -0.0235248;
    else return exp(par[0]+par[1]*x[0])+par[2]+par[3]*exp(-0.5*pow((x[0]-par[4])/par[5],2));
}

Double_t correction_muon(Double_t *x, Double_t *par)
{
    if(x[0]<=2.) return 0.0120827;
    else if(x[0]>=5.) return 0.0373432;
    else return par[0]*(x[0]-par[1])*(x[0]-par[1])+par[2]*x[0]+par[3];
}

Double_t correction_proton(Double_t *x, Double_t *par)
{
    if(x[0]<=4.2) return 0.0125211;
    else if(x[0]>=15.) return  -0.0288104;
    else return par[0]+par[1]*TMath::Log(1+par[2]*(x[0]-par[3])*(x[0]-par[3]));
}



void dEdxCorrection()
{
     gROOT->ProcessLine(".L libdaughtersInfo.so");
    //importo TGraph MPV 
    //TFile *f = TFile::Open("mpvPlot.root","UPDATE");
    TFile *f= TFile::Open("new_plotMPV_25only.root", "UPDATE");
    TDirectory *d = (TDirectory*)f->Get("muon");
    auto *gdiff_dEdx_muon = (TGraphErrors*)d->Get("gdiff_dEdx");
    auto *g_dati_mpvs_muon = (TGraphErrors*)d->Get("dati_mpvs_muon");
    auto *g_mc_mpvs_muon = (TGraphErrors*)d->Get("mc_mpvs_muon");
    TDirectory *d1 = (TDirectory*)f->Get("proton");
    auto *gdiff_dEdx_proton = (TGraphErrors*)d1->Get("gdiff_dEdx");
    auto *g_dati_mpvs_proton = (TGraphErrors*)d1->Get("dati_mpvs_proton");
    auto *g_mc_mpvs_proton = (TGraphErrors*)d1->Get("mc_mpvs_proton");


    f->cd();

    std::vector<double> x_val;
    std::vector<double> y_val;
    std::vector<double> x_err;
    std::vector<double> y_err;
    
    // dEdx MPVs difference in funcion of MC,DATA mean 
    for(int i=0; i<gdiff_dEdx_muon->GetN(); i++)
    {
        double x,y,ex,ey;
        gdiff_dEdx_muon->GetPoint(i,x,y);
        ex = gdiff_dEdx_muon->GetErrorX(i);
        ey = gdiff_dEdx_muon->GetErrorY(i);
        x_val.push_back(x);
        y_val.push_back(y);
        x_err.push_back(ex);
        y_err.push_back(ey);
    }

    for(int i=0; i<gdiff_dEdx_proton->GetN(); i++)
    {
        double x,y,ex,ey;
        gdiff_dEdx_proton->GetPoint(i,x,y);
        ex = gdiff_dEdx_proton->GetErrorX(i);
        ey = gdiff_dEdx_proton->GetErrorY(i);
        x_val.push_back(x);
        y_val.push_back(y);
        x_err.push_back(ex);
        y_err.push_back(ey);

    }
gdiff_dEdx_muon->Write("gdiff_dEdx_muon", TObject::kOverwrite);
gdiff_dEdx_proton->Write("gdiff_dEdx_proton", TObject::kOverwrite);

////////////////////// CORRECTION FUNCTIONS - RIASSUNTO FINALE DEI PARAMETRI FITTATI ////////////////////////// 
/*
    //exp([0]+[1]*x)+[2]+[3]*exp(-0.5*((x-[4])/[5])**2)

    //full range
    //TF1 * fitMCcorr = new TF1("fitMCcorr", "expo(0)+pol0(2)+gaus(3)",0.,14.39137);
    //fitMCcorr->SetParameters(4.42099e-01, -2.19646e+00, -1.30966e-02, 4.39112e-02, 7.16513e+00, -2.79178e+00);
    //estremo sinistro: 0.0354575 (1.64)
    //estremo destro: -0.0164076 (15)


    //solo 25 cm
    //fitMCcorr->SetParameters(-1.08107e+00, -1.38302e+00, -2.84635e-02, 5.87297e-02, 7.17814e+00, -3.51461e+00);
    //estremo sinistro: 0.0127159 (2)
    //estremo destro: -0.0235248 (15)

    //solo muoni
    //fitMCcorr->SetParameters(7.65506e-03, 1.59985e+00, -2.06713e-02, 5.21996e-02);
    //estremo sinistro: 0.0120827(2)
    //estremo destro: 0.0373432 (5)

    //solo protoni
    //fitMCcorr->SetParameters(3.04412e-02, -3.33271e-02, 8.03482e-02, 7.17696e+00);
    //estremo sinistro: 0.0125211 (4.2)
    //estremo destro: -0.0288104  (15.)
*/    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////// fit correzione dEdx MC - here we perform the fit ///////////////////////////////////////////////////////////////////

    TGraphErrors *g = new TGraphErrors(x_val.size(), x_val.data(), y_val.data(), x_err.data(), y_err.data());
    g->SetLineWidth(2);
    g->SetMarkerStyle(7);
    g->SetMarkerColor(kBlue);
    g->SetLineColor(kAzure-4);

    //solo muoni
    //TF1 *fitfuncMU_fit = new TF1("fitfun", "[0]*(x-[1])*(x-[1])+[2]*x+[3] ", 2.,5.);
    //fitfuncMU_fit->SetParameters(0.006, 1.9, -0.02, 0.056);
    //gdiff_dEdx_muon->Fit(fitfuncMU_fit, "RE");

    //solo protoni
    TF1 *fitfuncPRO_fit = new TF1("fitfunPRO", "[0]+[1]*TMath::Log(1+[2]*(x-[3])*(x-[3]))", 4.2, 15.);
    fitfuncPRO_fit->SetParameters(0.03,-0.04,0.05,6.8);
    gdiff_dEdx_proton->Fit(fitfuncPRO_fit, "RE");


    //fit
    
    //TF1 *fitExpPolGauss = new TF1("fitExpPolGauss", "expo(0)+pol0(2)+gaus(3)",0.4821658,15.);
    
    //fitExpPolGauss->SetParameters(0.548, -2.236, -0.01751, 0.04862, 7.14, -3.001 );
    //g->Fit(fitExpPolGauss, "RE");
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //this are the already fitted functions
    //totale
    TF1 *fitExpPolGauss_tot = new TF1("fitExpPolGauss", correction_function, 0., 16., 6);
    fitExpPolGauss_tot->SetParameters(3.13069e-01, -2.08971e+00, -1.81601e-02, 4.87814e-02, 7.17860e+00, -3.03200e+00); //full range
    

    //25only;
    
    TF1 *fitExpPolGauss_25only = new TF1("fitExpPolGauss", correction_function_25only, 0., 16., 6);
    fitExpPolGauss_25only->SetParameters(-1.08107e+00, -1.38302e+00, -2.84635e-02, 5.87297e-02, 7.17814e+00, -3.51461e+00); //only 25 centimeters

    TF1* fitfuncMU = new TF1("fitfuncMU", correction_muon, 0., 16. , 4);
    fitfuncMU->SetParameters(7.65506e-03, 1.59985e+00, -2.06713e-02, 5.21996e-02);

    TF1 *fitfuncPRO = new TF1("fitfuncPRO", correction_proton, 0., 16., 4);
    fitfuncPRO->SetParameters(3.04412e-02, -3.33271e-02, 8.03482e-02, 7.17696e+00);
    

    TGraph *function_tot = new TGraph();
    TGraph *function_25only = new TGraph();
    TGraph *function_mu = new TGraph();
    TGraph *function_pro = new TGraph();
    for(int i =0; i<16000; i++ )
    {
        int point = function_tot->GetN();
        function_tot->SetPoint(point,i/1000.,fitExpPolGauss_tot->Eval(i/1000.) );
        function_mu->SetPoint(point, i/1000., fitfuncMU->Eval(i/1000.));
        function_pro->SetPoint(point, i/1000., fitfuncPRO->Eval(i/1000.));
        function_25only->SetPoint(point, i/1000., fitExpPolGauss_25only->Eval(i/1000.) );
    }

    function_tot->Write("totale",TObject::kOverwrite);
    function_25only->Write("25only",TObject::kOverwrite);
    function_mu->Write("muoni_specific",TObject::kOverwrite);
    function_pro->Write("protoni_specific",TObject::kOverwrite);

    function_tot->SetMarkerColor(kRed);
    function_tot->SetMarkerStyle(7);
    function_tot->SetLineColor(kRed);
    function_mu->SetMarkerColor(kGreen);
    function_mu->SetMarkerStyle(7);
    function_mu->SetLineColor(kGreen);
    function_pro->SetMarkerColor(kBlack);
    function_pro->SetMarkerStyle(7);
    function_pro->SetLineColor(kBlack);

    g->Write("g_diff_DATI_MC_dEdx",TObject::kOverwrite);

    TCanvas *c = new TCanvas();
    g->Draw("A P");
    g->GetXaxis()->SetTitle("(dati+mc)/2 [MeV/cm]");
    g->GetYaxis()->SetTitle("(dati-mc)/media");
    fitExpPolGauss_tot->SetLineColor(kRed);
    fitExpPolGauss_tot->Draw("same");
    c->Write(0,TObject::kOverwrite);


    TCanvas *c_Fit_colorata = new TCanvas("c_Fit_colorata");
    gdiff_dEdx_proton->GetXaxis()->SetLimits(0.,16.);
    gdiff_dEdx_proton->Draw("AP");
    gdiff_dEdx_proton->SetLineColor(kOrange);
    gdiff_dEdx_proton->SetMarkerColor(kOrange-3);
    gdiff_dEdx_proton->SetMarkerStyle(7);
    gdiff_dEdx_muon->Draw("P same");
    gdiff_dEdx_muon->SetLineColor(kAzure-4);
    gdiff_dEdx_muon->SetMarkerStyle(7);
    gdiff_dEdx_muon->SetMarkerColor(kAzure);
    function_tot->Draw("p same");
    function_mu->Draw("p same");
    function_pro->Draw("p same");    
    TLegend *l2 = new TLegend();
    l2->AddEntry(gdiff_dEdx_muon, "muons");
    l2->AddEntry(gdiff_dEdx_proton, "protons");
    l2->AddEntry(function_tot, "global fit");
    l2->AddEntry(function_mu, "muon fit");
    l2->AddEntry(function_pro, "proton fit");
    //l2->AddEntry(fitExpPolGauss, "fit");
    l2->Draw("same");
    c_Fit_colorata->Write(0,TObject::kOverwrite);
    
////////////////////////////// Calcolo dei Residui separati per mu e pro //////////////////////////////////////////////////////////////////////////

    int nm = gdiff_dEdx_muon->GetN();
    std::vector<double> x_vals_muon(nm), res_vals_muon(nm), ex_vals_muon(nm), ey_vals_muon(nm);

    int np = gdiff_dEdx_proton->GetN();
    std::vector<double> x_vals_proton(np), res_vals_proton(np), ex_vals_proton(np), ey_vals_proton(np);

for (int i = 0; i < nm; ++i) {
        double x, y;
        gdiff_dEdx_muon->GetPoint(i, x, y);
        double yfit = fitExpPolGauss_tot->Eval(x);
        double yerr = gdiff_dEdx_muon->GetErrorY(i);
        double residual = y - yfit;

        x_vals_muon[i] = x;
        res_vals_muon[i] = residual;
        ex_vals_muon[i] = 0.;
        ey_vals_muon[i] = yerr;

    }

    for (int i = 0; i < np; ++i) {
        double x, y;
        gdiff_dEdx_proton->GetPoint(i, x, y);
        double yfit = fitExpPolGauss_tot->Eval(x);
        double yerr = gdiff_dEdx_proton->GetErrorY(i);
        double residual = y - yfit;

        x_vals_proton[i] = x;
        res_vals_proton[i] = residual;
        ex_vals_proton[i] = 0.;
        ey_vals_proton[i] = yerr;

    }

    // Residuals TGraphErrors
    TGraphErrors* gResiduals_muon = new TGraphErrors(nm, &x_vals_muon[0], &res_vals_muon[0], &ex_vals_muon[0], &ey_vals_muon[0]);
    TGraphErrors* gResiduals_proton = new TGraphErrors(np, &x_vals_proton[0], &res_vals_proton[0], &ex_vals_proton[0], &ey_vals_proton[0]);

    gResiduals_muon->Write("residui_muone_totale",TObject::kOverwrite);
    gResiduals_proton->Write("residui_protone_totale",TObject::kOverwrite);

    TCanvas *c_residual = new TCanvas();

    gResiduals_muon->Draw("AP");
    gResiduals_muon->GetXaxis()->SetLimits(0.,15.);
    gResiduals_muon->SetMarkerStyle(7);
    gResiduals_muon->SetLineColor(kAzure-4);
    gResiduals_muon->SetMarkerColor(kAzure);
    gResiduals_proton->SetMarkerStyle(7);
    gResiduals_proton->SetLineColor(kOrange);
    gResiduals_proton->SetMarkerColor(kOrange-3);
    gResiduals_proton->Draw("P same");
    TLegend *l = new TLegend();
    l->AddEntry(gResiduals_muon, "muons");
    l->AddEntry(gResiduals_proton, "protons");
    l->Draw("same");

    // Zero line for visual reference
    TLine* zero = new TLine(0, 0, 15, 0);
    zero->SetLineStyle(2);
    zero->SetLineColor(kGray+2);
    zero->Draw("same");

    c_residual->Write(0,TObject::kOverwrite);

 
/////////  Ricalcolo Eint con la correzzione del dEdx sul MC /////////////////////////////////////////////////////////////////////
   
    bool ricalcolo_Eint = true;
    if(ricalcolo_Eint){
    //Eint, Eint_condizione_sulla_mediana, mediana
    std::vector<TH1D*> histovec_dati_mu;
    std::vector<TH1D*> histovec_mc_mu_pre;
    std::vector<TH1D*> histovec_mc_mu_post;

    histovec_mc_mu_post=confrontoEintDatiMC("mc", "muon", f, fitfuncMU, true, 3.2);
    histovec_mc_mu_pre=confrontoEintDatiMC("mc", "muon", f, fitfuncMU, false, 3.2);
    
    
    histovec_dati_mu=confrontoEintDatiMC("dati", "muon", f, fitfuncMU, false, 3.2);
    
    f->cd();

    TCanvas *canvas1 = new TCanvas("Confronto_Eint_dati_mc_mu", " ");
    histovec_dati_mu[0]->Draw("P");
    histovec_dati_mu[0]->SetLineWidth(2);
    histovec_mc_mu_pre[0]->Draw("hist same");
    histovec_mc_mu_pre[0]->SetLineWidth(2);
    histovec_mc_mu_pre[0]->SetLineColor(kRed);
    histovec_mc_mu_pre[0]->SetLineStyle(kDashed);
    histovec_mc_mu_post[0]->Draw("hist same");
    histovec_mc_mu_post[0]->SetLineColor(kGreen);
    histovec_mc_mu_post[0]->SetLineWidth(2);
    TLegend *l1 = new TLegend();
    l1->AddEntry(histovec_dati_mu[0], "MUON DATA");
    l1->AddEntry(histovec_mc_mu_pre[0], "MUON MC"); 
    l1->AddEntry(histovec_mc_mu_post[0], "MUON MC CORRECTED");
    l1->Draw("same");
    canvas1->Write(0, TObject::kOverwrite);

    histovec_dati_mu[0]->Scale(1./histovec_dati_mu[0]->Integral());
    histovec_mc_mu_pre[0]->Scale(1./histovec_mc_mu_pre[0]->Integral());
    histovec_mc_mu_post[0]->Scale(1./histovec_mc_mu_post[0]->Integral());

    histovec_dati_mu[0]->Write("Eint_dati_nomediana", TObject::kOverwrite);
    histovec_mc_mu_pre[0]->Write("Eint_mc_nomediana", TObject::kOverwrite);
    histovec_mc_mu_post[0]->Write("Eint_mc_corretto_nomediana", TObject::kOverwrite);

    TCanvas *canvas2 = new TCanvas("Confronto_Eint_dati_mc_mu_mediana", " ");
    histovec_dati_mu[1]->Draw("P");
    histovec_dati_mu[1]->SetLineWidth(2);
    histovec_mc_mu_pre[1]->Draw("hist same");
    histovec_mc_mu_pre[1]->SetLineWidth(2);
    histovec_mc_mu_pre[1]->SetLineColor(kRed);
    histovec_mc_mu_pre[1]->SetLineStyle(kDashed);
    histovec_mc_mu_post[1]->Draw("hist same");
    histovec_mc_mu_post[1]->SetLineColor(kGreen);
    histovec_mc_mu_post[1]->SetLineWidth(2);
    TLegend *l5 = new TLegend();
    l5->AddEntry(histovec_dati_mu[1], "MUON DATA");
    l5->AddEntry(histovec_mc_mu_pre[1], "MUON MC"); 
    l5->AddEntry(histovec_mc_mu_post[1], "MUON MC CORRECTED");
    l5->Draw("same");
    canvas2->Write(0, TObject::kOverwrite);

    histovec_dati_mu[1]->Scale(1./histovec_dati_mu[1]->Integral());
    histovec_mc_mu_pre[1]->Scale(1./histovec_mc_mu_pre[1]->Integral());
    histovec_mc_mu_post[1]->Scale(1./histovec_mc_mu_post[1]->Integral());  
    
    histovec_dati_mu[1]->Write("Eint_dati_mediana", TObject::kOverwrite);
    histovec_mc_mu_pre[1]->Write("Eint_mc_mediana", TObject::kOverwrite);
    histovec_mc_mu_post[1]->Write("Eint_mc_corretto_mediana", TObject::kOverwrite);
    }
    


////////////////  Trasformazione delle reference curves di ArgoNeut da TProfile a TH1D /////////////////////////////////////////////////////////////////////////////////////////////////////////
   int update_reference = false;

   if(update_reference){
    TFile *fileReference = TFile::Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/THdedx_copy.root", "UPDATE");
    auto dedx_range_pro_AN = (TProfile*)fileReference->Get("dedx_range_pro");
    auto dedx_range_mu_AN  = (TProfile*)fileReference->Get("dedx_range_mu");
    
    TH1D *hdedx_range_pro_AN = new TH1D("hdedx_range_pro_AN", "", dedx_range_pro_AN->GetNbinsX(), dedx_range_pro_AN->GetXaxis()->GetXmin(), dedx_range_pro_AN->GetXaxis()->GetXmax());
    TH1D *hdedx_range_mu_AN = new TH1D("hdedx_range_mu_AN", "", dedx_range_mu_AN->GetNbinsX(), dedx_range_mu_AN->GetXaxis()->GetXmin(), dedx_range_mu_AN->GetXaxis()->GetXmax());

    for(int i=1; i<=hdedx_range_pro_AN->GetNbinsX(); i++)
    {
        hdedx_range_pro_AN->SetBinContent(i,dedx_range_pro_AN->GetBinContent(i));
        hdedx_range_pro_AN->SetBinError(i,dedx_range_pro_AN->GetBinError(i));
    }

    for(int i=1; i<=hdedx_range_mu_AN->GetNbinsX(); i++)
    {
        hdedx_range_mu_AN->SetBinContent(i,dedx_range_mu_AN->GetBinContent(i));
        hdedx_range_mu_AN->SetBinError(i,dedx_range_mu_AN->GetBinError(i));
    }
    hdedx_range_mu_AN->Write(0,TObject::kOverwrite);
    hdedx_range_pro_AN->Write(0,TObject::kOverwrite);

    TFile *file_new_ref = TFile::Open("hReference.root", "UPDATE");
    file_new_ref->cd();
    hdedx_range_pro_AN->Write(0,TObject::kOverwrite);
    hdedx_range_mu_AN->Write(0,TObject::kOverwrite);

    file_new_ref->Close();
    }
}




void plotMPVcorrected()
{
gStyle->SetOptStat(0);

//TFile *f= new TFile("mpvPlot_corrected.root", "RECREATE");
TFile *f= new TFile("mpvPlot_corrected_new.root", "RECREATE");
TDirectory *d = (TDirectory*)f->mkdir("muon");
TDirectory *d1 = (TDirectory*)f->mkdir("proton");

std::array<std::string,3> variFiletype={"", "25only_", "specific_"};

for(std::string filetype : variFiletype)
{

std::array<std::string,2> particles={"muon", "proton"};
for(std::string particle : particles){

 
    int Nbins;
    double x_high;
    if(filetype=="" && particle=="muon" ){Nbins=500; x_high=150;}
    else if(filetype=="" && particle=="proton"){Nbins=200; x_high=60;}
    else if(filetype=="25only_" || filetype=="specific_"){Nbins=85; x_high=25.5;}


/////////////////////////////////// CREAZIONE DELLE NUOVE REFERENCE CURVES di ICARUS con dEdx corretto ///////////////////////////////////////////////////
    TH1D *hmpv_mc_muon_corrected_fwhm = new TH1D(Form("%smc_muon_corrected_fwhm", filetype.c_str()), "", Nbins, 0, x_high);
    TH1D *hmpv_mc_proton_corrected_fwhm = new TH1D(Form("%smc_proton_corrected_fwhm", filetype.c_str()), "", Nbins, 0, x_high);

//////////////////// importo risultati fit con il dEdx corretto //////////////////////////////////////////////////////////////

std::string path; 
  



std::vector<size_t> colonne = {0, 1, 2};
path=Form("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/%sconfronto_%s_dati_corrected.txt", filetype.c_str(), particle.c_str());
auto dati = LeggiColonneDaFile(path.c_str(), colonne, "all");

std::vector<double>& x    = dati[0];
std::vector<double>& y    = dati[1];
std::vector<double>& yerr = dati[2];
std::vector<double> xerr ;

std::vector<size_t> colonne_ref={0,1,2,9};
path=Form("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/%sconfronto_%s_mc_corrected.txt", filetype.c_str(), particle.c_str());
auto ref = LeggiColonneDaFile(path.c_str(), colonne_ref, "all");

std::vector<double>& x_ref    = ref[0];
std::vector<double>& y_ref    = ref[1];
std::vector<double>& yerr_ref = ref[2];
std::vector<double> xerr_ref ;
std::vector<double>& fwhm = ref[3];

/////////////////////////////////   updating reference file with FWHM/2.354 errors for various MC dEdx corrected ///////////////////////////////////////////////////////////////////
/*
if(particle=="muon")
{
    for(int bin=3; bin<hmpv_mc_muon_corrected_fwhm->GetNbinsX(); bin++)
    {
        if(filetype=="25only_" && x_ref[bin-3]==20.85){y_ref[bin-3]=2.11906;} //corregge un punto anomalo 
        if(filetype=="specific_" && x_ref[bin-3]==10.65){y_ref[bin-3]=2.55294;} //corregge un punto anomalo
        hmpv_mc_muon_corrected_fwhm->SetBinContent(bin, y_ref[bin-3]);
        hmpv_mc_muon_corrected_fwhm->SetBinError(bin, fwhm[bin-3]/2.354 );
    }
}
else if(particle=="proton")
{
    for(int bin=3; bin<hmpv_mc_proton_corrected_fwhm->GetNbinsX(); bin++)
    {
        hmpv_mc_proton_corrected_fwhm->SetBinContent(bin, y_ref[bin-3]);
        hmpv_mc_proton_corrected_fwhm->SetBinError(bin, fwhm[bin-3]/2.354 );
    }
}
TFile *file_ref = TFile::Open("hReference.root", "UPDATE");
file_ref->cd();
if(particle=="muon")
{
    hmpv_mc_muon_corrected_fwhm->Write(0,TObject::kOverwrite);
}
else if(particle=="proton")
{   
    hmpv_mc_proton_corrected_fwhm->Write(0,TObject::kOverwrite);
}

file_ref->Close();

*/
///////////////////////////// costruisco i TGraph con gli MPV DATI e MC corretto ///////////////////////////////////////////////////

TDirectory *d3 = (TDirectory*)f->Get(particle.c_str());
d3->cd();

xerr.clear();
xerr_ref.clear();
for(int i=0; i<x.size(); i++)
{
    xerr.push_back(0.15);
    xerr_ref.push_back(0.15);
}

TGraphErrors *g = new TGraphErrors(x.size(), x.data(), y.data(), xerr.data(), yerr.data());

TGraphErrors *g_ref = new TGraphErrors(x_ref.size(), x_ref.data(), y_ref.data(), xerr_ref.data(), yerr_ref.data());


TCanvas *canvas = new TCanvas();
g_ref->GetXaxis()->SetTitle("Residual Range [cm]");
g_ref->GetYaxis()->SetTitle("MPV dE/dx [MeV/cm]");
g_ref->Draw("P");
g_ref->SetMarkerStyle(7);
g_ref->SetMarkerColor(kRed);
g_ref->SetLineColor(kRed);
g->Draw("P same");
g->SetMarkerStyle(7);
g->SetMarkerColor(kBlack);
g->SetLineColor(kBlack);
TLegend* l = new TLegend(0.65, 0.75, 0.85, 0.89);
l->AddEntry(g, "DATI");
l->AddEntry(g_ref, "MC");
l->Draw("same");
canvas->Write(Form("%sConfrontoDatiMC", filetype.c_str()));

//g_ref->Write("mc_mpvs");
//g->Write("dati_mpvs");


///////////////////////////////////////////// differenza DATI MC con il dEdx corretto ////////////////////////////////////////////
std::vector<double> differenza;
std::vector<double> err_differenza;


TGraphErrors *g_diff = new TGraphErrors();

TGraphErrors *trend = new TGraphErrors();

TGraphErrors *gdiff_dEdx = new TGraphErrors();


std::vector<double> dummy_valori;
std::vector<double> dummy_pesi;
std::vector<double> dummy_rr;



for(int i=2; i<x.size(); i++)//si skippano rr=0.75 e rr=1.05
{

    double diff = 2*(y[i]-y_ref[i])/(y[i]+y_ref[i]); //dati - mc
    double err_diff=propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] );

    int n = g_diff->GetN();
    g_diff->SetPoint(n, x[i], diff );
    g_diff->SetPointError(n,xerr[i],err_diff );

    int n_dEdx = gdiff_dEdx->GetN();
    gdiff_dEdx->SetPoint(n_dEdx, (y[i]+y_ref[i])/2, diff );
    gdiff_dEdx->SetPointError(n_dEdx, 0, err_diff );
    

    //RAGGRUPPARE I PUNTI IN CLUSTER
    if(i>=2)
    {
        dummy_valori.push_back(2.*(y[i]-y_ref[i])/(y[i]+y_ref[i]));
        dummy_pesi.push_back(propagateError(y[i],yerr[i],y_ref[i],yerr_ref[i] ));
        dummy_rr.push_back(x[i]);
    }

    
    if(int(dummy_valori.size())==3)
    {
        int nn = trend->GetN();
        std::pair<double,double> media = mediaPesata(dummy_valori, dummy_pesi);
        trend->SetPoint(nn, dummy_rr[1], media.first);
        trend->SetPointError(nn, 0.,  media.second);

        dummy_valori.clear();
        dummy_pesi.clear();
        dummy_rr.clear();

    }
    
}

//gdiff_dEdx->Write("gdiff_dEdx");

//differenza media in funzione del dEdx
TCanvas *diff_dEdx = new TCanvas();
gdiff_dEdx->SetMarkerStyle(7);
gdiff_dEdx->GetXaxis()->SetLimits(0., 16.);
gdiff_dEdx->GetXaxis()->SetTitle("(dati+mc)/2 [MeV/cm]");
gdiff_dEdx->GetYaxis()->SetTitle("(dati-mc)/media");
gdiff_dEdx->Draw("AP");
TLine* zero = new TLine(0, 0, 16, 0);
zero->SetLineStyle(2);
zero->SetLineColor(kGray+2);
zero->Draw("same");
diff_dEdx->Write(Form("%sdiff_dEdx",filetype.c_str()));
gdiff_dEdx->Write(Form("h%sdiff_dEdx",filetype.c_str()));

//differenza/media in funzione del rr
TCanvas *diffDatiMC = new TCanvas();
g_diff->Draw("A p");
g_diff->GetXaxis()->SetLimits(0.,26.);
g_diff->SetMarkerStyle(7);
TLine* zero1 = new TLine(0, 0, 26, 0);
zero1->SetLineStyle(2);
zero1->SetLineColor(kGray+2);
zero1->Draw("same");
diffDatiMC->Write(Form("%sdifferenza_media", filetype.c_str()));
g_diff->Write(Form("h%sdifferenza_media", filetype.c_str()));

//differenza/media in funzione del rr rebinnato
TCanvas *diffSpaced = new TCanvas();
trend->Draw("A P");
trend->GetXaxis()->SetLimits(0.,26.);
trend->SetMarkerStyle(7);
TLine* zero2 = new TLine(0, 0, 26, 0);
zero2->SetLineStyle(2);
zero2->SetLineColor(kGray+2);
zero2->Draw("same");
diffSpaced->Write(Form("%strend",filetype.c_str()));
trend->Write(Form("h%strend",filetype.c_str()));

}

}

}

/*
#include "broken_simul_xz.C"
void chi2()
{

  //importing reference curves
  TFile *f = TFile::Open("hReference.root", "READ");
  auto muon_mc_corr = (TH1D*)f->Get("mc_muon_corrected_fwhm");  
  auto proton_mc_corr = (TH1D*)f->Get("mc_proton_corrected_fwhm");
  auto muon_mc = (TH1D*)f->Get("mc_muon_fwhm");
  auto proton_mc = (TH1D*)f->Get("mc_pro_fwhm");
  auto muon_ArgoNeut = (TH1D*)f->Get("hdedx_range_mu_AN");
  auto proton_ArgoNeut = (TH1D*)f->Get("hdedx_range_pro_AN");
  auto muon_specific = (TH1D*)f->Get("specific_mc_muon_corrected_fwhm");
  auto proton_specific = (TH1D*)f->Get("specific_mc_proton_corrected_fwhm");

TFile *outfile = new TFile("Chi2NEW.root", "RECREATE");
outfile->cd();


///////////// quante hit sotto 1,5 MeV/cm //////////////////////////////////////////////////////////////////
std::array<std::string,2> particles ={"muon", "proton"};
std::array<std::string,2> dms = {"dati", "mc"};
for(std::string particle : particles)
{
    for(std::string dm : dms)
    {
        TH1D *h_dEdxtot = new TH1D(Form("h_dEdxtot_%s_%s", particle.c_str(), dm.c_str() ), "", 100, 0, 50);

        EventsData dat = load_data(dm, particle, "Np", "full", false, true, "");
        cout << dat.tree->GetEntries() << " total tracks" << endl;

        for(int track=0; track<dat.tree->GetEntries(); track++ )
        {
            dat.tree->GetEntry(track);
            for(int hit=0; hit<int(dat.track.dE->size()); hit++)
            {
                h_dEdxtot->Fill(dat.track.dE->at(hit));
            }
        }
        outfile->cd();
        h_dEdxtot->Write();

    }
}
outfile->cd();
///////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////// plot errori e reference curve ArgoNeut e ICARUS a confronto /////////////////////////////////////////////
  bool plot=true;

  if(plot){
  std::vector<double> err_mu_ICARUS;
  std::vector<double> x_mu_ICARUS;
  std::vector<double> err_pro_ICARUS;
  std::vector<double> x_pro_ICARUS;
  std::vector<double> err_mu_ArgoNeut;
  std::vector<double> x_mu_ArgoNeut;
  std::vector<double> err_pro_ArgoNeut;
  std::vector<double> x_pro_ArgoNeut;
  std::vector<double> err_mu_ArgoNeut_all;
  std::vector<double> err_pro_ArgoNeut_all;

  //0.04231 + 0.0001783 * trkdedx[i] * trkdedx[i] errore sulla singola hit 

  for(int i=0; i<muon_ArgoNeut->GetNbinsX(); i++)
  {
    double err = muon_ArgoNeut->GetBinError(i)/muon_ArgoNeut->GetBinContent(i)*100.;
    err_mu_ArgoNeut.push_back(err);
    x_mu_ArgoNeut.push_back(muon_ArgoNeut->GetBinCenter(i));
  }
  for(int i=0; i<proton_ArgoNeut->GetNbinsX(); i++)
  {
        double err = proton_ArgoNeut->GetBinError(i)/proton_ArgoNeut->GetBinContent(i)*100.;
        err_pro_ArgoNeut.push_back(err);
        x_pro_ArgoNeut.push_back(proton_ArgoNeut->GetBinCenter(i));
    }

  for(int i=0; i<muon_ArgoNeut->GetNbinsX(); i++)
  {
    double err = muon_ArgoNeut->GetBinError(i)/muon_ArgoNeut->GetBinContent(i)*100.;
    double err_hit = (0.04231 + 0.0001783 * muon_ArgoNeut->GetBinContent(i) * muon_ArgoNeut->GetBinContent(i))*100;
    err_mu_ArgoNeut_all.push_back(std::sqrt(err*err + err_hit*err_hit));
    x_mu_ArgoNeut.push_back(muon_ArgoNeut->GetBinCenter(i));
  }
  for(int i=0; i<proton_ArgoNeut->GetNbinsX(); i++)
  {
    double err = proton_ArgoNeut->GetBinError(i)/proton_ArgoNeut->GetBinContent(i)*100.;
    double err_hit = (0.04231 + 0.0001783 * proton_ArgoNeut->GetBinContent(i) * proton_ArgoNeut->GetBinContent(i))*100;
    err_pro_ArgoNeut_all.push_back(std::sqrt(err*err + err_hit*err_hit));
    x_pro_ArgoNeut.push_back(proton_ArgoNeut->GetBinCenter(i));
  }

  for(int i=0; i<muon_mc_corr->GetNbinsX(); i++)
  {
    err_mu_ICARUS.push_back(muon_mc_corr->GetBinError(i)/muon_mc_corr->GetBinContent(i)*100.);
    x_mu_ICARUS.push_back(muon_mc_corr->GetBinCenter(i));
  }
  for(int i=0; i<proton_mc_corr->GetNbinsX(); i++)
  {
    err_pro_ICARUS.push_back(proton_mc_corr->GetBinError(i)/proton_mc_corr->GetBinContent(i)*100.);
    x_pro_ICARUS.push_back(proton_mc_corr->GetBinCenter(i));
  }

  TGraph *g_err_mu_ArgoNeut = new TGraph(err_mu_ArgoNeut.size(), &x_mu_ArgoNeut[0], &err_mu_ArgoNeut[0]);
  TGraph *g_err_pro_ArgoNeut = new TGraph(err_pro_ArgoNeut.size(), &x_pro_ArgoNeut[0], &err_pro_ArgoNeut[0]);
  TGraph *g_err_mu_ICARUS = new TGraph(err_mu_ICARUS.size(), &x_mu_ICARUS[0], &err_mu_ICARUS[0]);
  TGraph *g_err_pro_ICARUS = new TGraph(err_pro_ICARUS.size(), &x_pro_ICARUS[0], &err_pro_ICARUS[0]);
  TGraph *g_err_mu_ArgoNeut_all = new TGraph(err_mu_ArgoNeut_all.size(), &x_mu_ArgoNeut[0], &err_mu_ArgoNeut_all[0]);
  TGraph *g_err_pro_ArgoNeut_all = new TGraph(err_pro_ArgoNeut_all.size(), &x_pro_ArgoNeut[0], &err_pro_ArgoNeut_all[0]);

  TCanvas *canvas_err = new TCanvas("canvas_err");
  g_err_mu_ArgoNeut->SetMarkerStyle(7);
  g_err_mu_ArgoNeut->SetMarkerColor(kAzure);
  g_err_pro_ArgoNeut->SetMarkerStyle(7);
  g_err_pro_ArgoNeut->SetMarkerColor(kGreen+2);
  g_err_mu_ICARUS->SetMarkerStyle(7);
  g_err_mu_ICARUS->SetMarkerColor(kOrange-3);
  g_err_pro_ICARUS->SetMarkerStyle(7);  
  g_err_pro_ICARUS->SetMarkerColor(kRed);
  g_err_pro_ICARUS->Draw("AP");
  g_err_pro_ICARUS->GetXaxis()->SetTitle("Residual Range [cm]");
  g_err_pro_ICARUS->GetYaxis()->SetTitle("Relative Error [%]");
  g_err_mu_ICARUS->Draw("P same");
  g_err_mu_ArgoNeut->Draw("P same");
  g_err_pro_ArgoNeut->Draw("P same");
  g_err_pro_ICARUS->GetXaxis()->SetRangeUser(0.7, 25.5);
  g_err_pro_ICARUS ->GetYaxis()->SetRangeUser(0., 50.);
  g_err_mu_ArgoNeut_all->SetMarkerStyle(7);
  g_err_mu_ArgoNeut_all->SetMarkerColor(kAzure+8);
  g_err_mu_ArgoNeut_all->SetLineColor(kAzure+8);
  g_err_mu_ArgoNeut_all->Draw("P same");
  g_err_pro_ArgoNeut_all->SetMarkerStyle(7);
  g_err_pro_ArgoNeut_all->SetMarkerColor(kGreen-7);
  g_err_pro_ArgoNeut_all->SetLineColor(kGreen-7);
  g_err_pro_ArgoNeut_all->Draw("P same");
  TLegend *l_err = new TLegend();
  l_err->AddEntry(g_err_mu_ArgoNeut, "ArgoNeut muons");
  l_err->AddEntry(g_err_pro_ArgoNeut, "ArgoNeut protons");
  l_err->AddEntry(g_err_mu_ICARUS, "ICARUS MC corrected muons");
  l_err->AddEntry(g_err_pro_ICARUS, "ICARUS MC corrected protons");
    l_err->AddEntry(g_err_mu_ArgoNeut_all, "ArgoNeut muons (with hit error)");
    l_err->AddEntry(g_err_pro_ArgoNeut_all, "ArgoNeut protons (with hit error)");   
  l_err->Draw("same");
  canvas_err->Write(0, TObject::kOverwrite);

  TCanvas *c = new TCanvas("references_mu");
  muon_ArgoNeut->Draw("p");
  muon_ArgoNeut->GetXaxis()->SetRangeUser(0.7, 25.5);
  muon_ArgoNeut->GetYaxis()->SetRangeUser(1.5, 10.);
  muon_ArgoNeut->GetXaxis()->SetTitle("Residual Range [cm]");
  muon_ArgoNeut->GetYaxis()->SetTitle("dEdx MPV [MeV/cm]");
  muon_ArgoNeut->SetLineWidth(2);
  muon_ArgoNeut->SetMarkerStyle(7);
  muon_ArgoNeut->SetMarkerColor(kAzure);
  muon_ArgoNeut->SetLineColor(kAzure-4);
  muon_specific->Draw("p same");
  muon_specific->SetLineWidth(2);
  muon_specific->SetMarkerStyle(7);
  muon_specific->SetMarkerColor(kOrange-3);
  muon_specific->SetLineColor(kOrange);
  TLegend * l = new TLegend();
  l->AddEntry(muon_ArgoNeut, "ArgoNeut");
  l->AddEntry(muon_specific, "ICARUS MC corrected #mu specific");
  l->Draw("same");
  c->Write(0, TObject::kOverwrite);

  TCanvas *c1 = new TCanvas("references_pro");
  proton_ArgoNeut->Draw("p");
  proton_ArgoNeut->GetXaxis()->SetRangeUser(0.7, 25.5);
  proton_ArgoNeut->GetYaxis()->SetRangeUser(3.5, 24.);
  proton_ArgoNeut->GetXaxis()->SetTitle("Residual Range [cm]");
  proton_ArgoNeut->GetYaxis()->SetTitle("dEdx MPV [MeV/cm]");
  proton_ArgoNeut->SetLineWidth(2);
  proton_ArgoNeut->SetMarkerStyle(7);
  proton_ArgoNeut->SetMarkerColor(kGreen+2);
  proton_ArgoNeut->SetLineColor(kGreen-7);
  proton_specific->Draw("p same");
  proton_specific->SetLineWidth(2);
  proton_specific->SetMarkerStyle(7);
  proton_specific->SetMarkerColor(kRed);
  proton_specific->SetLineColor(kRed-7);
  TLegend * l1 = new TLegend();
  l1->AddEntry(proton_ArgoNeut, "ArgoNeut");
  l1->AddEntry(proton_specific, "ICARUS MC corrected p specific");
  l1->Draw("same");
  c1->Write(0, TObject::kOverwrite);
  }



///////////////////////// calcolo del chi2 ////////////////////////////////////////////////////////////////////////
std::array<std::string,2> data_o_mc ={"dati", "mc"};
for(std::string str : data_o_mc)
{

  std::array<std::string,2> corrected_or_not ={"", "corr"};
  for(std::string corr : corrected_or_not ){

    std::array<std::string,4> options = {"", "withCOS", "YZvar", "broken"};
    for(std::string option : options){

        if(str=="dati" && (corr!="" || option!=""))continue;

        bool isCorrected=true;
        bool isTrue1muNp=true;
        if(corr=="")isCorrected=false;

        cout << str << " " << corr << " " << option << endl; 

  std::string option_scritto, corr_scritto;
  if(option=="")option_scritto=option;
  else option_scritto=option+"_";

  if(corr=="")corr_scritto=corr;
  else corr_scritto=corr+"_";

  TH1D *hchi2_mu_ArgoNeut = new TH1D(Form("%s_%s%smu_ArgoNeut", str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);
  TH1D *hchi2_mu_mcICARUS = new TH1D(Form("%s_%s%smu_mcICARUS",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);
  TH1D *hchi2_mu_mcICARUScorr = new TH1D(Form("%s_%s%smu_mcICARUScorr",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);
  TH1D *hchi2_mu_mcICARUSspec = new TH1D(Form("%s_%s%smu_mcICARUSspec",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);

  TH1D *hchi2_pro_ArgoNeut = new TH1D(Form("%s_%s%spro_ArgoNeut",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);
  TH1D *hchi2_pro_mcICARUS = new TH1D(Form("%s_%s%spro_mcICARUS",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);
  TH1D *hchi2_pro_mcICARUScorr = new TH1D(Form("%s_%s%spro_mcICARUScorr",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);
  TH1D *hchi2_pro_mcICARUSspec = new TH1D(Form("%s_%s%spro_mcICARUSspec",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 150, 0, 150);

  TH1D *hchi2_mu_as_pro_ArgoNeut = new TH1D(Form("%s_%s%smu_as_pro_ArgoNeut",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);
  TH1D *hchi2_mu_as_pro_mcICARUS = new TH1D(Form("%s_%s%smu_as_pro_mcICARUS",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);
  TH1D *hchi2_mu_as_pro_mcICARUScorr = new TH1D(Form("%s_%s%smu_as_pro_mcICARUScorr",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);
  TH1D *hchi2_mu_as_pro_mcICARUSspec = new TH1D(Form("%s_%s%smu_as_pro_mcICARUSspec",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);

  TH1D *hchi2_pro_as_mu_ArgoNeut = new TH1D(Form("%s_%s%spro_as_mu_ArgoNeut",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);
  TH1D *hchi2_pro_as_mu_mcICARUS = new TH1D(Form("%s_%s%spro_as_mu_mcICARUS",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);
  TH1D *hchi2_pro_as_mu_mcICARUScorr = new TH1D(Form("%s_%s%spro_as_mu_mcICARUScorr",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);
  TH1D *hchi2_pro_as_mu_mcICARUSspec = new TH1D(Form("%s_%s%spro_as_mu_mcICARUSspec",str.c_str(), option_scritto.c_str(), corr_scritto.c_str()), "", 300, 0, 300);


  if(option=="broken")
  {

    EventsData dat_mu_broken;
    if(corr=="corr")dat_mu_broken = load_data("mc", "muon", "Np", "25only", true, false, "YZvar");
    if(corr=="")dat_mu_broken = load_data("mc", "muon", "Np", "25only", false, false, "YZvar");
    cout << dat_mu_broken.tree->GetEntries() << " total tracks" << endl;
    outfile->cd();
    std::vector<std::vector<double>> dEdx_broken;
    std::vector<std::vector<double>> rr_broken;
    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> calo_broken = simul_broken_xz(dat_mu_broken);
    dEdx_broken=calo_broken.first;
    rr_broken=calo_broken.second;

    for(int track=0; track<int(dEdx_broken.size()); track++)
    {
        std::vector<double> dEdx = dEdx_broken[track];
        std::vector <double> RR = rr_broken[track];
        hchi2_mu_ArgoNeut->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_ArgoNeut, proton_ArgoNeut)[0] );
        hchi2_mu_mcICARUS->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc, proton_mc)[0] );
        hchi2_mu_mcICARUScorr->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc_corr, proton_mc_corr)[0] );
        hchi2_mu_mcICARUSspec->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_specific, proton_specific)[0] );
        hchi2_mu_as_pro_ArgoNeut->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_ArgoNeut, proton_ArgoNeut)[1] );
        hchi2_mu_as_pro_mcICARUS->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc, proton_mc)[1] );
        hchi2_mu_as_pro_mcICARUScorr->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc_corr, proton_mc_corr)[1] );
        hchi2_mu_as_pro_mcICARUSspec->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_specific, proton_specific)[1] );
    }

        hchi2_mu_ArgoNeut->Scale(1./hchi2_mu_ArgoNeut->Integral());
        hchi2_mu_mcICARUS->Scale(1./hchi2_mu_mcICARUS->Integral());
        hchi2_mu_mcICARUScorr->Scale(1./hchi2_mu_mcICARUScorr->Integral());
        hchi2_mu_mcICARUSspec->Scale(1./hchi2_mu_mcICARUSspec->Integral());
        hchi2_mu_as_pro_ArgoNeut->Scale(1./hchi2_mu_as_pro_ArgoNeut->Integral());
        hchi2_mu_as_pro_mcICARUS->Scale(1./hchi2_mu_as_pro_mcICARUS->Integral());
        hchi2_mu_as_pro_mcICARUScorr->Scale(1./hchi2_mu_as_pro_mcICARUScorr->Integral());
        hchi2_mu_as_pro_mcICARUSspec->Scale(1./hchi2_mu_as_pro_mcICARUSspec->Integral());

        hchi2_mu_ArgoNeut->Write(0, TObject::kOverwrite);
        hchi2_mu_mcICARUS->Write(0, TObject::kOverwrite);
        hchi2_mu_mcICARUScorr->Write(0, TObject::kOverwrite);
        hchi2_mu_mcICARUSspec->Write(0, TObject::kOverwrite);
        hchi2_mu_as_pro_ArgoNeut->Write(0, TObject::kOverwrite);
        hchi2_mu_as_pro_mcICARUS->Write(0, TObject::kOverwrite);
        hchi2_mu_as_pro_mcICARUScorr->Write(0, TObject::kOverwrite);
        hchi2_mu_as_pro_mcICARUSspec->Write(0, TObject::kOverwrite);

    continue;

  }


  EventsData dat = load_data(str, "muon", "Np", "full", isCorrected, isTrue1muNp, option);
  cout << dat.tree->GetEntries() << " total tracks" << endl;
  outfile->cd();
  cout << "chi2 of muons using " << str << " " << option <<  " " << corr << endl;
  for(int track=0; track<dat.tree->GetEntries(); track++)
  {
        dat.tree->GetEntry(track);

        std::vector<double> dEdx;
        std::vector<double> RR;

        for(int hit=0; hit<int(dat.track.rr->size()); hit++)
        {
            dEdx.push_back(dat.track.dE->at(hit));
            RR.push_back(dat.track.rr->at(hit));
        }

        hchi2_mu_ArgoNeut->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_ArgoNeut, proton_ArgoNeut)[0] );
        hchi2_mu_mcICARUS->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc, proton_mc)[0] );
        hchi2_mu_mcICARUScorr->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc_corr, proton_mc_corr)[0] );
        hchi2_mu_mcICARUSspec->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_specific, proton_specific)[0] );
        hchi2_mu_as_pro_ArgoNeut->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_ArgoNeut, proton_ArgoNeut)[1] );
        hchi2_mu_as_pro_mcICARUS->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc, proton_mc)[1] );
        hchi2_mu_as_pro_mcICARUScorr->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc_corr, proton_mc_corr)[1] );
        hchi2_mu_as_pro_mcICARUSspec->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_specific, proton_specific)[1] );
  }

    
  
  EventsData dat_pro = load_data(str, "proton", "Np", "full", isCorrected, isTrue1muNp, option);
  cout << dat_pro.tree->GetEntries() << " total tracks" << endl;
  outfile->cd();
  cout << "chi2 of protons using " << str << " " << option << " " << corr << endl;
  for(int track=0; track<dat_pro.tree->GetEntries(); track++)
  {
        dat_pro.tree->GetEntry(track);

        std::vector<double> dEdx;
        std::vector<double> RR;

        for(int hit=0; hit<int(dat_pro.track.rr->size()); hit++)
        {
            dEdx.push_back(dat_pro.track.dE->at(hit));
            RR.push_back(dat_pro.track.rr->at(hit));
        }

        hchi2_pro_ArgoNeut->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_ArgoNeut, proton_ArgoNeut)[1] );
        hchi2_pro_mcICARUS->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc, proton_mc)[1] );
        hchi2_pro_mcICARUScorr->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc_corr, proton_mc_corr)[1] );
        hchi2_pro_mcICARUSspec->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_specific, proton_specific)[1] );
        hchi2_pro_as_mu_ArgoNeut->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_ArgoNeut, proton_ArgoNeut)[0] );
        hchi2_pro_as_mu_mcICARUS->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc, proton_mc)[0] );
        hchi2_pro_as_mu_mcICARUScorr->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_mc_corr, proton_mc_corr)[0] );
        hchi2_pro_as_mu_mcICARUSspec->Fill( chi2_ALG(dEdx, RR, 0.6, 25., muon_specific, proton_specific)[0] );
  }
  
  //TDirectory *d = (TDirectory*)outfile->mkdir(str.c_str());
  //d->cd();

  hchi2_mu_ArgoNeut->Scale(1./hchi2_mu_ArgoNeut->Integral());
  hchi2_mu_mcICARUS->Scale(1./hchi2_mu_mcICARUS->Integral());
  hchi2_mu_mcICARUScorr->Scale(1./hchi2_mu_mcICARUScorr->Integral());
  hchi2_mu_mcICARUSspec->Scale(1./hchi2_mu_mcICARUSspec->Integral());
  hchi2_mu_ArgoNeut->Write();
  hchi2_mu_mcICARUS->Write();
  hchi2_mu_mcICARUScorr->Write();
  hchi2_mu_mcICARUSspec->Write();
  hchi2_mu_as_pro_ArgoNeut->Scale(1./hchi2_mu_as_pro_ArgoNeut->Integral());
  hchi2_mu_as_pro_mcICARUS->Scale(1./hchi2_mu_as_pro_mcICARUS->Integral());
  hchi2_mu_as_pro_mcICARUScorr->Scale(1./hchi2_mu_as_pro_mcICARUScorr->Integral());
  hchi2_mu_as_pro_mcICARUSspec->Scale(1./hchi2_mu_as_pro_mcICARUSspec->Integral());
  hchi2_mu_as_pro_ArgoNeut->Write();
  hchi2_mu_as_pro_mcICARUS->Write();
  hchi2_mu_as_pro_mcICARUScorr->Write();
  hchi2_mu_as_pro_mcICARUSspec->Write();

if(option!="broken"){
    hchi2_pro_ArgoNeut->Scale(1./hchi2_pro_ArgoNeut->Integral()); 
    hchi2_pro_mcICARUS->Scale(1./hchi2_pro_mcICARUS->Integral());
    hchi2_pro_mcICARUScorr->Scale(1./hchi2_pro_mcICARUScorr->Integral());
    hchi2_pro_mcICARUSspec->Scale(1./hchi2_pro_mcICARUSspec->Integral());
    hchi2_pro_as_mu_ArgoNeut->Scale(1./hchi2_pro_as_mu_ArgoNeut->Integral());
    hchi2_pro_as_mu_mcICARUS->Scale(1./hchi2_pro_as_mu_mcICARUS->Integral());
    hchi2_pro_as_mu_mcICARUScorr->Scale(1./hchi2_pro_as_mu_mcICARUScorr->Integral());
    hchi2_pro_as_mu_mcICARUSspec->Scale(1./hchi2_pro_as_mu_mcICARUSspec->Integral());
    hchi2_pro_as_mu_ArgoNeut->Write();
    hchi2_pro_as_mu_mcICARUS->Write();
    hchi2_pro_as_mu_mcICARUScorr->Write();
    hchi2_pro_as_mu_mcICARUSspec->Write();
    hchi2_pro_ArgoNeut->Write();
    hchi2_pro_mcICARUS->Write();
    hchi2_pro_mcICARUScorr->Write();
    hchi2_pro_mcICARUSspec->Write();
}
            }//normale, with cosmic o YZvar
        }//scaling or not
    }//dati o mc


    outfile->Close();
  
}

int findBinMpv(TH1D* h)
{
    int binMax=0;
    double maxVal=0;
    for(int j=1; j<h->GetNbinsX(); j++)
    {
        double y=h->GetBinContent(j);
        if(y>maxVal)
        {
            maxVal=y;
            binMax=j;
        }
    }
    return binMax;
}

void plotHisto(std::string name1, std::string name2, std::string name3="", std::string name4="", const char *canvasName="canvas")
{
    TFile *f = TFile::Open("Chi2NEW.root", "UPDATE");
    TH1D * h1 = (TH1D*)f->Get(name1.c_str());
    TH1D * h2 = (TH1D*)f->Get(name2.c_str());
    TH1D *h3;
    if(name3!=""){h3 = (TH1D*)f->Get(name3.c_str());}
    TH1D *h4;
    if(name4!=""){h4 = (TH1D*)f->Get(name4.c_str());}

    double h1height=h1->GetBinContent(findBinMpv(h1));
    double h2height=h2->GetBinContent(findBinMpv(h2));
    double h3height=0;
    if(name3!=""){h3height=h3->GetBinContent(findBinMpv(h3));}

    double maxheight=h1height;
    if(h2height>h1height)
    {
        maxheight=h2height;
        if(h3height>h2height){maxheight=h3height;}
    }
    else if(h3height>h1height){maxheight=h3height;}


    TCanvas *c = new TCanvas(canvasName);
    h1->Draw("hist");
    h1->GetYaxis()->SetRangeUser(0.,maxheight+0.2);
    h1->GetXaxis()->SetTitle("#chi2");
    h1->GetYaxis()->SetTitle("counts (area normalized)");
    h1->SetLineWidth(2);
    h1->SetLineColor(kOrange);
    h2->Draw("hist same");
    h2->SetLineWidth(2);
    h2->SetLineColor(kBlue);
    if(name3!="")
    {
        h3->Draw("hist same");
        h3->SetLineWidth(2);
        h3->SetLineColor(kBlack);
    }
    if(name4!="")
    {
        h4->Draw("hist same");
        h4->SetLineWidth(2);
        h4->SetLineColor(kRed);
    }
    TLegend *l = new TLegend();

    std::string nome1, nome2, nome3, nome4;
    if(name4!=""){nome4=name4;}
    nome1=name1;
    nome2=name2;
    if(name3!=""){nome3=name3;}

    std::replace(nome1.begin(), nome1.end(), '_', ' ');

    std::replace(nome2.begin(), nome2.end(), '_', ' ');

    if(name3!="")
    {
        std::replace(nome3.begin(), nome3.end(), '_', ' ');
    }
    if(name4!="")
    {
        std::replace(nome4.begin(), nome4.end(), '_', ' ');
    }

    l->AddEntry(h1, nome1.c_str());
    l->AddEntry(h2, nome2.c_str());
    if(name3!=""){l->AddEntry(h3,nome3.c_str());}
    if(name4!=""){l->AddEntry(h4,nome4.c_str());}
    l->Draw("same");
    c->Draw();
    c->Write(0, TObject::kOverwrite);

    f->Close();
}


void plotBananaConCurve()
{
    TH2D *dedx_range = new TH2D("dedx_range_mu", "", 500, 0, 25, 300, 0, 30);
    EventsData dat = load_data("mc", "muon", "Np", "25only", false, true, "");
    
    for(int track=0; track<dat.tree->GetEntries(); track++)
    {
        dat.tree->GetEntry(track);
        
        for(int hit=1; hit<int(dat.track.rr->size())-1; hit++)
        {
            dedx_range->Fill(dat.track.rr->at(hit), dat.track.dE->at(hit));
        }
    }

    TFile *f = TFile::Open("hReference.root", "READ");
    auto *reference_ArgoNeut_originale = (TH1D*)f->Get("hdedx_range_mu_AN");
    auto *reference_ICARUS_originale = (TH1D*)f->Get("25only_mc_muon_corrected_fwhm");
    auto *ICARUS_mc_originale = (TH1D*)f->Get("mc_muon_fwhm");

    TH1D *reference_ArgoNeut =  (TH1D*)reference_ArgoNeut_originale->Clone("ref_AN");
    TH1D *reference_ICARUS = (TH1D*)reference_ICARUS_originale->Clone("ref_ICARUS");
    TH1D *ICARUS_mc = (TH1D*)ICARUS_mc_originale->Clone("ICARUS_mc");

    for(int i=1; i<=reference_ArgoNeut->GetNbinsX(); i++)
    {
        reference_ArgoNeut->SetBinError(i,0);
    }
    for(int i=1; i<=reference_ICARUS->GetNbinsX(); i++)
    {
        reference_ICARUS->SetBinError(i,0);
    }
    for(int i=1; i<ICARUS_mc->GetNbinsX(); i++)
    {
        ICARUS_mc->SetBinError(i,0);
    }
    

    TCanvas *c = new TCanvas();
    dedx_range->Draw("colz");
    dedx_range->GetXaxis()->SetTitle("Residual Range [cm]");
    dedx_range->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
    reference_ArgoNeut->Draw("p same");
    reference_ArgoNeut->SetMarkerStyle(7);
    reference_ArgoNeut->SetMarkerColor(kRed);
    reference_ArgoNeut->SetLineColor(kRed);
    reference_ArgoNeut->SetLineWidth(3);
    reference_ICARUS->Draw("p same");
    reference_ICARUS->SetMarkerColor(kBlack);
    reference_ICARUS->SetLineColor(kBlack);
    reference_ICARUS->SetLineWidth(3);
    reference_ICARUS->SetMarkerStyle(7);
    ICARUS_mc->Draw("p same");
    ICARUS_mc->SetMarkerStyle(7);
    ICARUS_mc->SetMarkerColor(kMagenta);
    ICARUS_mc->SetLineColor(kMagenta);
    ICARUS_mc->SetLineWidth(3);
    TLegend *l = new TLegend();
    l->AddEntry(reference_ArgoNeut, "ArgoNeut");
    l->AddEntry(ICARUS_mc, "ICARUS MC");
    l->AddEntry(reference_ICARUS, "Corrected ICARUS MC");
    l->Draw("same");

    TFile *file = TFile::Open("controlli.root", "UPDATE");

    c->Write("banana", TObject::kOverwrite);



}





*/


