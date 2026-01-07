import ROOT
import numpy as np
import math


def correction_function_25only(x, par):
    if x[0] <= 2:
        return 0.0127159
    elif x[0] >= 15:
        return -0.0235248
    else:
        return math.exp(par[0] + par[1]*x[0]) + par[2] + par[3] * math.exp(-0.5 * ((x[0] - par[4]) / par[5])**2)

def check_factors():
    correction_function = ROOT.TF1("correction_function", correction_function_25only, 0., 16., 6)
    correction_function.SetParameters(-1.08107e+00, -1.38302e+00, -2.84635e-02, 5.87297e-02, 7.17814e+00, -3.51461e+00)

    graph = ROOT.TGraph()

    for i in range(0,16000):
        factor= (2+correction_function.Eval(i/1000.))/(2-correction_function.Eval(i/1000.))
        graph.SetPoint(i, i/1000., factor)

    factors_mu = np.loadtxt("dump_controllo_chi2_mu.txt")
    factors_pro = np.loadtxt("dump_controllo_chi2_pro.txt")

    graph_mu = ROOT.TGraph()
    graph_pro = ROOT.TGraph()
    h_dEdx_mu = ROOT.TH1D("h_dEdx_mu", "", 300, 0, 30)
    h_dEdx_pro = ROOT.TH1D("h_dEdx_pro", "", 300, 0, 30)


    for i in range(len(factors_mu[:,0])):
        graph_mu.SetPoint(i, factors_mu[i,0], factors_mu[i,1])
        h_dEdx_mu.Fill(factors_mu[i,0])

    for i in range(len(factors_pro[:,0])):
        graph_pro.SetPoint(i, factors_pro[i,0], factors_pro[i,1])
        h_dEdx_pro.Fill(factors_pro[i,0])

    canvas = ROOT.TCanvas()

    graph.SetMarkerStyle(7)
    graph.SetMarkerColor(ROOT.kRed)
    graph_mu.SetMarkerStyle(2)
    graph_mu.SetMarkerColor(ROOT.kGreen)
    graph_pro.SetMarkerStyle(3)
    graph.Draw("AP")
    graph_mu.Draw("P same")
    graph_pro.Draw("P same")
    canvas.Draw()
    file = ROOT.TFile("controlli.root", "RECREATE")
    canvas.Write()

    h_dEdx_mu.Write()
    h_dEdx_pro.Write()


def stampa_chiavi():
    f=ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/Chi2NEW.root")
    keys=f.GetListOfKeys()

    for key in keys:
        print(key.GetName())

def plot_dEdx():
    file_confronto = ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC.root")
    file_confronto_corretto = ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_corrected.root")
    file_controllo = ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/controlli.root", "RECREATE")

    proton_dir_confronto = file_confronto.Get("proton")
    proton_dir_confronto_mc = proton_dir_confronto.Get("mc")
    proton_dir_confronto_dati = proton_dir_confronto.Get("dati")

    muon_dir_confronto = file_confronto.Get("muon")
    muon_dir_confronto_mc = muon_dir_confronto.Get("mc")
    muon_dir_confronto_dati = muon_dir_confronto.Get("dati")

    proton_dir_confronto_corr = file_confronto_corretto.Get("proton")  
    proton_dir_confronto_corr_mc = proton_dir_confronto_corr.Get("mc")

    muon_dir_confronto_corr = file_confronto_corretto.Get("muon")
    muon_dir_confronto_corr_mc = muon_dir_confronto_corr.Get("mc")

    
    dEdx_mu_mc = muon_dir_confronto_mc.Get("dEdx_100_150_tot")
    dEdx_mu_dati = muon_dir_confronto_dati.Get("dEdx_100_150_tot")
    dEdx_mu_mc_corr = muon_dir_confronto_corr_mc.Get("dEdx_100_150_tot")
    dEdx_pro_mc = proton_dir_confronto_mc.Get("dEdx_100_150_tot")
    dEdx_pro_dati = proton_dir_confronto_dati.Get("dEdx_100_150_tot")
    dEdx_pro_mc_corr = proton_dir_confronto_corr_mc.Get("dEdx_100_150_tot")

    dEdx_mu_mc.Scale(1./dEdx_mu_mc.Integral())
    dEdx_mu_dati.Scale(1./dEdx_mu_dati.Integral())
    dEdx_mu_mc_corr.Scale(1./dEdx_mu_mc_corr.Integral())

    print(dEdx_mu_dati.GetMean(), dEdx_mu_mc.GetMean(), dEdx_mu_mc_corr.GetMean()  )

    file_controllo.cd()

    c = ROOT.TCanvas()
    dEdx_mu_mc.Draw("hist")
    dEdx_mu_mc.SetLineColor(ROOT.kBlue)
    dEdx_mu_dati.Draw("p same")
    dEdx_mu_dati.SetMarkerStyle(2)
    dEdx_mu_dati.SetLineColor(ROOT.kRed)
    dEdx_mu_dati.SetMarkerColor(ROOT.kRed)
    dEdx_mu_dati.SetMarkerSize(0.5)
    l = ROOT.TLegend()
    l.AddEntry(dEdx_mu_mc, "MC")
    l.AddEntry(dEdx_mu_dati, "DATA")
    l.Draw("same")
    c.Write()

    c1 =ROOT.TCanvas()
    dEdx_mu_mc_corr.Draw("hist")
    dEdx_mu_mc_corr.SetLineColor(ROOT.kBlue)
    dEdx_mu_dati.Draw("p same")
    dEdx_mu_dati.SetMarkerStyle(2)
    dEdx_mu_dati.SetLineColor(ROOT.kRed)
    dEdx_mu_dati.SetMarkerColor(ROOT.kRed)
    dEdx_mu_dati.SetMarkerSize(0.5)
    l1 = ROOT.TLegend()
    l1.AddEntry(dEdx_mu_mc_corr, "MC corrected (full RR)")
    l1.AddEntry(dEdx_mu_dati, "DATA")
    l1.Draw("same")
    c1.Write()


def plotDirections():
    f = ROOT.TFile.Open("directions.root", "UPDATE")
    dirx_mu_dati = f.Get("h_dirx_mu")
    dirx_pro_dati = f.Get("h_dirx_pro")
    dirz_mu_dati = f.Get("h_dirz_mu")
    dirz_pro_dati = f.Get("h_dirz_pro")
    dirx_mu_mc = f.Get("h_dirx_mu_mc")
    dirx_pro_mc = f.Get("h_dirx_pro_mc")
    dirz_mu_mc = f.Get("h_dirz_mu_mc")
    dirz_pro_mc = f.Get("h_dirz_pro_mc")

    c = ROOT.TCanvas("c", "c", 800, 600)
    dirx_mu_dati.Draw("p")
    dirx_mu_dati.SetTitle("MUONS direction x")
    dirx_mu_dati.GetXaxis().SetTitle("direction x")
    dirx_mu_dati.GetYaxis().SetTitle("counts (area normalized)")
    dirx_mu_dati.SetMarkerStyle(7)
    dirx_mu_dati.SetLineColor(ROOT.kOrange)
    dirx_mu_dati.SetLineWidth(2)
    dirx_mu_dati.SetMarkerColor(ROOT.kOrange)

    dirx_mu_mc.Draw("p same")
    dirx_mu_mc.SetMarkerStyle(7)
    dirx_mu_mc.SetMarkerColor(ROOT.kAzure)
    dirx_mu_mc.SetLineColor(ROOT.kAzure)
    dirx_mu_mc.SetLineWidth(2)

    l = ROOT.TLegend()
    l.AddEntry(dirx_mu_dati, "Data")
    l.AddEntry(dirx_mu_mc, "MC")
    l.Draw("same")
    c.Write("dirx_mu",ROOT.TObject.kOverwrite)
########################################################################################
    c1 = ROOT.TCanvas("c1", "c1", 800, 600)
    dirx_pro_dati.Draw("p")
    dirx_pro_dati.SetTitle("PROTONS direction x")
    dirx_pro_dati.GetXaxis().SetTitle("direction x")
    dirx_pro_dati.GetYaxis().SetTitle("counts (area normalized)")
    dirx_pro_dati.SetMarkerStyle(7)
    dirx_pro_dati.SetMarkerColor(ROOT.kOrange)
    dirx_pro_dati.SetLineColor(ROOT.kOrange)
    dirx_pro_dati.SetLineWidth(2)

    dirx_pro_mc.Draw("p same")
    dirx_pro_mc.SetMarkerStyle(7)
    dirx_pro_mc.SetMarkerColor(ROOT.kAzure)
    dirx_pro_mc.SetLineColor(ROOT.kAzure)
    dirx_pro_mc.SetLineWidth(2)

    l1 = ROOT.TLegend()
    l1.AddEntry(dirx_pro_dati, "Data")
    l1.AddEntry(dirx_pro_mc, "MC")
    l1.Draw("same")
    c1.Write("dirx_pro",ROOT.TObject.kOverwrite)
##########################################################################
    c2 = ROOT.TCanvas("c2", "c2", 800, 600)

    dirz_mu_mc.SetTitle("MUONS direction z")
    dirz_mu_mc.SetLineColor(ROOT.kAzure)
    dirz_mu_mc.SetLineWidth(2)
    dirz_mu_mc.Draw("p")
    dirz_mu_mc.GetXaxis().SetTitle("direction z")
    dirz_mu_mc.GetYaxis().SetTitle("counts (area normalized)")
    dirz_mu_mc.SetMarkerStyle(7)
    dirz_mu_mc.SetMarkerColor(ROOT.kAzure)
       
    dirz_mu_dati.Draw("p same")
    dirz_mu_dati.SetMarkerStyle(7)
    dirz_mu_dati.SetMarkerColor(ROOT.kOrange)
    dirz_mu_dati.SetLineColor(ROOT.kOrange)
    dirz_mu_dati.SetLineWidth(2)
 
    l2 = ROOT.TLegend()
    l2.AddEntry(dirz_mu_dati, "Data")
    l2.AddEntry(dirz_mu_mc, "MC")
    l2.Draw("same")
    c2.Write("dirz_mu",ROOT.TObject.kOverwrite)  
#############################################################################
    c3 = ROOT.TCanvas("c3", "c3", 800, 600)

    dirz_pro_mc.SetLineColor(ROOT.kAzure)
    dirz_pro_mc.SetTitle("PROTONS direction z")
    dirz_pro_mc.SetLineWidth(2)
    dirz_pro_mc.Draw("p")
    dirz_pro_mc.GetXaxis().SetTitle("direction z")
    dirz_pro_mc.GetYaxis().SetTitle("counts (area normalized)")
    dirz_pro_mc.SetMarkerStyle(7)
    dirz_pro_mc.SetMarkerColor(ROOT.kAzure)
    dirz_pro_mc.SetLineColor(ROOT.kAzure)

    dirz_pro_dati.Draw("p same")
    dirz_pro_dati.SetMarkerStyle(7)
    dirz_pro_dati.SetMarkerColor(ROOT.kOrange)
    dirz_pro_dati.SetLineColor(ROOT.kOrange)
    dirz_pro_dati.SetLineWidth(2)

    l3 = ROOT.TLegend()
    l3.AddEntry(dirz_pro_dati, "Data")
    l3.AddEntry(dirz_pro_mc, "MC")
    l3.Draw("same")
    c3.Write("dirz_pro",ROOT.TObject.kOverwrite)  


def plotHIsto():
    f = ROOT.TFile.Open("mpvPlot.root", "READ")
    f1 = ROOT.TFile.Open("mpvPlot_above30driftAngle.root", "READ")

    g = f.Get("proton/gtrend")
    g1 = f1.Get("proton/gtrend")

    c = ROOT.TCanvas()

    g.Draw("AP")
    g.GetXaxis().SetTitle("Residual Range [cm]")
    g.GetYaxis().SetTitle("(DATA-MC)/average")
    g.SetLineColor(ROOT.kRed)
    g.SetLineWidth(2)
    g.SetMarkerStyle(7)
    g.SetMarkerColor(ROOT.kRed)
    g1.Draw("P same")
    g1.SetLineWidth(2)
    g1.SetLineColor(ROOT.kAzure)
    g1.SetMarkerStyle(7)
    g1.SetMarkerColor(ROOT.kAzure)
    
    l = ROOT.TLegend()
    l.AddEntry(g, "no condition on angle")
    l.AddEntry(g1, "angle to drift above 30 degrees")
    l.Draw("same")

    out = ROOT.TFile("plot_output.root", "RECREATE")
    c.Write()
    

    


    
def check_smoothness():

    x_mu = np.loadtxt("correzionimuon.txt")[:,0]
    y_mu = np.loadtxt("correzionimuon.txt")[:,1]
    x_pro = np.loadtxt("correzioniproton.txt")[:,0]
    y_pro = np.loadtxt("correzioniproton.txt")[:,1]

    g_mu = ROOT.TGraph()
    for i in range(len(x_mu)):
        g_mu.SetPoint(i, x_mu[i], y_mu[i])
    g_pro = ROOT.TGraph()
    for i in range(len(x_pro)):
        g_pro.SetPoint(i, x_pro[i], y_pro[i])    

    f = ROOT.TFile.Open("probDensities.root", "UPDATE")
    f.cd()
    c = ROOT.TCanvas()
    g_mu.Draw("AP")
    g_pro.Draw("P same")
    g_mu.SetMarkerColor(ROOT.kOrange)
    g_pro.SetMarkerColor(ROOT.kGreen)
    c.Write()

    


if __name__ == "__main__":
    check_smoothness()