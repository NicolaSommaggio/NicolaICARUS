import ROOT
import uproot
import matplotlib.pyplot as plt
import numpy as np
import awkward as ak
import pandas
from matplotlib.colors import LogNorm

def root_hist_to_pyhist(th1d, colore, llabel, rebin=False, rebfactor=1):
    if rebin : th1d.Rebin(rebfactor)
    n_bins = th1d.GetNbinsX()
    x_edges = [th1d.GetBinLowEdge(i+1) for i in range(n_bins)]
    x_edges.append(th1d.GetBinLowEdge(n_bins) + th1d.GetBinWidth(n_bins))
    bin_contents = [th1d.GetBinContent(i+1) for i in range(n_bins)]
    bin_errors = [th1d.GetBinError(i+1) for i in range(n_bins)]
    bin_centers =[th1d.GetBinCenter(i+1) for i in range(n_bins)]
    x_bin_errors=[th1d.GetBinWidth(i+1)/2 for i in range(n_bins)]

    
    # Creazione step per linee orizzontali
    x_step = np.repeat(x_edges, 2)[1:-1]
    y_step = np.repeat(bin_contents, 2)
    
    # Aggiunta linee verticali agli estremi
    x_step = np.insert(x_step, 0, x_edges[0])
    y_step = np.insert(y_step, 0, 0)  # parte dal basso
    x_step = np.append(x_step, x_edges[-1])
    y_step = np.append(y_step, 0)     # torna a 0 alla fine
    
    plt.step(x_step, y_step, where='mid', linewidth=2, color=colore, label=llabel)
    
    #plt.errorbar(bin_centers,bin_contents,xerr=x_bin_errors,yerr=bin_errors, fmt='o', color=colore, label=llabel, markersize=1)

    plt.ylim(0, max(bin_contents)*1.1)
    #plt.savefig('histo_prova_bar.pdf', format='pdf', bbox_inches='tight')



def root_hist_to_pyhist_errors(th1d, colore, llabel, rebin=False, rebfactor=1):
    if rebin: th1d.Rebin(rebfactor)
    n_bins = th1d.GetNbinsX()
    x_edges = [th1d.GetBinLowEdge(i+1) for i in range(n_bins)]
    x_edges.append(th1d.GetBinLowEdge(n_bins) + th1d.GetBinWidth(n_bins))
    bin_contents = [th1d.GetBinContent(i+1) for i in range(n_bins)]
    bin_errors = [th1d.GetBinError(i+1) for i in range(n_bins)]
    bin_centers =[th1d.GetBinCenter(i+1) for i in range(n_bins)]
    x_bin_errors=[th1d.GetBinWidth(i+1)/2 for i in range(n_bins)]

    
    # Creazione step per linee orizzontali
    x_step = np.repeat(x_edges, 2)[1:-1]
    y_step = np.repeat(bin_contents, 2)
    
    # Aggiunta linee verticali agli estremi
    x_step = np.insert(x_step, 0, x_edges[0])
    y_step = np.insert(y_step, 0, 0)  # parte dal basso
    x_step = np.append(x_step, x_edges[-1])
    y_step = np.append(y_step, 0)     # torna a 0 alla fine
    
    plt.step(x_step, y_step, where='mid', linewidth=2, color=colore, label=llabel)

    plt.errorbar(
        bin_centers,
        bin_contents,
        yerr=bin_errors,
        fmt='',           
        markersize=0,
        color=colore,
        ecolor=colore,
        capsize=3,
        alpha=0.7,
        linestyle='none',
        label=None
    )

    #plt.grid(color='lightgray', linestyle='--', linewidth=0.5, alpha=0.7, axis='x')
    
    #plt.errorbar(bin_centers,bin_contents,xerr=x_bin_errors,yerr=bin_errors, fmt='o', color=colore, label=llabel, markersize=1)

    #plt.ylim(0, max(bin_contents)*1.1)
    plt.ylim(1e-4, max(bin_contents)*1.1)
    #plt.savefig('histo_prova_bar.pdf', format='pdf', bbox_inches='tight')

def root_hist_to_pyerrorbar(th1d, colore, llabel, rebin=False, rebfactor=1, marksize=1):
    if rebin: th1d.Rebin(rebfactor)
    n_bins = th1d.GetNbinsX()
    x_edges = [th1d.GetBinLowEdge(i+1) for i in range(n_bins)]
    x_edges.append(th1d.GetBinLowEdge(n_bins) + th1d.GetBinWidth(n_bins))
    bin_contents = [th1d.GetBinContent(i+1) for i in range(n_bins)]
    bin_errors = [th1d.GetBinError(i+1) for i in range(n_bins)]
    bin_centers =[th1d.GetBinCenter(i+1) for i in range(n_bins)]
    x_bin_errors=[th1d.GetBinWidth(i+1)/2 for i in range(n_bins)]
    
    plt.errorbar(bin_centers,bin_contents,yerr=bin_errors, fmt='', color=colore, label=llabel, markersize=marksize, alpha=0.7,capsize=3,linestyle='none')
    plt.scatter(bin_centers,bin_contents,color=colore,s=4,label=None)

def root_hist_to_pyerrorbar_mod(th1d, colore, llabel, rebin=False, rebfactor=1, marksize=1):
    if rebin:
        th1d.Rebin(rebfactor)

    n_bins = th1d.GetNbinsX()

    # usa un bin sì e uno no: i = 0, 2, 4, ...
    selected_bins = range(1, n_bins+1, 2)

    bin_centers  = [th1d.GetBinCenter(i)   for i in selected_bins]
    bin_contents = [th1d.GetBinContent(i)  for i in selected_bins]
    bin_errors   = [th1d.GetBinError(i)    for i in selected_bins]
    x_bin_errors = [th1d.GetBinWidth(i)/2  for i in selected_bins]

    plt.errorbar(
        bin_centers,
        bin_contents,
        yerr=bin_errors,
        fmt='',
        color=colore,
        label=llabel,
        alpha=0.7,
        capsize=3,
        linestyle='none'
    )

    plt.scatter(bin_centers, bin_contents, color=colore, s=4)


def tgraph_to_matplotlib(graph, llabel,ccolor,a=1):
    # Estrazione punti e errori
    n_points = graph.GetN()
    x = np.array([graph.GetPointX(i) for i in range(n_points)], dtype=float)
    y = np.array([graph.GetPointY(i) for i in range(n_points)], dtype=float)

    plt.scatter(x, y, s=2, color=ccolor,alpha=a)
    plt.plot(x,y,label=llabel,color=ccolor, lw=2,alpha=a)

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root','UPDATE')
    graph = file.Get('rocCurve_chi2_coll')
    n_points = graph.GetN()
    for point in range(n_points):
        if point==316 or point==320 or point==324:
            print(graph.GetPointX(point),graph.GetPointY(point),graph.GetPointY(0))
            graph.SetPoint(point,graph.GetPointX(point),graph.GetPointY(0))
    graph.Write('rocCurve_chi2_coll',ROOT.TObject.kOverwrite)

    

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root')
    muon_lkl_coll = file.Get("muon_likelihood_ratio_coll")
    muon_lkl_coll_ind1 = file.Get("muon_likelihood_ratio_coll_ind1")
    muon_lkl_coll_ind1_ind2 = file.Get("muon_likelihood_ratio_coll_ind1_ind2")
    proton_lkl_coll = file.Get("proton_likelihood_ratio_coll")
    proton_lkl_coll_ind1 = file.Get("proton_likelihood_ratio_coll_ind1")
    proton_lkl_coll_ind1_ind2 = file.Get("proton_likelihood_ratio_coll_ind1_ind2")
    #plt.plot([], [], ' ', label="Stats errors only")
    root_hist_to_pyhist(muon_lkl_coll, 'maroon', "muons LR COLL",False,2)
    root_hist_to_pyhist(muon_lkl_coll_ind1, 'red', "muons LR COLL+IND1",False,2)
    root_hist_to_pyhist(muon_lkl_coll_ind1_ind2, 'orange', "muons LR COLL+IND1+IND2",False,2)
    root_hist_to_pyhist(proton_lkl_coll, 'darkblue', "protons LR COLL",False,2)
    root_hist_to_pyhist(proton_lkl_coll_ind1, 'dodgerblue', "protons LR COLL+IND1",False,2)
    root_hist_to_pyhist(proton_lkl_coll_ind1_ind2, 'cyan', "protons LR COLL+IND1+IND2",False,2)
    
    plt.xlabel('likelihood ratio (LR)', fontsize=18)
    plt.ylabel('entries (area normalized)', fontsize=18)
    plt.title("MC 2D DECONVOLUTION", fontsize=18)

    leg = plt.legend(loc='upper center', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(-1,1.,0.1), fontsize=14, rotation=60)
    plt.yticks(fontsize=14)
    plt.ylim(1e-4,20)
    plt.yscale('log')
    plt.savefig('grafici/new_likelihood_ratios_log.svg', format='svg', bbox_inches='tight')
    plt.close()

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root')
    muon_lkl_coll = file.Get("muon_likelihood_ratio_coll")
    muon_lkl_coll_ind1 = file.Get("muon_likelihood_ratio_coll_ind1")
    muon_lkl_coll_ind1_ind2 = file.Get("muon_likelihood_ratio_coll_ind1_ind2")
    proton_lkl_coll = file.Get("proton_likelihood_ratio_coll")
    proton_lkl_coll_ind1 = file.Get("proton_likelihood_ratio_coll_ind1")
    proton_lkl_coll_ind1_ind2 = file.Get("proton_likelihood_ratio_coll_ind1_ind2")
    #plt.plot([], [], ' ', label="Stats errors only")
    root_hist_to_pyhist_errors(muon_lkl_coll, 'maroon', "muons LR COLL",False,2)
    root_hist_to_pyhist_errors(muon_lkl_coll_ind1, 'red', "muons LR COLL+IND1",False,2)
    root_hist_to_pyhist_errors(muon_lkl_coll_ind1_ind2, 'orange', "muons LR COLL+IND1+IND2",False,2)
    root_hist_to_pyhist_errors(proton_lkl_coll, 'darkblue', "protons LR COLL",False,2)
    root_hist_to_pyhist_errors(proton_lkl_coll_ind1, 'dodgerblue', "protons LR COLL+IND1",False,2)
    root_hist_to_pyhist_errors(proton_lkl_coll_ind1_ind2, 'cyan', "protons LR COLL+IND1+IND2",False,2)
    
    plt.xlabel('likelihood ratio (LR)', fontsize=18)
    plt.ylabel('entries (area normalized)', fontsize=18)
    plt.title("MC 2D DECONVOLUTION", fontsize=18)

    leg = plt.legend(loc='upper center', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(-1,1.,0.1), fontsize=14, rotation=60)
    plt.yticks(fontsize=14)
    plt.ylim(0,1.015)
    plt.savefig('grafici/new_likelihood_ratios.svg', format='svg', bbox_inches='tight')
    plt.close()

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root')
    roc_coll = file.Get("rocCurve_coll")
    roc_coll_ind1 = file.Get("rocCurve_coll_ind1")
    roc_coll_ind1_ind2 = file.Get("rocCurve_coll_ind1_ind2")

    roc_chi2_coll = file.Get("rocCurve_chi2_coll")
    roc_chi2_coll_ind1 = file.Get("rocCurve_chi2_coll_ind1")
    roc_chi2_coll_ind1_ind2 = file.Get("rocCurve_chi2_coll_ind1_ind2")

    tgraph_to_matplotlib(roc_chi2_coll,r"$\chi^2$ COLL","limegreen",0.7)
    tgraph_to_matplotlib(roc_chi2_coll_ind1, r"$\chi^2$ COLL+IND1","deepskyblue",0.7)
    tgraph_to_matplotlib(roc_chi2_coll_ind1_ind2, r"$\chi^2$ COLL+IND1+IND2", "orange",0.7)
    tgraph_to_matplotlib(roc_coll,"LR COLL","darkgreen",0.7)
    tgraph_to_matplotlib(roc_coll_ind1, "LR COLL+IND1","darkblue",0.7)
    tgraph_to_matplotlib(roc_coll_ind1_ind2, "LR COLL+IND1+IND2", "red",0.7)

    plt.xlabel(r'background rejection eff. ($1-\epsilon_{back.}$)', fontsize=18)
    plt.ylabel(r'signal eff. ($\epsilon_{sig.}$)', fontsize=18)
    plt.title("MC 2D DECONVOLUTION ROC Curves", fontsize=18)

    leg = plt.legend(loc='lower left', fontsize=14)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0.,1.,0.1), fontsize=14, rotation=60)
    plt.yticks(fontsize=14)
    plt.xlim(0.,1.01)
    plt.ylim(0.98,1.005)
    #plt.yscale('log')
    plt.savefig('grafici/rocCurves_total_zoom.svg', format='svg', bbox_inches='tight')
    plt.close()

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root')
    roc_chi2_coll = file.Get("rocCurve_chi2_coll")
    roc_chi2_coll_ind1 = file.Get("rocCurve_chi2_coll_ind1")
    roc_chi2_coll_ind1_ind2 = file.Get("rocCurve_chi2_coll_ind1_ind2")

    tgraph_to_matplotlib(roc_chi2_coll,r"$\chi^2$ COLL","green",0.7)
    tgraph_to_matplotlib(roc_chi2_coll_ind1, r"$\chi^2$ COLL+IND1","deepskyblue",0.7)
    tgraph_to_matplotlib(roc_chi2_coll_ind1_ind2, r"$\chi^2$ COLL+IND1+IND2", "orangered",0.7)

    plt.xlabel(r'background rejection eff. ($1-\epsilon_{back.}$)', fontsize=18)
    plt.ylabel(r'signal eff. ($\epsilon_{sig.}$)', fontsize=18)
    plt.title("MC 2D DECONVOLUTION ROC Curves", fontsize=18)

    leg = plt.legend(loc='center', fontsize=14)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,1.,0.1), fontsize=14, rotation=60)
    plt.yticks(fontsize=14)
    plt.ylim(0.,1.02)
    #plt.yscale('log')
    plt.savefig('grafici/rocCurves_chi2.svg', format='svg', bbox_inches='tight')
    plt.close()

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root')
    roc_chi2_coll = file.Get("rocCurve_chi2_coll")
    roc_chi2_coll_ind1 = file.Get("rocCurve_chi2_coll_ind1")
    roc_chi2_coll_ind1_ind2 = file.Get("rocCurve_chi2_coll_ind1_ind2")

    tgraph_to_matplotlib(roc_chi2_coll,r"$\chi^2$ COLL","green",0.7)
    tgraph_to_matplotlib(roc_chi2_coll_ind1, r"$\chi^2$ COLL+IND1","deepskyblue",0.7)
    tgraph_to_matplotlib(roc_chi2_coll_ind1_ind2, r"$\chi^2$ COLL+IND1+IND2", "orangered",0.7)

    plt.xlabel(r'background rejection eff. ($1-\epsilon_{back.}$)', fontsize=18)
    plt.ylabel(r'signal eff. ($\epsilon_{sig.}$)', fontsize=18)
    plt.title("MC 2D DECONVOLUTION ROC Curves", fontsize=18)

    leg = plt.legend(loc='center left', fontsize=14)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,1.,0.1), fontsize=14, rotation=60)
    plt.yticks(fontsize=14)
    plt.ylim(0.98,1.002)
    #plt.yscale('log')
    plt.savefig('grafici/rocCurves_chi2_zoom.svg', format='svg', bbox_inches='tight')
    plt.close()

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root')
    chi_mu_coll_asmu = file.Get("muon_hchi2_coll_asmu")
    chi_mu_coll_aspro = file.Get("muon_hchi2_coll_aspro")
    chi_mu_coll_ind1_asmu = file.Get("muon_hchi2_coll_ind1_asmu")
    chi_mu_coll_ind1_aspro = file.Get("muon_hchi2_coll_ind1_aspro")
    chi_mu_coll_ind1_ind2_asmu = file.Get("muon_hchi2_coll_ind1_ind2_asmu")
    chi_mu_coll_ind1_ind2_aspro = file.Get("muon_hchi2_coll_ind1_ind2_aspro")
    #plt.plot([], [], ' ', label="Stats errors only")
    root_hist_to_pyhist(chi_mu_coll_asmu, 'maroon', r"muons $\chi^2$ as $\mu$ COLL",True,2)
    root_hist_to_pyhist(chi_mu_coll_ind1_asmu, 'red', r"muons $\chi^2$ as $\mu$ COLL+IND1",True,2)
    root_hist_to_pyhist(chi_mu_coll_ind1_ind2_asmu, 'orange', r"muons $\chi^2$ as $\mu$ COLL+IND1+IND2",True,2)
    root_hist_to_pyhist(chi_mu_coll_aspro, 'darkblue', r"muons $\chi^2$ as $p$ COLL",True,2)
    root_hist_to_pyhist(chi_mu_coll_ind1_aspro, 'dodgerblue', r"muons $\chi^2$ as $p$ COLL+IND1",True,2)
    root_hist_to_pyhist(chi_mu_coll_ind1_ind2_aspro, 'cyan', r"muons $\chi^2$ as $p$ COLL+IND1+IND2",True,2)
    
    plt.xlabel(r'$\chi^2$', fontsize=18)
    plt.ylabel('entries (area normalized)', fontsize=18)
    plt.title(r"MUONS $\chi^2$, MC 2D DECONVOLUTION", fontsize=18)

    leg = plt.legend(loc='upper center', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,200.,10), fontsize=14, rotation=60)
    plt.xlim(0,200)
    plt.ylim(0,0.45)
    plt.yticks(fontsize=14)
    plt.savefig('grafici/chi2_adding_planes.svg', format='svg', bbox_inches='tight')
    plt.close()

if False:
    file = ROOT.TFile.Open('likelihood_ratios.root')
    chi_pro_coll_asmu = file.Get("proton_hchi2_coll_asmu")
    chi_pro_coll_aspro = file.Get("proton_hchi2_coll_aspro")
    chi_pro_coll_ind1_asmu = file.Get("proton_hchi2_coll_ind1_asmu")
    chi_pro_coll_ind1_aspro = file.Get("proton_hchi2_coll_ind1_aspro")
    chi_pro_coll_ind1_ind2_asmu = file.Get("proton_hchi2_coll_ind1_ind2_asmu")
    chi_pro_coll_ind1_ind2_aspro = file.Get("proton_hchi2_coll_ind1_ind2_aspro")
    #plt.plot([], [], ' ', label="Stats errors only")
    root_hist_to_pyhist(chi_pro_coll_asmu, 'maroon', r"protons $\chi^2$ as $\mu$ COLL",True,2)
    root_hist_to_pyhist(chi_pro_coll_ind1_asmu, 'red', r"protons $\chi^2$ as $\mu$ COLL+IND1",True,2)
    root_hist_to_pyhist(chi_pro_coll_ind1_ind2_asmu, 'orange', r"protons $\chi^2$ as $\mu$ COLL+IND1+IND2",True,2)
    root_hist_to_pyhist(chi_pro_coll_aspro, 'darkblue', r"protons $\chi^2$ as $p$ COLL",True,2)
    root_hist_to_pyhist(chi_pro_coll_ind1_aspro, 'dodgerblue', r"protons $\chi^2$ as $p$ COLL+IND1",True,2)
    root_hist_to_pyhist(chi_pro_coll_ind1_ind2_aspro, 'cyan', r"protons $\chi^2$ as $p$ COLL+IND1+IND2",True,2)
    
    plt.xlabel(r'$\chi^2$', fontsize=18)
    plt.ylabel('entries (area normalized)', fontsize=18)
    plt.title(r"PROTONS $\chi^2$, MC 2D DECONVOLUTION", fontsize=18)

    leg = plt.legend(loc='upper right', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,100.,10), fontsize=14, rotation=60)
    plt.yticks(fontsize=14)
    plt.xlim(0,100)
    plt.ylim(0,0.3)
    plt.savefig('grafici/chi2_adding_planes_protons.svg', format='svg', bbox_inches='tight')
    plt.close()


#PROFILES OF dEdx VS RR IN THE THREE DIFF PLANES
if False:
    file = uproot.open("../datafiles/mc2d.root")
    mutree = file['DATAtreeMU']
    protree = file['DATAtreePRO']
    branches_mu = mutree.arrays()
    branches_pro = protree.arrays()

    #print(mutree.keys())
    dEdx_mu = np.array(ak.flatten(branches_mu['dE_mu']))
    dEdx_pro = np.array(ak.flatten(branches_pro['dE_pro']))
    rr_mu = np.array(ak.flatten(branches_mu['rr_mu']))
    rr_pro = np.array(ak.flatten(branches_pro['rr_pro']))

    dEdx_mu_ind1 = np.array(ak.flatten(branches_mu['dE_mu_ind1']))
    dEdx_mu_ind2 = np.array(ak.flatten(branches_mu['dE_mu_ind2']))
    rr_mu_ind1 = np.array(ak.flatten(branches_mu['rr_mu_ind1']))
    rr_mu_ind2 = np.array(ak.flatten(branches_mu['rr_mu_ind2']))
    dEdx_pro_ind1 = np.array(ak.flatten(branches_pro['dE_pro_ind1']))
    dEdx_pro_ind2 = np.array(ak.flatten(branches_pro['dE_pro_ind2']))
    rr_pro_ind1 = np.array(ak.flatten(branches_pro['rr_pro_ind1']))
    rr_pro_ind2 = np.array(ak.flatten(branches_pro['rr_pro_ind2']))

    print(len(ak.to_list(branches_mu['dE_mu'])), "tracce di muone")
    print(len(ak.to_list(branches_pro['dE_pro'])), "tracce di protone")

    file_ref = ROOT.TFile.Open('../dEdx/THdedx.root', 'UPDATE')
    mu_curve = file_ref.Get('dedx_range_mu')
    pro_curve = file_ref.Get('dedx_range_pro')
    #for i in range(1,mu_curve.GetNbinsX()+1):
        #print(mu_curve.GetNbinsX(),mu_curve.GetBinCenter(i), mu_curve.GetBinWidth(i)) 
    profile_pro_coll = ROOT.TProfile('profile_pro_coll','',325,0,26,0,40)
    profile_pro_ind1 = ROOT.TProfile('profile_pro_ind1','',325,0,26,0,40)
    profile_pro_ind2 = ROOT.TProfile('profile_pro_ind2','',325,0,26,0,40)
    #profile_pro_coll = ROOT.TProfile('profile_pro_coll','',100,0,25,0,40)
    #profile_pro_ind1 = ROOT.TProfile('profile_pro_ind1','',100,0,25,0,40)
    #profile_pro_ind2 = ROOT.TProfile('profile_pro_ind2','',100,0,25,0,40)

    for track in range(len(branches_pro['dE_pro'])):
        for hit_coll in range(len(branches_pro['dE_pro'][track])):
            profile_pro_coll.Fill(branches_pro['rr_pro'][track][hit_coll], branches_pro['dE_pro'][track][hit_coll]) 
        for hit_ind1 in range(len(branches_pro['dE_pro_ind1'][track])):
            profile_pro_ind1.Fill(branches_pro['rr_pro_ind1'][track][hit_ind1],branches_pro['dE_pro_ind1'][track][hit_ind1])
        for hit_ind2 in range(len(branches_pro['dE_pro_ind2'][track])):
            profile_pro_ind2.Fill(branches_pro['rr_pro_ind2'][track][hit_ind2],branches_pro['dE_pro_ind2'][track][hit_ind2])

    profile_pro_coll.Write('profile_pro_coll', ROOT.TObject.kOverwrite)
    profile_pro_ind1.Write('profile_pro_ind1', ROOT.TObject.kOverwrite)
    profile_pro_ind2.Write('profile_pro_ind2', ROOT.TObject.kOverwrite)

    rr_coll = []
    rr_ind1 = []
    rr_ind2 = []
    rr_reference = []
    err_coll = []
    err_ind1 = []
    err_ind2 = []
    err_reference = []

    
    for i in range(1,profile_pro_coll.GetNbinsX()+1):
        rr_coll.append(profile_pro_coll.GetBinCenter(i))
        err_coll.append(profile_pro_coll.GetBinError(i))
        rr_ind1.append(profile_pro_ind1.GetBinCenter(i))
        err_ind1.append(profile_pro_ind1.GetBinError(i))
        rr_ind2.append(profile_pro_ind2.GetBinCenter(i))
        err_ind2.append(profile_pro_ind2.GetBinError(i))
    for i in range(1,pro_curve.GetNbinsX()+1):
        rr_reference.append(pro_curve.GetBinCenter(i))
        err_reference.append(pro_curve.GetBinError(i))
    

    #root_hist_to_pyerrorbar_mod(pro_curve, 'red', 'Reference Curve', False,4)
    #root_hist_to_pyerrorbar(profile_pro_coll, 'cornflowerblue', 'Profile COLL ICARUS-MC 2D',False,4)
    #root_hist_to_pyerrorbar(profile_pro_ind1, 'orange','Profile IND1 ICARUS-MC 2D' ,False,4)
    #root_hist_to_pyerrorbar(profile_pro_ind2, 'green', 'Profile IND2 ICARUS-MC 2D',False,4)

    plt.scatter(rr_reference, err_reference, color='red', label='Reference curve', s=4)
    plt.scatter(rr_coll, err_coll, color='cornflowerblue', label='Profile COLL ICARUS-MC 2D', s=4)
    plt.scatter(rr_ind1, err_ind1, color='orange', label='Profile IND1 ICARUS-MC 2D', s=4)
    plt.scatter(rr_ind2, err_ind2, color='green', label='Profile IND2 ICARUS-MC 2D', s=4)
    
    
    
    plt.xlabel('Residual Range [cm]', fontsize=18)
    #plt.ylabel('bin average dE/dx [MeV/cm]', fontsize=18)
    plt.ylabel('bin dE/dx std err of the mean [MeV/cm]', fontsize=18)
    plt.title("PROTONS, MC 2D DECONVOLUTION", fontsize=18)

    leg = plt.legend(loc='upper right', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,25.,1), fontsize=14, rotation=60)
    plt.xlim(0,25)
    plt.ylim(0,3)
    plt.yticks(fontsize=14)
    plt.savefig('grafici/errori.pdf', format='pdf', bbox_inches='tight')
    plt.close()



#NUMBER OF HITS PER TRACK IN THE THREE DIFF PLANES
if False:
    n_hits_mu_coll=ROOT.TH1D("n_hits_mu_coll","",200,0,200)
    n_hits_mu_ind1=ROOT.TH1D("n_hits_mu_ind1","",200,0,200)
    n_hits_mu_ind2=ROOT.TH1D("n_hits_mu_ind2","",200,0,200)

    n_hits_pro_coll=ROOT.TH1D("n_hits_pro_coll","",200,0,200)
    n_hits_pro_ind1=ROOT.TH1D("n_hits_pro_ind1","",200,0,200)
    n_hits_pro_ind2=ROOT.TH1D("n_hits_pro_ind2","",200,0,200)

    #print(branches_pro['rr_pro'])

    for i in range(len(branches_mu['dE_mu'])):
        n_hits_mu_coll.Fill(len(branches_mu['dE_mu'][i]))
        n_hits_mu_ind1.Fill(len(branches_mu['dE_mu_ind1'][i]))
        n_hits_mu_ind2.Fill(len(branches_mu['dE_mu_ind2'][i]))

    for i in range(len(branches_pro['dE_pro'])):
        n_hits_pro_coll.Fill(len(branches_pro['dE_pro'][i]))
        n_hits_pro_ind1.Fill(len(branches_pro['dE_pro_ind1'][i]))
        n_hits_pro_ind2.Fill(len(branches_pro['dE_pro_ind2'][i]))


    root_hist_to_pyhist_errors(n_hits_mu_coll, 'cornflowerblue','#hits per track COLL' ,True,4)
    root_hist_to_pyhist_errors(n_hits_mu_ind1, 'orange','#hits per track IND1' ,True,4)
    root_hist_to_pyhist_errors(n_hits_mu_ind2, 'green','#hits per track IND2' ,True,4)
    
    plt.xlabel('number of hits per track', fontsize=18)
    plt.ylabel('entries (area normalized)', fontsize=18)
    plt.title(r"MUONS true 1$\mu$Np - MC 2D DECONVOLUTION", fontsize=18)

    leg = plt.legend(loc='upper center', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,200.,20), fontsize=14, rotation=60)
    plt.xlim(0,140)
    plt.ylim(0,100)
    plt.yticks(fontsize=14)
    plt.savefig('grafici/number_of_hits_muon.pdf', format='pdf', bbox_inches='tight')
    plt.close()

    root_hist_to_pyhist_errors(n_hits_pro_coll, 'cornflowerblue','#hits per track COLL' ,True,4)
    root_hist_to_pyhist_errors(n_hits_pro_ind1, 'orange','#hits per track IND1' ,True,4)
    root_hist_to_pyhist_errors(n_hits_pro_ind2, 'green','#hits per track IND2' ,True,4)
    
    plt.xlabel('number of hits per track', fontsize=18)
    plt.ylabel('entries (area normalized)', fontsize=18)
    plt.title(r"PROTONS true 1$\mu$Np - MC 2D DECONVOLUTION", fontsize=18)

    leg = plt.legend(loc='upper right', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,200.,20), fontsize=14, rotation=60)
    plt.xlim(0,140)
    plt.ylim(0, 210)
    plt.yticks(fontsize=14)
    plt.savefig('grafici/number_of_hits_protons.pdf', format='pdf', bbox_inches='tight')
    plt.close()



#PLOT dEdx vs RR
if False:

    colore_titolo = (34,73,120)
    colore_titolo_mpl = tuple(c / 255 for c in colore_titolo)

    file = uproot.open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2d_general_contained.root")
    mutree = file['DATAtreeMU']
    protree = file['DATAtreePRO']
    pitree = file['DATAtreePI']

    branches_mu = mutree.arrays()
    branches_pro = protree.arrays()
    branches_pi = pitree.arrays()

    print(len(ak.to_list(branches_mu['dE_mu'])), "tracce di muone")
    print(len(ak.to_list(branches_pro['dE_pro'])), "tracce di protone")
    print(len(ak.to_list(branches_pi['dE_pi'])), "tracce di pione")

    #PROTONS
    dEdx_pro = np.array(ak.flatten(branches_pro['dE_pro']))
    rr_pro = np.array(ak.flatten(branches_pro['rr_pro']))
    dEdx_pro_ind1 = np.array(ak.flatten(branches_pro['dE_pro_ind1']))
    dEdx_pro_ind2 = np.array(ak.flatten(branches_pro['dE_pro_ind2']))
    rr_pro_ind1 = np.array(ak.flatten(branches_pro['rr_pro_ind1']))
    rr_pro_ind2 = np.array(ak.flatten(branches_pro['rr_pro_ind2']))

    plt.hist2d(rr_pro, dEdx_pro, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("PROTONS MC 2D OVERLAYS - COLLECTION PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/pro_dEdx_rr_mc_coll.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    plt.hist2d(rr_pro_ind1, dEdx_pro_ind1, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("PROTONS MC 2D OVERLAYS - INDUCTION 1 PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/pro_dEdx_rr_mc_ind1.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    plt.hist2d(rr_pro_ind2, dEdx_pro_ind2, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("PROTONS MC 2D OVERLAYS - INDUCTION 2 PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/pro_dEdx_rr_mc_ind2.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    #MUONS
    dEdx_mu = np.array(ak.flatten(branches_mu['dE_mu']))
    rr_mu = np.array(ak.flatten(branches_mu['rr_mu']))
    dEdx_mu_ind1 = np.array(ak.flatten(branches_mu['dE_mu_ind1']))
    dEdx_mu_ind2 = np.array(ak.flatten(branches_mu['dE_mu_ind2']))
    rr_mu_ind1 = np.array(ak.flatten(branches_mu['rr_mu_ind1']))
    rr_mu_ind2 = np.array(ak.flatten(branches_mu['rr_mu_ind2']))

    plt.hist2d(rr_mu, dEdx_mu, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("MUONS MC 2D OVERLAYS - COLLECTION PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/mu_dEdx_rr_mc_coll.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    plt.hist2d(rr_mu_ind1, dEdx_mu_ind1, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("MUONS MC 2D OVERLAYS - INDUCTION 1 PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/mu_dEdx_rr_mc_ind1.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    plt.hist2d(rr_mu_ind2, dEdx_mu_ind2, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("MUONS MC 2D OVERLAYS - INDUCTION 2 PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/mu_dEdx_rr_mc_ind2.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    #PIONS
    #mask = branches_pi['end_process_pi'] == 3
    dEdx = np.array(ak.flatten(branches_pi['dE_pi']))
    rr = np.array(ak.flatten(branches_pi['rr_pi']))
    dEdx_ind1 = np.array(ak.flatten(branches_pi['dE_pi_ind1']))
    dEdx_ind2 = np.array(ak.flatten(branches_pi['dE_pi_ind2']))
    rr_ind1 = np.array(ak.flatten(branches_pi['rr_pi_ind1']))
    rr_ind2 = np.array(ak.flatten(branches_pi['rr_pi_ind2']))

    plt.hist2d(rr, dEdx, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - COLLECTION PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/pi_dEdx_rr_mc_coll.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    plt.hist2d(rr_ind1, dEdx_ind1, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - INDUCTION 1 PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/pi_dEdx_rr_mc_ind1.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    plt.hist2d(rr_ind2, dEdx_ind2, bins=(100,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm())
    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - INDUCTION 2 PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/pi_dEdx_rr_mc_ind2.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()




def copy_dir(src, dst):
    for key in src.GetListOfKeys():
        obj = key.ReadObj()

        if obj.InheritsFrom("TDirectory"):
            newdir = dst.mkdir(obj.GetName())
            copy_dir(obj, newdir)
        else:
            dst.cd()
            obj.Write()

if False:
    f = ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/newHISTOdistro2dprova.root", "UPDATE")

    src = f.Get("muon")
    dst = f.mkdir("pion")

    copy_dir(src, dst)

    f.Close()

if False:
    f_src = ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2d_general.root", "READ")
    f_dst = ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2dprova.root", "UPDATE")

    t_src = f_src.Get("DATAtreePI")

    f_dst.cd()

    t_dst = t_src.CloneTree(-1)
    t_dst.SetName("DATAtreePI")

    t_dst.Write()

    f_dst.Close()
    f_src.Close()

if False:
    file_origine = ROOT.TFile.Open("../dEdx/THdedx.root")
    file_destinazione = ROOT.TFile("../dEdx/RefCurvesChi2.root", "RECREATE")
    mu = file_origine.Get("dedx_range_mu")
    pro = file_origine.Get("dedx_range_pro")
    pi = file_origine.Get("dedx_range_pi")
    ka = file_origine.Get("dedx_range_ka")
    file_destinazione.cd()
    mu.Write()
    pro.Write()
    pi.Write()
    ka.Write()

#PLOT ROLLING MEDIAN
if False:
    file = uproot.open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2d_general_contained.root")

    #MUONS
    mu_tree = file['DATAtreeMU']
    mu_branches = mu_tree.arrays()

    #mask = (mu_branches["end_process_mu"]==3)
    rr_rm_mu = np.array(ak.flatten(mu_branches["rr_rm_mu"]))
    rm_mu = np.array(ak.flatten(mu_branches["rm_mu"]))

    colore_titolo = (34,73,120)
    colore_titolo_mpl = tuple(c / 255 for c in colore_titolo)
    
    h=plt.hist2d(rr_rm_mu, rm_mu, bins=(300,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm(), cmin=1)
    plt.xlabel('Residual Range [cm]', fontsize=17)
    plt.ylabel('dE/dx Rolling Median [MeV/cm]', fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\mu$ MC 2D OVERLAYS - COLLECTION PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/rolling_median_mu.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    #PROTONS
    pro_tree = file['DATAtreePRO']
    pro_branches = pro_tree.arrays()

    #mask = (pro_branches["end_process_pro"]==3)
    rr_rm_pro = np.array(ak.flatten(pro_branches["rr_rm_pro"]))
    rm_pro = np.array(ak.flatten(pro_branches["rm_pro"]))
    
    h=plt.hist2d(rr_rm_pro, rm_pro, bins=(300,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm(), cmin=1)
    plt.xlabel('Residual Range [cm]', fontsize=17)
    plt.ylabel('dE/dx Rolling Median [MeV/cm]', fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$p$ MC 2D OVERLAYS - COLLECTION PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/rolling_median_pro.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    #PIONS
    pi_tree = file['DATAtreePI']
    pi_branches = pi_tree.arrays()

    #mask = (pro_branches["end_process_pro"]==3)
    rr_rm_pi = np.array(ak.flatten(pi_branches["rr_rm_pi"]))
    rm_pi = np.array(ak.flatten(pi_branches["rm_pi"]))
    
    h=plt.hist2d(rr_rm_pi, rm_pi, bins=(300,300), range=[(0,30),(0,30)], cmap='viridis', norm=LogNorm(), cmin=1)
    plt.xlabel('Residual Range [cm]', fontsize=17)
    plt.ylabel('dE/dx Rolling Median [MeV/cm]', fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - COLLECTION PLANE", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici'
    '/rolling_median_pi.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

if False:
    reference = ROOT.TFile.Open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/dEdx/THdedx.root","READ")
    dedx_range_mu = reference.Get("dedx_range_mu")
    scatter_dedx = []
    scatter_rr = []
    for i in range(1,dedx_range_mu.GetNbinsX()+1):
        if dedx_range_mu.GetBinContent(i)!=0:
            scatter_rr.append(dedx_range_mu.GetBinCenter(i))
            scatter_dedx.append(dedx_range_mu.GetBinContent(i))
    scatter_dedx = np.array(scatter_dedx)
    scatter_rr = np.array(scatter_rr)
    mip_dedx =np.array([2,2])
    mip_rr = np.array([0,25])



    file = uproot.open("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/datafiles/mc2d_general_contained.root")
    mu_tree = file['DATAtreeMU']
    mu_branches = mu_tree.arrays()

    #increasing
    #run = 9960
    #evt = 52964
    #chi2_mip = 2.78
    #chi2_mip_dedx= 6.56
    #chi2 = 0.36
    #chi2_dedx = 22.80
    #range(0,16)

    #mip
    #run = 9388
    #evt = 65257
    #chi2_mip = 0.06
    #chi2_mip_dedx= 0.36
    #chi2 = 1.89
    #chi2_dedx = 5.03
    #range(0,9)

    #MICHEL
    #run = 9717
    #evt = 10481
    #chi2_mip = 2.22
    #chi2_mip_dedx= 5.90
    #chi2 = 1.79
    #chi2_dedx = 12.91
    #range(0,20)

    #PION ClASSIFICATO MIP CON PROTON
    #run = 9448
    #evt = 57075
    #chi2_mip = 0.05
    #chi2_mip_dedx= 1.89
    #chi2 = 2.38
    #chi2_dedx = 5.48

    #PION CLASSIFICATO STOPPING CON PROTON
    #run = 9384
    #evt = 82560
    #chi2_mip = 3.59
    #chi2_mip_dedx= 47.03
    #chi2 = 0.17
    #chi2_dedx = 25.48

    #MUONE ESEMPIO
    #run = 9594
    #evt = 70308

    run = 9942
    evt = 2734

    mask = (mu_branches["run"] == run) & (mu_branches["evt"] == evt)
    bestplane = mu_branches["bestplane_mu"][mask]
    str_bestplane = ""
    if bestplane==0 : 
        print("the best plane is INDUCTION 1")
        str_bestplane="INDUCTION 1"
    if bestplane==1 : 
        print("the best plane is INDUCTION 2")
        str_bestplane="INDUCTION 2"
    if bestplane==2 : 
        print("the best plane is COLLECTION")
        str_bestplane="COLLECTION"
    
    rr = np.array(ak.flatten(mu_branches["rr_mu_bestplane"][mask]))
    dE = np.array(ak.flatten(mu_branches["dE_mu_bestplane"][mask]))
    rr_rm = np.array(ak.flatten(mu_branches["rr_rm_mu_bestplane"][mask]))
    rm = np.array(ak.flatten(mu_branches["rm_mu_bestplane"][mask]))

    plt.plot(scatter_rr,scatter_dedx,lw=2,color="blue",label="Muons dE/dx vs RR Reference Curve", alpha=0.7)
    plt.plot(mip_rr,mip_dedx,lw=2,color='red',label="dE/dx = 2 MeV/cm", alpha=0.7)
    #plt.scatter(rr,dE,marker="o",s=40,color='cornflowerblue', edgecolors="black", linewidths=1.2, label=r"dE/dx vs RR, $\chi^2$ Stopping={:.2f} $\chi^2$ MIP={:.2f}".format(chi2_dedx,chi2_mip_dedx))
    #plt.scatter(rr_rm,rm,marker="^", s=40, color='orange', edgecolors="black", linewidths=1.2, label=r"Rolling Median vs RR, $\chi^2$ Stopping={:.2f} $\chi^2$ MIP={:.2f}".format(chi2,chi2_mip))
    
    plt.scatter(rr,dE,marker="o",s=40,color='cornflowerblue', edgecolors="black", linewidths=1.2, label="dE/dx vs RR")
    plt.scatter(rr_rm,rm,marker="^", s=40, color='orange', edgecolors="black", linewidths=1.2, label="Rolling Median vs RR")

    plt.xlabel('Residual Range [cm]', fontsize=18)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=18)
    colore_titolo = (34,73,120)
    colore_titolo_mpl = tuple(c / 255 for c in colore_titolo)
    plt.title(("RUN={:d} EVT={:d} - MC 2D OVERLAY - {:s}".format(run,evt,str_bestplane)), fontsize=18, color = colore_titolo_mpl)

    leg = plt.legend(loc='upper center', fontsize=10)
    leg.get_frame().set_facecolor('white')  # sfondo bianco della legenda
    leg.get_frame().set_alpha(1.0)          # opaco
    plt.setp(leg.get_title(), fontweight='bold')

    plt.setp(plt.gca().get_legend().get_title(), fontweight='bold')
    plt.xticks(np.arange(0,25.,1), fontsize=14, rotation=60)
    plt.yticks(fontsize=14)
    plt.ylim(1,15)
    plt.xlim(0,25)
    plt.savefig('grafici/new_grafici/example_muone_RUN_{:d}_EVT_{:d}.pdf'.format(run,evt), format='pdf', bbox_inches='tight')
    plt.close()

if False:
    increasing_dedx = np.loadtxt("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/dump_incresing_dedx.txt")
    incresing_rm = np.loadtxt("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/dump_incresing_rm.txt")
    NOT_increasing_dedx = np.loadtxt("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/dump_NONincresing_dedx.txt")
    NOT_increasing_rm = np.loadtxt("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/dump_NONincresing_rm.txt")

    rr_increasing_dedx = increasing_dedx[:,0]
    dedx_increasing_dedx = increasing_dedx[:,1]
    rr_increasing_rm = incresing_rm[:,0]
    dedx_increasing_rm = incresing_rm[:,1]
    rr_NONincreasing_dedx = NOT_increasing_dedx[:,0]
    dedx_NONincreasing_dedx = NOT_increasing_dedx[:,1]
    rr_NONincreasing_rm = NOT_increasing_rm[:,0]
    dedx_NONincreasing_rm = NOT_increasing_rm[:,1]

    colore_titolo = (34,73,120)
    colore_titolo_mpl = tuple(c / 255 for c in colore_titolo)

    h=plt.hist2d(rr_increasing_dedx, dedx_increasing_dedx, bins=(250,300), range=[(0,25),(0,30)], cmap='viridis', norm=LogNorm(), cmin=1)
    plt.xlabel('Residual Range [cm]', fontsize=17)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - COLLECTION PLANE - $\chi^2_{STOP.}\leq\chi^2_{MIP}$", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/incresing_dedx_pi.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    h=plt.hist2d(rr_NONincreasing_dedx, dedx_NONincreasing_dedx, bins=(250,300), range=[(0,25),(0,30)], cmap='viridis', norm=LogNorm(), cmin=1)
    plt.xlabel('Residual Range [cm]', fontsize=17)
    plt.ylabel('dE/dx [MeV/cm]', fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - COLLECTION PLANE - $\chi^2_{STOP.}>\chi^2_{MIP}$", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/NONincresing_dedx_pi.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    h=plt.hist2d(rr_increasing_rm, dedx_increasing_rm, bins=(250,300), range=[(0,25),(0,30)], cmap='viridis', norm=LogNorm(), cmin=1)
    plt.xlabel('Residual Range [cm]', fontsize=17)
    plt.ylabel('dE/dx Rolling Median [MeV/cm]', fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - COLLECTION PLANE - $\chi^2_{STOP.}\leq\chi^2_{MIP}$", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/incresing_rm_pi.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()

    h=plt.hist2d(rr_NONincreasing_rm, dedx_NONincreasing_rm, bins=(250,300), range=[(0,25),(0,30)], cmap='viridis', norm=LogNorm(), cmin=1)
    plt.xlabel('Residual Range [cm]', fontsize=17)
    plt.ylabel('dE/dx Rolling Median [MeV/cm]', fontsize=17)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r"$\pi^{\pm}$ MC 2D OVERLAYS - COLLECTION PLANE - $\chi^2_{STOP.}>\chi^2_{MIP}$", fontsize=18, color=colore_titolo_mpl)
    cbar=plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=14)  # cambia la dimensione dei numeri
    cbar.set_label('Frequency', fontsize=16)  # etichetta della barra
    plt.savefig('/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/grafici/new_grafici/NONincresing_rm_pi.png', format='png',dpi=600, bbox_inches='tight')
    plt.close()




