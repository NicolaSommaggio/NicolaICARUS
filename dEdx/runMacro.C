void runMacro(int func)
{
    /*
    if (!gSystem->AccessPathName("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/FitResults.root"))
    {
        gSystem->Exec("rm \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/FitResults.root\"");
    } 
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("FitFileMaker(\"FitResults.root\")");
    
    if (!gSystem->AccessPathName("/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/GRAFICI.root"))
    {
        gSystem->Exec("rm \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/GRAFICI.root\"");
    } 
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("GeneralFileMaker(\"GRAFICI.root\")");

    gROOT->ProcessLine(".L macro.C");
    gROOT->ProcessLine("macro(\"mc\", \"muon\", \"Np\", \"25only\")");
    gROOT->ProcessLine("macro(\"mc\", \"proton\", \"Np\", \"25only\")");
    gROOT->ProcessLine(".q");
    */

    /*
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/OUTPUT_likelihood.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("likelihood_file()");
    gROOT->ProcessLine(".L macroTPC2.C");
    gROOT->ProcessLine("likelihood(\"muon\", \"muon\")");
    gROOT->ProcessLine(".q");
    */
   
    
    if(func==0){
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("confronto_file(\"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC.root\")");
    gROOT->ProcessLine(".L macroTPC2.C");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"mc\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"dati\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"mc\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"dati\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC.root\", \"\")");
    gROOT->ProcessLine(".q");
    }

        if(func==-1){
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_below30driftAngle.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("confronto_file(\"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_below30driftAngle.root\")");
    gROOT->ProcessLine(".L macroTPC2.C");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"mc\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_below30driftAngle.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"dati\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_below30driftAngle.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"mc\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_below30driftAngle.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"dati\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_below30driftAngle.root\", \"\")");
    gROOT->ProcessLine(".q");
    }

            if(func==-11){
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_above30driftAngle.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("confronto_file(\"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_above30driftAngle.root\")");
    gROOT->ProcessLine(".L macroTPC2.C");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"mc\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_above30driftAngle.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"dati\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_above30driftAngle.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"mc\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_above30driftAngle.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"dati\", \"\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_above30driftAngle.root\", \"\")");
    gROOT->ProcessLine(".q");
    }


    if(func==1){
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_corrected.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("confronto_file(\"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_corrected.root\")");
    gROOT->ProcessLine(".L macroTPC2.C");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"mc\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_corrected.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"dati\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_corrected.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"mc\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_corrected.root\", \"\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"dati\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/ConfrontoDatiMC_corrected.root\", \"\")");
    gROOT->ProcessLine(".q");
    }

 

    if(func==11){
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/25only_ConfrontoDatiMC_corrected.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("confronto_file(\"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/25only_ConfrontoDatiMC_corrected.root\")");
    gROOT->ProcessLine(".L macroTPC2.C");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"mc\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/25only_ConfrontoDatiMC_corrected.root\", \"25only\")");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"dati\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/25only_ConfrontoDatiMC_corrected.root\", \"25only\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"mc\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/25only_ConfrontoDatiMC_corrected.root\", \"25only\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"dati\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/25only_ConfrontoDatiMC_corrected.root\", \"25only\")");
    gROOT->ProcessLine(".q");
    }

    if(func==111){
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/specific_ConfrontoDatiMC_corrected.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("confronto_file(\"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/specific_ConfrontoDatiMC_corrected.root\")");
    gROOT->ProcessLine(".L macroTPC2.C");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"mc\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/specific_ConfrontoDatiMC_corrected.root\", \"specific\")");
    gROOT->ProcessLine("confrontoDatiMC(\"muon\", \"dati\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/specific_ConfrontoDatiMC_corrected.root\", \"specific\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"mc\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/specific_ConfrontoDatiMC_corrected.root\", \"specific\")");
    gROOT->ProcessLine("confrontoDatiMC(\"proton\", \"dati\", \"corrected\", \"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/MPV/specific_ConfrontoDatiMC_corrected.root\", \"specific\")");
    gROOT->ProcessLine(".q");
    }

    if(func==2){
    gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/mpvPlot.root");
    //gSystem->Exec("rm /storage/gpfs_data/icarus/local/users/sommaggio/simul_z/hReference.root");
    gROOT->ProcessLine(".L fileMaker.C");
    gROOT->ProcessLine("plotMPV_file(\"mpvPlot.root\")");
    //gROOT->ProcessLine("hReference_file()");
    gROOT->ProcessLine(".L plotMPV.C");
    gROOT->ProcessLine("plotMPV()");
    gROOT->ProcessLine("plotMPVpro()");
    gROOT->ProcessLine("dEdxCorrection()");
    //gROOT->ProcessLine("plotMPVcorrected()");
    //gSystem->Exec("python script.py");
    gROOT->ProcessLine(".q");
    }

    if(func==4)
    {
        gROOT->ProcessLine(".L plotMPV.C");
        gROOT->ProcessLine("chi2()");
        gSystem->Exec("python macro.py");
        gROOT->ProcessLine("plotHisto(\"dati_mu_ArgoNeut\", \"mc_mu_ArgoNeut\", \"mc_corr_mu_ArgoNeut\", \"\", \"canvas1\" )");
        gROOT->ProcessLine("plotHisto(\"dati_mu_ArgoNeut\", \"mc_broken_mu_ArgoNeut\", \"mc_broken_corr_mu_ArgoNeut\", \"\",  \"canvas2\")" );
        gROOT->ProcessLine("plotHisto(\"dati_pro_ArgoNeut\", \"mc_pro_ArgoNeut\", \"mc_corr_pro_ArgoNeut\", \"\",\"canvas3\" )");
        gROOT->ProcessLine("plotHisto(\"dati_mu_as_pro_ArgoNeut\", \"mc_mu_as_pro_ArgoNeut\", \"mc_corr_mu_as_pro_ArgoNeut\", \"\",\"canvas4\" )");
        gROOT->ProcessLine("plotHisto(\"dati_mu_as_pro_ArgoNeut\", \"mc_broken_mu_as_pro_ArgoNeut\", \"mc_broken_corr_mu_as_pro_ArgoNeut\",\"\", \"canvas5\")" );
        gROOT->ProcessLine("plotHisto(\"dati_mu_ArgoNeut\", \"dati_mu_mcICARUS\", \"dati_mu_mcICARUSspec\", \"\",\"canvas6\" )");
        gROOT->ProcessLine("plotHisto(\"mc_mu_ArgoNeut\", \"mc_mu_mcICARUS\", \"mc_mu_mcICARUSspec\", \"\",\"canvas7\")" );
        gROOT->ProcessLine("plotHisto(\"dati_pro_ArgoNeut\", \"dati_pro_mcICARUS\", \"dati_pro_mcICARUSspec\", \"\",\"canvas8\" )");
        gROOT->ProcessLine("plotHisto(\"mc_pro_ArgoNeut\", \"mc_pro_mcICARUS\", \"mc_pro_mcICARUSspec\", \"\",\"canvas9\" )");
        gROOT->ProcessLine("plotHisto(\"dati_mu_as_pro_ArgoNeut\", \"dati_mu_as_pro_mcICARUS\", \"dati_mu_as_pro_mcICARUSspec\", \"\",\"canvas10\")" );
        gROOT->ProcessLine("plotHisto(\"mc_mu_as_pro_ArgoNeut\", \"mc_mu_as_pro_mcICARUS\", \"mc_mu_as_pro_mcICARUSspec\", \"\",\"canvas11\")" );
        gROOT->ProcessLine("plotHisto(\"dati_pro_as_mu_ArgoNeut\", \"dati_pro_as_mu_mcICARUS\", \"dati_pro_as_mu_mcICARUSspec\", \"\",\"canvas12\")" );
        gROOT->ProcessLine("plotHisto(\"mc_pro_as_mu_ArgoNeut\", \"mc_pro_as_mu_mcICARUS\", \"mc_pro_as_mu_mcICARUSspec\", \"\",\"canvas13\")" );
        gROOT->ProcessLine("plotHisto(\"mc_broken_mu_ArgoNeut\", \"mc_broken_mu_mcICARUS\", \"mc_broken_mu_mcICARUSspec\", \"\",\"canvas14\")" );
        gROOT->ProcessLine("plotHisto(\"mc_broken_mu_as_pro_ArgoNeut\", \"mc_broken_mu_as_pro_mcICARUS\", \"mc_broken_mu_as_pro_mcICARUSspec\", \"\",\"canvas15\")" );

        gROOT->ProcessLine("plotHisto(\"mc_mu_ArgoNeut\", \"mc_corr_mu_ArgoNeut\", \"mc_mu_mcICARUSspec\", \"mc_corr_mu_mcICARUSspec\",\"canvas17\")" );
        gROOT->ProcessLine("plotHisto(\"mc_pro_ArgoNeut\", \"mc_corr_pro_ArgoNeut\", \"mc_pro_mcICARUSspec\", \"mc_corr_pro_mcICARUSspec\",\"canvas18\" )");
        gROOT->ProcessLine("plotHisto(\"mc_mu_as_pro_ArgoNeut\", \"mc_corr_mu_as_pro_ArgoNeut\", \"mc_mu_as_pro_mcICARUSspec\", \"mc_corr_mu_as_pro_mcICARUSspec\",\"canvas19\")" );
        gROOT->ProcessLine("plotHisto(\"mc_pro_as_mu_ArgoNeut\", \"mc_corr_pro_as_mu_ArgoNeut\", \"mc_pro_as_mu_mcICARUSspec\", \"mc_corr_pro_as_mu_mcICARUSspec\",\"canvas20\")" );
        gROOT->ProcessLine("plotHisto(\"mc_broken_mu_ArgoNeut\", \"mc_broken_corr_mu_ArgoNeut\", \"mc_broken_mu_mcICARUSspec\", \"mc_broken_corr_mu_mcICARUSspec\",\"canvas21\")" );
        gROOT->ProcessLine("plotHisto(\"mc_broken_mu_as_pro_ArgoNeut\", \"mc_broken_corr_mu_as_pro_ArgoNeut\", \"mc_broken_mu_as_pro_mcICARUSspec\", \"mc_broken_corr_mu_as_pro_mcICARUSspec\",\"canvas22\")" );

        gROOT->ProcessLine("plotHisto(\"dati_mu_ArgoNeut\", \"dati_mu_as_pro_ArgoNeut\", \"dati_mu_mcICARUSspec\", \"dati_mu_as_pro_mcICARUSspec\",\"canvas23\")" );
        gROOT->ProcessLine("plotHisto(\"dati_pro_ArgoNeut\", \"dati_pro_as_mu_ArgoNeut\", \"dati_pro_mcICARUSspec\", \"dati_pro_as_mu_mcICARUSspec\",\"canvas24\")" );
        gROOT->ProcessLine(".q");

        //gROOT->ProcessLine("plotHisto(\"dati_mu_ArgoNeut\", \"mc_mu_ArgoNeut\", \"\", \"canvas1\")");
        //gROOT->ProcessLine("plotHisto(\"dati_mu_ArgoNeut\", \"mc_corr_mu_ArgoNeut\", \"\", \"canvas2\")");
    }

    
}
