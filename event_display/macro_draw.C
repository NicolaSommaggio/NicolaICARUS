#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"


#include "helper_draw.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "THStack.h"
#include <TLine.h>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>


using namespace ana;

std::string pdg_to_string(int pdg)
{
    if(std::abs(pdg)==13){return "muon";}
    else if(pdg==2212){return "proton";}
    else if(std::abs(pdg)==211){return "pion";}
    else if(std::abs(pdg)==321){return "kaon";}
    else if(pdg==22){return "photon";}
    else if(abs(pdg)==11){return "electron";}
    else if(pdg==2112){return "neutron";}
    else if(pdg==111){return "neutral_pion";}
    else {return "other";}
}

void WriteTrackDrawInfo(SliceDrawInfo thislice, std::ofstream &file, int idx)
{
    file << endl;
    file << "#slice numero " << idx << " | nuE " << thislice.neutrino_interaction[0] << " parentPDG " << thislice.neutrino_interaction[1] << " Q2 " << thislice.neutrino_interaction[2] << endl; 
    file << "VERTEX " << thislice.vertex[0] << " " << thislice.vertex[1] << " " << thislice.vertex[2] << endl <<endl;
    for(int track=0; track < int(thislice.tracks_coordinate.size()); track++)
    {
        if(thislice.tracks_coordinate[track].size()==0)continue;
        bool is_muon=false;
        if(thislice.ipfp[track]==thislice.ipfp_mu)
        {
            is_muon=true;
            file << "#muone" << endl;
            for(int hit=0; hit<int(thislice.tracks_coordinate[track].size()); hit++)
            {
                file << thislice.tracks_coordinate[track][hit][0] << " " << thislice.tracks_coordinate[track][hit][1] << " " << thislice.tracks_coordinate[track][hit][2] << " " << thislice.tracks_dedx[track][hit] /**/<< " muon muon"/**/ << endl;
            }
        }
        bool is_proton=false;
        for(int ipro=0; ipro<int(thislice.v_ipfp_pro.size()); ipro++)
        {    
            if(thislice.ipfp[track]==thislice.v_ipfp_pro[ipro])
            {
                is_proton=true;
                file << "#protone" << endl;
                for(int hit=0; hit<int(thislice.tracks_coordinate[track].size()); hit++)
                {
                    file << thislice.tracks_coordinate[track][hit][0] << " " << thislice.tracks_coordinate[track][hit][1] << " " << thislice.tracks_coordinate[track][hit][2] << " " << thislice.tracks_dedx[track][hit] /**/<< Form(" proton%d",ipro) << " proton"/**/ << endl;
                }
            }
        }
        bool is_pion = false;
        for(int ipi=0; ipi<int(thislice.v_ipfp_pi.size()); ipi++)
        {    
            if(thislice.ipfp[track]==thislice.v_ipfp_pi[ipi])
            {
                is_pion=true;
                file << "#pione" << endl;
                for(int hit=0; hit<int(thislice.tracks_coordinate[track].size()); hit++)
                {
                    file << thislice.tracks_coordinate[track][hit][0] << " " << thislice.tracks_coordinate[track][hit][1] << " " << thislice.tracks_coordinate[track][hit][2] << " " << thislice.tracks_dedx[track][hit] /**/<< Form(" pion%d",ipi) << " pion"/**/ << endl;
                }
            }
        }
        if(!is_muon && !is_proton && !is_pion)
        {
            file << endl;
            for(int hit=0; hit<int(thislice.tracks_coordinate[track].size()); hit++)
            {
                file << thislice.tracks_coordinate[track][hit][0] << " " << thislice.tracks_coordinate[track][hit][1] << " " << thislice.tracks_coordinate[track][hit][2] << " " << thislice.tracks_dedx[track][hit] /**/<< Form(" track%d",track) << " other"/**/ << endl;
            }
        }


    }  

    /**/
    int true_part_count=0;
    for(auto const &itrue : thislice.true_particles)
    {
        file << endl;
        file << itrue.second[0] << " " << itrue.second[1] << " " << itrue.second[2] << " " << Form("true_part%d",true_part_count) << " " << pdg_to_string(itrue.first)  << " mc" << endl;
        file << itrue.second[3] << " " << itrue.second[4] << " " << itrue.second[5] << " " << Form("true_part%d",true_part_count) << " " << pdg_to_string(itrue.first)  << " mc" << endl;

        true_part_count++;
    }
    /**/

    file << "#new slice" << endl; 
}


ofstream dumpRollingMedian("dumpRollingMedian.txt");
void macro_draw(){
  
const std::string fdata = "/storage/gpfs_data/icarus/local/users/cfarnese/production_mc_2025A_ICARUS_Overlays_BNB_MC_RUN2_summer_2025_v10_06_00_04p04/MC_overlay_neutrino_stage1_flat_cafs_v10_06_00_04p04_concat.root";

//const std::string fdata= "/storage/gpfs_data/icarus/plain/user/cfarnese/test_genie_v0984_largestat_nuonly/concat_singleneu_2.flat.caf.root";

SpectrumLoader loader(fdata);      

const Binning kBinz     = Binning::Simple(300,0,30);

Spectrum s1("", kBinz, loader, disegna ,kCRTPMTNeutrino );

loader.Go();

double factor = s1.POT();  

TH1D* h1 = s1.ToTH1(factor);

ofstream tracce3D("tracce3Ddedx.txt");

for(int i=0; i<int(slices3D.size()); i++)
{
    WriteTrackDrawInfo(slices3D[i],tracce3D,i); 
}

TFile * f = new TFile("dumpRM.root","RECREATE");
TH2D * h_rolling_median = new TH2D("h_rolling_median","",200, 0, 20, 300, 0, 30);
TH2D * h_rolling_median_broken = new TH2D("h_rolling_median_broken","",200, 0, 20, 300, 0, 30);
TH1D *h_median = new TH1D("hmedian","",300,0,30);

cout << dumpRM.size() << " " << dumpMedian.size() << " " << dedx.size() << " " << conta_muoni << endl;
int conta=-1;
for(int i=0; i<int(dumpRM.size()); i++)
{
    for(int j=0; j<int(dumpRM[i].size()); j++) h_rolling_median->Fill(dumpRM[i][j].first,dumpRM[i][j].second);
    for(int j=0; j<int(dumpRM[i].size()); j++) {h_rolling_median->Fill(dumpRM[i][j].first,dumpRM[i][j].second);}
    //if(dumpMedian[i]<3.2 && dumpMedian[i]>1)
    //{
        conta++;
        for(int j=0; j<int(dumpRM[i].size()); j++) {h_rolling_median_broken->Fill(dumpRM[i][j].first,dumpRM[i][j].second);}
        dumpRollingMedian << "track " << conta << endl;
        for(int j=0; j<int(dedx[i].size()); j++) dumpRollingMedian << rr[i][j] << " ";
        dumpRollingMedian << endl;
        for(int j=0; j<int(dedx[i].size()); j++) dumpRollingMedian << dedx[i][j] << " ";
        dumpRollingMedian << endl;
        for(int j=0; j<int(dumpRM[i].size()); j++) dumpRollingMedian << dumpRM[i][j].first << " ";
        dumpRollingMedian << endl;
        for(int j=0; j<int(dumpRM[i].size()); j++) dumpRollingMedian << dumpRM[i][j].second << " ";
        dumpRollingMedian << endl;
        for(int j=0; j<int(x[i].size()); j++) dumpRollingMedian << x[i][j] << " ";
        dumpRollingMedian << endl;
    //}
    h_median->Fill(dumpMedian[i]);
}

h_rolling_median->Write(0,TObject::kOverwrite);
h_rolling_median_broken->Write(0,TObject::kOverwrite);
h_median->Write(0,TObject::kOverwrite);



}