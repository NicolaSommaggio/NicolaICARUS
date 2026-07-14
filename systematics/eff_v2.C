#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TGraph2D.h>
#include <TCanvas.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"          // per TH2D
#include "TProfile.h"      // per TProfile
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"  // per TTreeReaderArray
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TVector3.h"      // per TVector3
#include <fstream>         // per ofstream, ifstream
#include <string>
#include <vector>
#include <cmath>           // per isnan

std::pair<int,int> wireToHit(std::vector<std::pair<int,int>> wires_hits, std::vector<TVector3> coordinates, int first_wire, int second_wire)
{
	std::vector<int> hits_first_wire;
	std::vector<int> hits_second_wire;
	for(int w = 0; w<wires_hits.size(); w++)
	{
		if(wires_hits[w].first == first_wire)
		{
			hits_first_wire.push_back(w);
		}
		if(wires_hits[w].first == second_wire)
		{
			hits_second_wire.push_back(w);
		}
	}

	if(hits_first_wire.size() == 1 && hits_second_wire.size() == 1)return {wires_hits[hits_first_wire[0]].second,wires_hits[hits_second_wire[0]].second};
	
	
	double best_distance = 1e6;
	std::pair<int,int> best_couple = {-1, -1};
	for(int i = 0; i < hits_first_wire.size(); i++)
	{
		for(int j = 0; j < hits_second_wire.size(); j++)
		{
			double distance = (coordinates[hits_first_wire[i]] - coordinates[hits_second_wire[j]]).Mag();
			if(distance < best_distance)
			{
				best_distance = distance;
				best_couple = {hits_first_wire[i],hits_second_wire[j]};
			}
		}
	}

	return{wires_hits[best_couple.first].second,wires_hits[best_couple.second].second};
	

}

//define main function
void efficiency(string file_list, int file_number, std::string sample_name) {


	bool PRINT_WEST_IND1_TPC0 = false;
	bool PRINT_WEST_IND1_TPC1 = false;
	bool PRINT_WEST_IND1_TPC2 = false;
	bool PRINT_WEST_IND1_TPC3 = false;
	bool PRINT_WEST_IND2_TPC01 = false;
	bool PRINT_WEST_IND2_TPC23 = false;
	bool PRINT_WEST_COLL_TPC01 = false;
	bool PRINT_WEST_COLL_TPC23 = false;

	bool PRINT_EAST_IND1_TPC0 = false;
	bool PRINT_EAST_IND1_TPC1 = false;
	bool PRINT_EAST_IND1_TPC2 = false;
	bool PRINT_EAST_IND1_TPC3 = false;
	bool PRINT_EAST_IND2_TPC01 = false;
	bool PRINT_EAST_IND2_TPC23 = false;
	bool PRINT_EAST_COLL_TPC01 = false;
	bool PRINT_EAST_COLL_TPC23 = false;


  	TH2D *hind1buchi=new TH2D("hind1buchi","",400,0,4.,50,-0.5,49.5);
  	TH2D *hind2buchi=new TH2D("hind2buchi","",400,0,4.,50,-0.5,49.5);
  	TH2D *hcollbuchi=new TH2D("hcollbuchi","",400,0,4.,50,-0.5,49.5);

	TH1F *h0_0_w=new TH1F("h0_0_w","",1056,-0.5,1055.5);
	TH1F *h0_1_w=new TH1F("h0_1_w","",1056,-0.5,1055.5);
	TH1F *h0_2_w=new TH1F("h0_2_w","",1056,-0.5,1055.5);
	TH1F *h0_3_w=new TH1F("h0_3_w","",1056,-0.5,1055.5);
	TH1F *h0_0_e=new TH1F("h0_0_e","",1056,-0.5,1055.5);
	TH1F *h0_1_e=new TH1F("h0_1_e","",1056,-0.5,1055.5);
	TH1F *h0_2_e=new TH1F("h0_2_e","",1056,-0.5,1055.5);
	TH1F *h0_3_e=new TH1F("h0_3_e","",1056,-0.5,1055.5);

	TH1F *h1_01_w=new TH1F("h1_01_w","",5600,-0.5,5599.5);
	TH1F *h1_23_w=new TH1F("h1_23_w","",5600,-0.5,5599.5);
	TH1F *h1_01_e=new TH1F("h1_01_e","",5600,-0.5,5599.5);
	TH1F *h1_23_e=new TH1F("h1_23_e","",5600,-0.5,5599.5);

	TH1F *h2_01_w=new TH1F("h2_01_w","",5600,-0.5,5599.5);
	TH1F *h2_23_w=new TH1F("h2_23_w","",5600,-0.5,5599.5);
	TH1F *h2_01_e=new TH1F("h2_01_e","",5600,-0.5,5599.5);
	TH1F *h2_23_e=new TH1F("h2_23_e","",5600,-0.5,5599.5);


	//generic counter
	float count_0 = 0;
	float count_1 = 0;
	float count_2 = 0;
	float count_3 = 0;
	float count_01 = 0;
	float count_23 = 0;

    float wcount_0 = 0;
    float wcount_1 = 0;
    float wcount_2 = 0;
    float wcount_3 = 0;
    float wcount_01 = 0;
    float wcount_23 = 0;

	//bins for histograms
	int bin;

	//generic holders
	float holder_0 = 0;
	float holder_1 = 0;
	float holder_2 = 0;
	float holder_3 = 0;
	float holder_01 = 0;
	float holder_23 = 0;


//		//
//		// DEFINE TTRE 
//		//


	int _run = -1;
	int _evt = -1;
	short _cryo = -1; //0 --> WEST		1 --> EAST
	short _tpc = -1; //0 e 2 --> IND2/COLL, 0,1,2,3 --> IND1 
	short _plane = -1;
	int _min_wire = -1;
	int _max_wire = -1;
	TVector3 _start;
	TVector3 _end;
	float _hit_wires = -1;
	float _tot_wires = -1;
	float _avg_pitch = -1;
	float _trk_length = -1;
	int _max_buco = -1;

	short _which_t0 = -1;
	float _t0PFP = -1;
	float _t0CRTTrack = -1;
	float _t0CRTHit = -1;

	int _nholes = -1;
	std::vector<std::vector<double>> _wire_holes; // --> it contains:
									 			  // _hole_dimension
												  // _hit_pos_before_hole_x,y,z
									 			  // _hit_pos_after_hole_x,y,z
											      // _last_wire_before_hole
									 			  // _first_wire_before_hole
									 			  // _hit_dir_before_buco
												  // _hit_dir_after_buco 
												  // _hit_time_before_buco
												  // _hit_time_after_buco
												  // _hit_mult_before_buco
												  // _hit_mult_after_buco

	int _G4ID = -999;
	float _truth_fraction = -1;

	TFile * tree_outfile = new TFile(Form("tree_outfile_%s.root",sample_name.c_str()),"RECREATE");
	tree_outfile -> cd();
	TTree * outree = new TTree("outree","outree");
	outree -> Branch("run",&_run);
	outree -> Branch("evt",&_evt);
	outree -> Branch("cryo",&_cryo);
	outree -> Branch("tpc",&_tpc);
	outree -> Branch("plane",&_plane);
	outree -> Branch("min_wire",&_min_wire);
	outree -> Branch("max_wire",&_max_wire);
	outree -> Branch("start",&_start);
	outree -> Branch("end",&_end);
	outree -> Branch("hit_wires",&_hit_wires);
	outree -> Branch("tot_wires",&_tot_wires);
	outree -> Branch("avg_pitch",&_avg_pitch);
	outree -> Branch("trk_length",&_trk_length);
	outree -> Branch("max_buco",&_max_buco);
	outree -> Branch("whicht0",&_which_t0);
	outree -> Branch("t0PFP", &_t0PFP);
	outree -> Branch("t0CRTTrack", &_t0CRTTrack);
	outree -> Branch("t0CRTHit", &_t0CRTHit);
	outree -> Branch("nholes",&_nholes);
	outree -> Branch("wire_holes",&_wire_holes);
	outree -> Branch("G4ID",&_G4ID);
	outree -> Branch("truth_fraction",&_truth_fraction);

//		//
//		//	DEFINE VARIABLES NEEDED FOR EFFICIENCY vs AVERAGE PITCH HISTOGRAMS
//		//

	//vectors to store efficiency and average pitch for each plane
	vector<float> hist_eff_0, hist_eff_1, hist_eff_2;
	vector<float> hist_average_pitch_0, hist_average_pitch_1, hist_average_pitch_2;

	//hold max and min wires for all tpcs
	int max_wire_0, max_wire_1, max_wire_2, max_wire_3, max_wire_01, max_wire_23;
	int min_wire_0, min_wire_1, min_wire_2, min_wire_3, min_wire_01, min_wire_23;

	//min amount of wires in each logic TPC needed to perform computations
	int wire_threshold = 25;

	//print traces with eff lower than upper_eff_print 
	float upper_eff_print = 0.5;
	//output file for traces with eff lower than upper_eff_print 
	//ofstream oFile("histograms_all_muons/low_eff_traces.txt");
    //ofstream outfiletracks("tracks_west_9435.txt");
    //ofstream outfiletracksE("tracks_east_9435.txt");


	//wire array format: wire_array_plane_logicTPC
	//fstream oFile("histograms_all_muons/low_eff_traces.txt");
	//
	
	//ind1
	short int wire_array_0_0[1056] = {0};
	short int wire_array_0_1[1056] = {0};
	short int wire_array_0_2[1056] = {0};
	short int wire_array_0_3[1056] = {0};
	
	//ind2
	short int wire_array_1_01[5600] = {0};
	short int wire_array_1_23[5600] = {0};
	
	//coll
	short int wire_array_2_01[5600] = {0};
	short int wire_array_2_23[5600] = {0};


    
    //ind1
    short int w_0_0_w[1056] = {0};
    short int w_0_1_w[1056] = {0};
    short int w_0_2_w[1056] = {0};
    short int w_0_3_w[1056] = {0};
    
    //ind2
    short int w_1_01_w[5600] = {0};
    short int w_1_23_w[5600] = {0};
    
    //coll
    short int w_2_01_w[5600] = {0};
    short int w_2_23_w[5600] = {0};

    //ind1
    short int w_0_0_e[1056] = {0};
    short int w_0_1_e[1056] = {0};
    short int w_0_2_e[1056] = {0};
    short int w_0_3_e[1056] = {0};
    
    //ind2
    short int w_1_01_e[5600] = {0};
    short int w_1_23_e[5600] = {0};
    
    //coll
    short int w_2_01_e[5600] = {0};
    short int w_2_23_e[5600] = {0};

   
	for(int ijk=0;ijk<1056;ijk++)w_0_0_w[ijk]=1;
	for(int ijk=0;ijk<1056;ijk++)w_0_1_w[ijk]=1;
	for(int ijk=0;ijk<1056;ijk++)w_0_2_w[ijk]=1;
	for(int ijk=0;ijk<1056;ijk++)w_0_3_w[ijk]=1;
	for(int ijk=0;ijk<1056;ijk++)w_0_0_e[ijk]=1;
	for(int ijk=0;ijk<1056;ijk++)w_0_1_e[ijk]=1;
	for(int ijk=0;ijk<1056;ijk++)w_0_2_e[ijk]=1;
	for(int ijk=0;ijk<1056;ijk++)w_0_3_e[ijk]=1;

	for(int ijk=0;ijk<5600;ijk++)w_1_01_w[ijk]=1;
	for(int ijk=0;ijk<5600;ijk++)w_1_23_w[ijk]=1;
	for(int ijk=0;ijk<5600;ijk++)w_2_01_w[ijk]=1;
	for(int ijk=0;ijk<5600;ijk++)w_2_23_w[ijk]=1;
	for(int ijk=0;ijk<5600;ijk++)w_1_01_e[ijk]=1;
	for(int ijk=0;ijk<5600;ijk++)w_1_23_e[ijk]=1;
	for(int ijk=0;ijk<5600;ijk++)w_2_01_e[ijk]=1;
	for(int ijk=0;ijk<5600;ijk++)w_2_23_e[ijk]=1;

 
//lista dead wires
    w_0_1_w[4]=0;
    w_0_1_w[5]=0;
    w_0_1_w[6]=0;
    w_0_1_w[7]=0;
    w_0_1_w[8]=0;
    w_0_1_w[9]=0;
    w_0_1_w[10]=0;
    w_0_2_w[29]=0;
    w_0_0_e[167]=0;
    w_0_2_e[204]=0;
    w_0_3_w[227]=0;
    w_0_0_w[319]=0;
    w_0_3_e[331]=0;
    w_0_2_w[503]=0;
    w_0_2_e[579]=0;
    w_0_3_w[607]=0;
    w_0_1_e[643]=0;
    w_0_1_e[644]=0;
    w_0_1_e[672]=0;
    w_0_1_e[674]=0;
    w_0_1_e[676]=0;
    w_0_1_e[678]=0;
    w_0_1_e[708]=0;
    w_0_1_e[736]=0;
    w_0_1_e[737]=0;
    w_0_1_e[738]=0;
    w_0_1_e[739]=0;
    w_0_3_e[761]=0;
    w_0_1_e[770]=0;
    w_0_1_e[861]=0;
    w_0_3_w[894]=0;
    w_0_0_w[992]=0;
    w_0_3_e[1041]=0;
    w_0_3_e[1042]=0;
    w_0_3_e[1043]=0;
    w_0_3_e[1044]=0;
    w_0_3_e[1045]=0;
    w_0_3_e[1046]=0;
    w_0_3_e[1047]=0;
    w_0_3_e[1048]=0;
    w_0_3_e[1049]=0;
    w_0_3_e[1050]=0;
    w_0_3_e[1051]=0;
    w_0_3_e[1052]=0;
    w_0_3_e[1053]=0;
    w_0_3_e[1054]=0;
    w_0_3_e[1055]=0;
    w_2_23_w[160]=0;
    w_2_23_w[161]=0;
    w_2_23_w[162]=0;
    w_2_23_w[163]=0;
    w_2_23_w[164]=0;
    w_2_23_w[165]=0;
    w_2_23_w[166]=0;
    w_2_23_w[167]=0;
    w_2_23_w[168]=0;
    w_2_23_w[169]=0;
    w_2_23_w[170]=0;
    w_2_23_w[171]=0;
    w_2_23_w[172]=0;
    w_2_23_w[173]=0;
    w_2_23_w[174]=0;
    w_2_23_w[175]=0;
    w_2_23_w[176]=0;
    w_2_23_w[177]=0;
    w_2_23_w[178]=0;
    w_2_23_w[179]=0;
    w_2_23_w[180]=0;
    w_2_23_w[181]=0;
    w_2_23_w[182]=0;
    w_2_23_w[183]=0;
    w_2_23_w[184]=0;
    w_2_23_w[185]=0;
    w_2_23_w[186]=0;
    w_2_23_w[187]=0;
    w_2_23_w[188]=0;
    w_2_23_w[189]=0;
    w_2_23_w[190]=0;
    w_2_23_w[191]=0;
    w_1_23_w[211]=0;
    w_2_01_e[395]=0;
    w_2_01_e[396]=0;
    w_2_01_e[397]=0;
    w_2_01_e[398]=0;
    w_2_01_e[399]=0;
    w_1_01_e[406]=0;
    w_1_01_w[448]=0;
    w_2_23_w[448]=0;
    w_1_01_w[449]=0;
    w_2_23_w[449]=0;
    w_1_01_w[450]=0;
    w_2_23_w[450]=0;
    w_1_01_w[451]=0;
    w_2_23_w[451]=0;
    w_1_01_w[452]=0;
    w_2_23_w[452]=0;
    w_1_01_w[453]=0;
    w_2_23_w[453]=0;
    w_1_01_w[454]=0;
    w_2_23_w[454]=0;
    w_1_01_w[455]=0;
    w_2_23_w[455]=0;
    w_1_01_w[456]=0;
    w_2_23_w[456]=0;
    w_1_01_w[457]=0;
    w_2_23_w[457]=0;
    w_1_01_w[458]=0;
    w_2_23_w[458]=0;
    w_1_01_w[459]=0;
    w_2_23_w[459]=0;
    w_1_01_w[460]=0;
    w_2_23_w[460]=0;
    w_1_01_w[461]=0;
    w_2_23_w[461]=0;
    w_1_01_w[462]=0;
    w_2_23_w[462]=0;
    w_1_01_w[463]=0;
    w_2_23_w[463]=0;
    w_1_01_w[464]=0;
    w_2_23_w[464]=0;
    w_1_01_w[465]=0;
    w_2_23_w[465]=0;
    w_1_01_w[466]=0;
    w_2_23_w[466]=0;
    w_1_01_w[467]=0;
    w_2_23_w[467]=0;
    w_1_01_w[468]=0;
    w_2_23_w[468]=0;
    w_1_01_w[469]=0;
    w_2_23_w[469]=0;
    w_1_01_w[470]=0;
    w_2_23_w[470]=0;
    w_1_01_w[471]=0;
    w_2_23_w[471]=0;
    w_1_01_w[472]=0;
    w_2_23_w[472]=0;
    w_1_01_w[473]=0;
    w_2_23_w[473]=0;
    w_1_01_w[474]=0;
    w_2_23_w[474]=0;
    w_1_01_w[475]=0;
    w_2_23_w[475]=0;
    w_1_01_w[476]=0;
    w_2_23_w[476]=0;
    w_1_01_w[477]=0;
    w_2_23_w[477]=0;
    w_1_01_w[478]=0;
    w_2_23_w[478]=0;
    w_1_01_w[479]=0;
    w_2_23_w[479]=0;
    w_1_23_w[570]=0;
    w_1_01_w[657]=0;
    w_1_23_e[745]=0;
    w_2_01_e[817]=0;
    w_1_23_e[1009]=0;
    w_2_23_w[1041]=0;
    w_2_01_e[1052]=0;
    w_2_23_e[1241]=0;
    w_2_01_e[1259]=0;
    w_2_01_e[1263]=0;
    w_1_01_w[1268]=0;
    w_1_01_w[1269]=0;
    w_1_01_e[1325]=0;
    w_1_01_e[1343]=0;
    w_2_23_e[1473]=0;
    w_2_01_e[1536]=0;
    w_1_01_w[1678]=0;
    w_1_01_e[1862]=0;
    w_2_23_e[2000]=0;
    w_2_23_e[2001]=0;
    w_2_01_e[2010]=0;
    w_2_23_w[2268]=0;
    w_2_01_w[2518]=0;
    w_1_01_e[3172]=0;
    w_1_23_e[3388]=0;
    w_2_23_e[3678]=0;
    w_2_01_e[3923]=0;
    w_1_01_w[4172]=0;
    w_2_23_e[5235]=0;
    w_2_23_e[5247]=0;
    w_2_23_e[5351]=0;
    w_1_23_w[5376]=0;
    w_0_0_w[1]=0;
    w_0_1_w[1]=0;
    w_0_2_w[1]=0;
    w_0_2_e[1]=0;
    w_0_3_e[1]=0;
    w_0_1_w[2]=0;
    w_0_1_e[2]=0;
    w_0_1_e[4]=0;
    w_0_0_e[5]=0;
    w_0_1_e[6]=0;
    w_0_0_e[7]=0;
    w_0_1_e[7]=0;
    w_0_0_e[8]=0;
    w_0_1_e[8]=0;
    w_0_3_e[18]=0;
    w_0_1_w[221]=0;
    w_0_0_w[354]=0;
    w_0_2_e[402]=0;
    w_0_2_e[530]=0;
    w_0_2_e[531]=0;
    w_0_2_w[677]=0;
    w_0_1_w[701]=0;
    w_0_2_e[812]=0;
    w_0_2_e[954]=0;
    w_0_2_e[956]=0;
    w_0_1_w[989]=0;
    w_0_3_e[1019]=0;
    w_0_0_e[1046]=0;
    w_0_1_e[1046]=0;
    w_0_0_e[1047]=0;
    w_0_1_e[1047]=0;
    w_0_0_e[1048]=0;
    w_0_1_e[1048]=0;
    w_0_1_e[1049]=0;
    w_0_0_w[1053]=0;
    w_0_0_e[1053]=0;
    w_0_1_e[1053]=0;
    w_0_2_e[1053]=0;
    w_0_1_w[1054]=0;
    w_0_3_w[1054]=0;
    w_0_2_e[1054]=0;
    w_1_23_e[31]=0;
    w_1_01_e[61]=0;
    w_2_23_e[64]=0;
    w_1_01_e[118]=0;
    w_1_01_e[121]=0;
    w_1_01_e[123]=0;
    w_1_01_e[156]=0;
    w_1_01_e[194]=0;
    w_1_01_e[197]=0;
    w_1_01_e[198]=0;
    w_1_01_e[201]=0;
    w_1_01_e[246]=0;
    w_1_01_e[251]=0;
    w_1_01_e[255]=0;
    w_1_01_e[261]=0;
    w_1_01_e[293]=0;
    w_2_23_e[355]=0;
    w_1_01_e[361]=0;
    w_1_01_e[451]=0;
    w_2_23_w[643]=0;
    w_2_01_w[674]=0;
    w_2_23_e[867]=0;
    w_2_23_e[892]=0;
    w_2_01_w[962]=0;
    w_2_01_w[1156]=0;
    w_2_01_w[1757]=0;
    w_2_23_e[2255]=0;
    w_2_23_w[2865]=0;
    w_1_01_w[2946]=0;
    w_1_01_w[3043]=0;
    w_2_23_e[3077]=0;
    w_1_01_w[3234]=0;
    w_1_01_w[3269]=0;
    w_1_23_e[3362]=0;
    w_1_01_w[3364]=0;
    w_2_01_w[3746]=0;
    w_1_23_w[3865]=0;
    w_2_01_w[4273]=0;
    w_1_01_w[4578]=0;
    w_1_01_w[4581]=0;
    w_1_23_e[5036]=0;
    w_2_01_w[5091]=0;
    w_2_01_w[5119]=0;
    w_1_23_e[5120]=0;
    w_1_23_e[5122]=0;
    w_1_01_e[5125]=0;
    w_1_23_e[5126]=0;
    w_1_23_e[5128]=0;
    w_1_23_e[5129]=0;
    w_1_23_e[5131]=0;
    w_1_01_w[5132]=0;
    w_1_23_e[5134]=0;
    w_1_23_e[5136]=0;
    w_1_23_e[5138]=0;
    w_1_23_e[5142]=0;
    w_1_23_e[5144]=0;
    w_1_23_e[5145]=0;
    w_1_23_e[5148]=0;
    w_1_23_e[5149]=0;
    w_1_23_e[5150]=0;
    w_2_01_e[5182]=0;
    w_1_23_e[5187]=0;
    w_1_23_w[5199]=0;
    w_1_23_e[5204]=0;
    w_1_23_e[5211]=0;
    w_1_23_e[5220]=0;
    w_1_23_w[5247]=0;
    w_2_01_w[5249]=0;
    w_1_23_e[5279]=0;
    w_1_23_e[5305]=0;
    w_1_23_e[5309]=0;
    w_1_23_e[5312]=0;
    w_2_01_e[5342]=0;
    w_1_23_e[5345]=0;
    w_1_23_e[5346]=0;
    w_2_01_e[5348]=0;
    w_1_23_e[5430]=0;
    w_1_23_e[5433]=0;
    w_1_23_e[5460]=0;
//
    
	int buco=0;
	int max_buco=0;


				
//      	//
//      	//	ADJUST NUMBER OF INPUT FILES
//      	// 

	string file_name;
	ifstream ifile(file_list.c_str());
	//read how many TFile are present in the list
	int line_count = 0;
    string line;
    while (getline(ifile, line)) 
	{
        line_count++;
    }
	//if input number is greater than total number of files, readjust 
	if(file_number > line_count) 
	{
		//warning message
		cout << endl << "### In " << file_list << " there are only " << line_count << " files ###" << endl << endl;	
		file_number = line_count;
	}
	//start again from the beginning of the file
	ifile.clear();
	ifile.seekg(0, ios::beg);

		



	
//      	//
//      	//	N-TUPLE DATA INITIALIZATION
//      	// 

	TFile *myFile;
    for(int ka = 0; ka < file_number; ka++) { 
    	
    //print name of the file to open
    ifile >> file_name;
	cout << " _______________________________________________________";
	cout << "____________" << endl;
    cout << endl << endl << file_name << " ## " << ka+1 << " ## " << endl;
    
	const char *c = file_name.c_str();
        	
    //open the file
    myFile = TFile::Open(c);
	if (!myFile || myFile->IsZombie()) 
	{
		return;
    }
    
	//define sudirectories
        	
	TTreeReader myReaderTPCW("caloskimW/TrackCaloSkim", myFile);
    TTreeReader myReaderTPCE("caloskimE/TrackCaloSkim", myFile);

		



//      	//
//      	//	TPC WEST INITIALIZATION
//      	//        	
        	
    //define TPC WEST variables for all planes
    TTreeReaderValue<int> runTPCW(myReaderTPCW, "meta.run"),
    					  evtTPCW(myReaderTPCW, "meta.evt"),
        				  cryostatTPCW(myReaderTPCW, "cryostat");
        	

	TTreeReaderValue<int> whicht0TPCW (myReaderTPCW, "whicht0");
    TTreeReaderValue<float>	lengthTPCW(myReaderTPCW, "length"), 
							t0TPCW(myReaderTPCW, "t0PFP"),
        					startxTPCW(myReaderTPCW, "start.x"),
        					startyTPCW(myReaderTPCW, "start.y"),
        					startzTPCW(myReaderTPCW, "start.z"),
							endxTPCW(myReaderTPCW, "end.x"),
							endyTPCW(myReaderTPCW, "end.y"),
							endzTPCW(myReaderTPCW, "end.z");
        				
    TTreeReaderValue<bool> clear_cosmic_muonTPCW_0(myReaderTPCW, "clear_cosmic_muon");

	TTreeReaderValue<float> t0_PFP_TPCW(myReaderTPCW, "t0PFP");
	TTreeReaderValue<float> t0_CRT_Track_TPCW(myReaderTPCW, "t0CRTTrack");
	TTreeReaderValue<float> t0_CRT_Hit_TPCW(myReaderTPCW, "t0CRTHit");

	TTreeReaderValue<int>   truthG4W(myReaderTPCW, "truth.p.G4ID");
        	
    //define TPC WEST variables for plane 0 (ind1)
    TTreeReaderArray<unsigned short> tpcTPCW_0(myReaderTPCW, "hits0.h.tpc");
        				
    TTreeReaderArray<bool> hasSPTPCW_0(myReaderTPCW, "hits0.h.hasSP");

	TTreeReaderArray<float> hitx_TPCW_0(myReaderTPCW, "hits0.h.sp.x");
	TTreeReaderArray<float> hity_TPCW_0(myReaderTPCW, "hits0.h.sp.y");
	TTreeReaderArray<float> hitz_TPCW_0(myReaderTPCW, "hits0.h.sp.z");

	TTreeReaderArray<float> hit_dirx_TPCW_0(myReaderTPCW, "hits0.dir.x");
	TTreeReaderArray<float> hit_diry_TPCW_0(myReaderTPCW, "hits0.dir.y");
	TTreeReaderArray<float> hit_dirz_TPCW_0(myReaderTPCW, "hits0.dir.z");

	TTreeReaderArray<unsigned short> hit_mult_TPCW_0(myReaderTPCW, "hits0.h.mult");

	TTreeReaderArray<float> hit_time_TPCW_0(myReaderTPCW, "hits0.h.time");
        				
    TTreeReaderArray<unsigned short> wireTPCW_0(myReaderTPCW, "hits0.h.wire");     				
        				
    TTreeReaderArray<float>	dqdxTPCW_0(myReaderTPCW, "hits0.dqdx"), 
							pitchTPCW_0(myReaderTPCW, "hits0.pitch");

	TTreeReaderArray<float>    h0_truth_ne_W  (myReaderTPCW, "hits0.h.truth.nelec");
    TTreeReaderArray<float>    h0_truth_e_W  (myReaderTPCW, "hits0.h.truth.e");
       					
    //define TPC WEST variables for plane 1 (ind2)
    TTreeReaderArray<unsigned short> tpcTPCW_1(myReaderTPCW, "hits1.h.tpc");
        				
    TTreeReaderArray<bool> hasSPTPCW_1(myReaderTPCW, "hits1.h.hasSP");

	TTreeReaderArray<float> hitx_TPCW_1(myReaderTPCW, "hits1.h.sp.x");
	TTreeReaderArray<float> hity_TPCW_1(myReaderTPCW, "hits1.h.sp.y");
	TTreeReaderArray<float> hitz_TPCW_1(myReaderTPCW, "hits1.h.sp.z");

	TTreeReaderArray<float> hit_dirx_TPCW_1(myReaderTPCW, "hits1.dir.x");
	TTreeReaderArray<float> hit_diry_TPCW_1(myReaderTPCW, "hits1.dir.y");
	TTreeReaderArray<float> hit_dirz_TPCW_1(myReaderTPCW, "hits1.dir.z");

	TTreeReaderArray<unsigned short> hit_mult_TPCW_1(myReaderTPCW, "hits1.h.mult");

	TTreeReaderArray<float> hit_time_TPCW_1(myReaderTPCW, "hits1.h.time");
        				
    TTreeReaderArray<unsigned short> wireTPCW_1(myReaderTPCW, "hits1.h.wire");     				
        				
    TTreeReaderArray<float>	dqdxTPCW_1(myReaderTPCW, "hits1.dqdx"),
							pitchTPCW_1(myReaderTPCW, "hits1.pitch");

	TTreeReaderArray<float>    h1_truth_ne_W  (myReaderTPCW, "hits1.h.truth.nelec");
    TTreeReaderArray<float>    h1_truth_e_W  (myReaderTPCW, "hits1.h.truth.e");
       					
    //define TPC WEST variables for plane 2 (coll)
    TTreeReaderArray<unsigned short> tpcTPCW_2(myReaderTPCW, "hits2.h.tpc");
        				
    TTreeReaderArray<bool> hasSPTPCW_2(myReaderTPCW, "hits2.h.hasSP");

	TTreeReaderArray<float> hitx_TPCW_2(myReaderTPCW, "hits2.h.sp.x");
	TTreeReaderArray<float> hity_TPCW_2(myReaderTPCW, "hits2.h.sp.y");
	TTreeReaderArray<float> hitz_TPCW_2(myReaderTPCW, "hits2.h.sp.z");

	TTreeReaderArray<float> hit_dirx_TPCW_2(myReaderTPCW, "hits2.dir.x");
	TTreeReaderArray<float> hit_diry_TPCW_2(myReaderTPCW, "hits2.dir.y");
	TTreeReaderArray<float> hit_dirz_TPCW_2(myReaderTPCW, "hits2.dir.z");

	TTreeReaderArray<unsigned short> hit_mult_TPCW_2(myReaderTPCW, "hits2.h.mult");

	TTreeReaderArray<float> hit_time_TPCW_2(myReaderTPCW, "hits2.h.time");
        				
    TTreeReaderArray<unsigned short> wireTPCW_2(myReaderTPCW, "hits2.h.wire");     				
        				
	TTreeReaderArray<float>	dqdxTPCW_2(myReaderTPCW, "hits2.dqdx"),
       					 	pitchTPCW_2(myReaderTPCW, "hits2.pitch");

	TTreeReaderArray<float>    h2_truth_ne_W  (myReaderTPCW, "hits2.h.truth.nelec");
    TTreeReaderArray<float>    h2_truth_e_W  (myReaderTPCW, "hits2.h.truth.e");

		
//      	//
//      	//	TPC EAST INITIALIZATION
//      	//        	
        	
    //define TPC EST variables for all planes
    TTreeReaderValue<int> runTPCE(myReaderTPCE, "meta.run"),
        				  evtTPCE(myReaderTPCE, "meta.evt"),
        				  cryostatTPCE(myReaderTPCE, "cryostat");
        	
    TTreeReaderValue<int> whicht0TPCE (myReaderTPCE, "whicht0");
    TTreeReaderValue<float>	lengthTPCE(myReaderTPCE, "length"),
	    					t0TPCE(myReaderTPCE, "t0PFP"),
        					startxTPCE(myReaderTPCE, "start.x"),
        					startyTPCE(myReaderTPCE, "start.y"),
        					startzTPCE(myReaderTPCE, "start.z"),
							endxTPCE(myReaderTPCE, "end.x"),
							endyTPCE(myReaderTPCE, "end.y"),
							endzTPCE(myReaderTPCE, "end.z");
        				
    TTreeReaderValue<bool>	clear_cosmic_muonTPCE_0(myReaderTPCE, "clear_cosmic_muon");

	TTreeReaderValue<float> t0_PFP_TPCE(myReaderTPCE, "t0PFP");
	TTreeReaderValue<float> t0_CRT_Track_TPCE(myReaderTPCE, "t0CRTTrack");
	TTreeReaderValue<float> t0_CRT_Hit_TPCE(myReaderTPCE, "t0CRTHit");

	TTreeReaderValue<int>   truthG4E(myReaderTPCE, "truth.p.G4ID");
        	
    //define TPC EAST variables for plane 0 (ind1)
    TTreeReaderArray<unsigned short> tpcTPCE_0(myReaderTPCE, "hits0.h.tpc");
        				
    TTreeReaderArray<bool> hasSPTPCE_0(myReaderTPCE, "hits0.h.hasSP");

	TTreeReaderArray<float> hitx_TPCE_0(myReaderTPCE, "hits0.h.sp.x");
	TTreeReaderArray<float> hity_TPCE_0(myReaderTPCE, "hits0.h.sp.y");
	TTreeReaderArray<float> hitz_TPCE_0(myReaderTPCE, "hits0.h.sp.z");

	TTreeReaderArray<float> hit_dirx_TPCE_0(myReaderTPCE, "hits0.dir.x");
	TTreeReaderArray<float> hit_diry_TPCE_0(myReaderTPCE, "hits0.dir.y");
	TTreeReaderArray<float> hit_dirz_TPCE_0(myReaderTPCE, "hits0.dir.z");

	TTreeReaderArray<unsigned short> hit_mult_TPCE_0(myReaderTPCE, "hits0.h.mult");

	TTreeReaderArray<float> hit_time_TPCE_0(myReaderTPCE, "hits0.h.time");
        				
    TTreeReaderArray<unsigned short> wireTPCE_0(myReaderTPCE, "hits0.h.wire");     				
        				
    TTreeReaderArray<float>	dqdxTPCE_0(myReaderTPCE, "hits0.dqdx"),
       						pitchTPCE_0(myReaderTPCE, "hits0.pitch");

	TTreeReaderArray<float>    h0_truth_ne_E  (myReaderTPCE, "hits0.h.truth.nelec");
    TTreeReaderArray<float>    h0_truth_e_E  (myReaderTPCE, "hits0.h.truth.e");
       					
    //define TPC EAST variables for plane 1 (ind2)
    TTreeReaderArray<unsigned short> tpcTPCE_1(myReaderTPCE, "hits1.h.tpc");
        				
    TTreeReaderArray<bool>	hasSPTPCE_1(myReaderTPCE, "hits1.h.hasSP");

	TTreeReaderArray<float> hitx_TPCE_1(myReaderTPCE, "hits1.h.sp.x");
	TTreeReaderArray<float> hity_TPCE_1(myReaderTPCE, "hits1.h.sp.y");
	TTreeReaderArray<float> hitz_TPCE_1(myReaderTPCE, "hits1.h.sp.z");

	TTreeReaderArray<float> hit_dirx_TPCE_1(myReaderTPCE, "hits1.dir.x");
	TTreeReaderArray<float> hit_diry_TPCE_1(myReaderTPCE, "hits1.dir.y");
	TTreeReaderArray<float> hit_dirz_TPCE_1(myReaderTPCE, "hits1.dir.z");

	TTreeReaderArray<unsigned short> hit_mult_TPCE_1(myReaderTPCE, "hits1.h.mult");

	TTreeReaderArray<float> hit_time_TPCE_1(myReaderTPCE, "hits1.h.time");
        				
    TTreeReaderArray<unsigned short> wireTPCE_1(myReaderTPCE, "hits1.h.wire");     				
        				
    TTreeReaderArray<float>	dqdxTPCE_1(myReaderTPCE, "hits1.dqdx"),
       						pitchTPCE_1(myReaderTPCE, "hits1.pitch");

	TTreeReaderArray<float>    h1_truth_ne_E  (myReaderTPCE, "hits1.h.truth.nelec");
    TTreeReaderArray<float>    h1_truth_e_E  (myReaderTPCE, "hits1.h.truth.e");
       					
    //define TPC EAST variables for plane 2 (coll)
    TTreeReaderArray<unsigned short> tpcTPCE_2(myReaderTPCE, "hits2.h.tpc");
        				
    TTreeReaderArray<bool> hasSPTPCE_2(myReaderTPCE, "hits2.h.hasSP");

	TTreeReaderArray<float> hitx_TPCE_2(myReaderTPCE, "hits2.h.sp.x");
	TTreeReaderArray<float> hity_TPCE_2(myReaderTPCE, "hits2.h.sp.y");
	TTreeReaderArray<float> hitz_TPCE_2(myReaderTPCE, "hits2.h.sp.z");

	TTreeReaderArray<float> hit_dirx_TPCE_2(myReaderTPCE, "hits2.dir.x");
	TTreeReaderArray<float> hit_diry_TPCE_2(myReaderTPCE, "hits2.dir.y");
	TTreeReaderArray<float> hit_dirz_TPCE_2(myReaderTPCE, "hits2.dir.z");

	TTreeReaderArray<unsigned short> hit_mult_TPCE_2(myReaderTPCE, "hits2.h.mult");

	TTreeReaderArray<float> hit_time_TPCE_2(myReaderTPCE, "hits2.h.time");
        				
    TTreeReaderArray<unsigned short> wireTPCE_2(myReaderTPCE, "hits2.h.wire");     				
        				
    TTreeReaderArray<float>	dqdxTPCE_2(myReaderTPCE, "hits2.dqdx"),
       						pitchTPCE_2(myReaderTPCE, "hits2.pitch");

	TTreeReaderArray<float>    h2_truth_ne_E  (myReaderTPCE, "hits2.h.truth.nelec");
    TTreeReaderArray<float>    h2_truth_e_E  (myReaderTPCE, "hits2.h.truth.e");

//
//	EFFICIENCY OVER THE RR LIMIT FOR ALL PLANES
//

//
// WEST
//

myReaderTPCW.Restart();
while(myReaderTPCW.Next())
{	
	//if(*truthG4W == -1)continue; // --> MC in Overlays --> HERE
	//if(*truthG4W != -1)continue; // --> DATA in Overlays

	int bestplane = -1;
	std::array<float,3> nhit = {0,0,0};
	std::array<float,3> n_valid_hit = {0,0,0};
	float truth_fraction = 0;

	for(int hit = 0; hit < tpcTPCW_0.GetSize(); hit++)
	{
		nhit[0] ++;
		if(h0_truth_ne_W[hit] > 0 && h0_truth_e_W[hit] > 0){n_valid_hit[0]++;}
	}
	for(int hit = 0; hit < tpcTPCW_1.GetSize(); hit++)
	{
		nhit[1] ++;
		if(h1_truth_ne_W[hit] > 0 && h1_truth_e_W[hit] > 0){n_valid_hit[1]++;}
	}
	for(int hit = 0; hit < tpcTPCW_2.GetSize(); hit++)
	{
		nhit[2] ++;
		if(h2_truth_ne_W[hit] > 0 && h2_truth_e_W[hit] > 0){n_valid_hit[2]++;}
	}

	bestplane = 0;
	for (int i = 1; i < 3; ++i) 
	{
    	if (nhit[i] > nhit[bestplane])
        bestplane = i;
	}

	if (nhit[bestplane] > 0){truth_fraction = n_valid_hit[bestplane] / nhit[bestplane];}

	//if(truth_fraction < 0.8)continue; --> HERE

	if(false)
	{
		cout << "WEST *** NEW TRACK, G4ID: " << *truthG4W << ", TRUTH FRACTION: " << truth_fraction << endl << endl;
		cout << "plane0 " << tpcTPCW_0.GetSize() << " " << nhit[0] << " | plane1 " << tpcTPCW_1.GetSize() << " " << nhit[1] << " | plane2 " << tpcTPCW_2.GetSize() << " " << nhit[2] << " " << endl << endl;
		for(int hit = 0; hit < h2_truth_ne_W.GetSize(); hit++)
		{
			cout << h2_truth_ne_W[hit] << " ";
		}
		cout << endl;
		for(int hit = 0; hit < h2_truth_e_W.GetSize(); hit++)
		{
			cout << h2_truth_e_W[hit] << " ";
		}
		cout << endl << endl;
	}

	//muon length must be greater than 1 meter
	if(*lengthTPCW > 100 && ( (*whicht0TPCW)==0 || (*whicht0TPCW)==2 )) 
	{
		//outfiletracks << (*runTPCW) << " " << (*evtTPCW) <<  " 1 " << (*startxTPCW) << " " << (*startyTPCW) << " " << (*startzTPCW) << " " << (*endxTPCW) << " " << (*endyTPCW) << " " << (*endzTPCW) << " " << (*lengthTPCW) << endl;
		
		//
		//	PLANE IND1 WEST
		//
				
		//loop over all hits of the trace
		for(int j = 0; j < wireTPCW_0.GetSize(); ++j) 
		{
			//conditions to consider a valid hit, can change
			if(dqdxTPCW_0[j] > 0 && pitchTPCW_0[j] < 4 && pitchTPCW_0[j] > 0 && !isnan(pitchTPCW_0[j])) 
			{
				if(tpcTPCW_0[j] == 0) 
				{
					holder_0 += pitchTPCW_0[j];
					++count_0;
				}
				else if(tpcTPCW_0[j] == 1) 
				{
					holder_1 += pitchTPCW_0[j];
					++count_1;
				}
				else if(tpcTPCW_0[j] == 2) 
				{
					holder_2 += pitchTPCW_0[j];
					++count_2;
				}
				else 
				{
					holder_3 += pitchTPCW_0[j];		
					++count_3;
				}
			}
				
			//write down valid wire hits, less restrictions
			if(hasSPTPCW_0[j])
			{
				if(tpcTPCW_0[j] == 0) wire_array_0_0[wireTPCW_0[j]] = 1;
				else if(tpcTPCW_0[j] == 1) wire_array_0_1[wireTPCW_0[j]] = 1;
				else if(tpcTPCW_0[j] == 2) wire_array_0_2[wireTPCW_0[j]] = 1;
				else wire_array_0_3[wireTPCW_0[j]] = 1;
            }
		}

		//check if there are at least the minimum amount of wires in the TPCs	
		/*
            for(int k = 0; k < 1056; ++k) 
			{

				if(wire_array_0_0[k] == 1)h0_0_w->Fill(k);
				if(wire_array_0_1[k] == 1)h0_1_w->Fill(k);
				if(wire_array_0_2[k] == 1)h0_2_w->Fill(k);
				if(wire_array_0_3[k] == 1)h0_3_w->Fill(k);

 			}
		*/				
		
		// TPC0 IND1
		if(count_0 > wire_threshold) 
		{
			//average pitch
			holder_0 = holder_0/count_0;
			//find min wire
			for(int k = 0; k < 1056; ++k) 
			{
				if(wire_array_0_0[k] == 1) 
				{
					min_wire_0 = k;
					break;		
				}
			}

			//find max wire
			//for(int k = 1056; k > 0; --k) 
			for(int k = 1055; k > 0; --k) 
			{
				if(wire_array_0_0[k] == 1) 
				{
					max_wire_0 = k;
					break;		
				}
			}
					
			//find total amount of wires in the trace
            count_0 = 0; buco=0; wcount_0=0; max_buco=0;


			//cout << "wire <-> hit" << endl;
			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_0.GetSize(); hit++)
			{
				if(tpcTPCW_0[hit] == 0 /*hasSPTPCW_0[hit]*/) // SP request already in wire_array... == 1
				{
					if(PRINT_WEST_IND1_TPC0)
					{
						cout << wireTPCW_0[hit] << " " << hit_time_TPCW_0[hit] << " " << hasSPTPCW_0[hit] << " " << (wireTPCW_0[hit] >= min_wire_0 && wireTPCW_0[hit] <= max_wire_0) << " " << w_0_0_w[wireTPCW_0[hit]] << " ";
						if(!isnan(hitx_TPCW_0[hit])){ cout << hitx_TPCW_0[hit] << " " << hity_TPCW_0[hit] << " " << hitz_TPCW_0[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_0[hit] >= min_wire_0 && wireTPCW_0[hit] <= max_wire_0)
					{
						wire_hits.push_back({wireTPCW_0[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_0[hit],hity_TPCW_0[hit],hitz_TPCW_0[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}

			double start_buco = -1;
			for(int k = min_wire_0; k <= max_wire_0; ++k) 
			{ 
                wcount_0+=w_0_0_w[k];
                if(w_0_0_w[k]>0 && wire_array_0_0[k] == 1 && buco>0) 
				{
					++count_0; 
					
					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_0_0[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_0[hh.first]);
					this_hole.push_back(hity_TPCW_0[hh.first]);
					this_hole.push_back(hitz_TPCW_0[hh.first]);
					this_hole.push_back(hitx_TPCW_0[hh.second]);
					this_hole.push_back(hity_TPCW_0[hh.second]);
					this_hole.push_back(hitz_TPCW_0[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_0[hh.first]);
					this_hole.push_back(hit_diry_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_0[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_0[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.second]);
					this_hole.push_back(hit_time_TPCW_0[hh.first]);
					this_hole.push_back(hit_time_TPCW_0[hh.second]);
					this_hole.push_back(hit_mult_TPCW_0[hh.first]);
					this_hole.push_back(hit_mult_TPCW_0[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hind1buchi->Fill(holder_0,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;
				}
				else if(w_0_0_w[k]>0 && wire_array_0_0[k] == 1 && buco==0) ++count_0;
				else if(w_0_0_w[k]>0 && wire_array_0_0[k] == 0)
				{ 
					if(buco == 0)start_buco = k;
					++buco;
				}
			}

			wire_hits.clear();
					
			//store efficiency and average pitch in the ind1 vectors
			if(wcount_0 > (max_wire_0-min_wire_0+1) || count_0 > wcount_0)cout << "IND 0 0 W " << count_0 << " " << wcount_0 << " " << min_wire_0 << " " << max_wire_0 << endl; 
			if(max_buco<10)hist_eff_0.push_back(count_0/(wcount_0));
			if(max_buco<10)hist_average_pitch_0.push_back(holder_0);
			if(holder_0 < 0.2) 
			{
				for(int i = 0; i < pitchTPCW_0.GetSize(); ++i) 
				{
					cout << pitchTPCW_0[i] << endl;
				}
				return;
			}	
			
			cout << "WEST CRYO, PLANE IND1, LOGIC TPC 0 ";
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE IND1, LOGIC TPC 0
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 0;
			_plane = 0;
			_min_wire = min_wire_0;
			_max_wire = max_wire_0;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_0;
			_tot_wires = wcount_0;
			_avg_pitch = holder_0;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
			_which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			//cout << _start.X() << " " << _start.Y() << " " << _start.Z() << " " << _end.X() << " " << _end.Y() << " " << _end.Z() << " | " << " " << _nholes << " holes: ";
			if(PRINT_WEST_IND1_TPC0)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();
		}
				
		// TPC1 IND1
		if(count_1 > wire_threshold) 
		{
			//average pitch
			holder_1 = holder_1/count_1;
			//find min wire
			for(int k = 0; k < 1056; ++k) 
			{
				if(wire_array_0_1[k] == 1) 
				{
					min_wire_1 = k;
					break;		
				}
			}

			//find max wire
			//for(int k = 1056; k > 0; --k) 
			for(int k = 1055; k > 0; --k) 
			{
				if(wire_array_0_1[k] == 1) 
				{
					max_wire_1 = k;
					break;		
				}
			}

			//find total amount of wires in the trace
			count_1 = 0; wcount_1 = 0; buco=0;max_buco=0;

			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_0.GetSize(); hit++)
			{
				if(tpcTPCW_0[hit] == 1 /*&& hasSPTPCW_0[hit]*/)
				{
					if(PRINT_WEST_IND1_TPC1)
					{
						cout << wireTPCW_0[hit] << " " << hit_time_TPCW_0[hit] << " " << hasSPTPCW_0[hit] << " " << (wireTPCW_0[hit] >= min_wire_1 && wireTPCW_0[hit] <= max_wire_1) << " " << w_0_1_w[wireTPCW_0[hit]] << " ";
						if(!isnan(hitx_TPCW_0[hit])){ cout << hitx_TPCW_0[hit] << " " << hity_TPCW_0[hit] << " " << hitz_TPCW_0[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_0[hit] >= min_wire_1 && wireTPCW_0[hit] <= max_wire_1)
					{
						wire_hits.push_back({wireTPCW_0[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_0[hit],hity_TPCW_0[hit],hitz_TPCW_0[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}

			double start_buco = -1;
			for(int k = min_wire_1; k <= max_wire_1; ++k) 
			{

                wcount_1+=w_0_1_w[k];
                if(w_0_1_w[k]>0 && wire_array_0_1[k] == 1 && buco>0) 
				{
					++count_1;

					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_0_1[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_0[hh.first]);
					this_hole.push_back(hity_TPCW_0[hh.first]);
					this_hole.push_back(hitz_TPCW_0[hh.first]);
					this_hole.push_back(hitx_TPCW_0[hh.second]);
					this_hole.push_back(hity_TPCW_0[hh.second]);
					this_hole.push_back(hitz_TPCW_0[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_0[hh.first]);
					this_hole.push_back(hit_diry_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_0[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_0[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.second]);
					this_hole.push_back(hit_time_TPCW_0[hh.first]);
					this_hole.push_back(hit_time_TPCW_0[hh.second]);
					this_hole.push_back(hit_mult_TPCW_0[hh.first]);
					this_hole.push_back(hit_mult_TPCW_0[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hind1buchi->Fill(holder_1,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;
				}
			  	else if(w_0_1_w[k]>0 && wire_array_0_1[k] == 1 && buco==0) ++count_1;
			  	else if(w_0_1_w[k]>0 && wire_array_0_1[k] == 0)
				{ 
					if(buco == 0)start_buco = k;
					++buco;
				}
			}

			wire_hits.clear();

			//count_1 = 0;
			//for(int k = min_wire_1; k <= max_wire_1; ++k) {
			//	if(wire_array_0_1[k] == 1) ++count_1;
			//}
			//store efficiency and average pitch in the ind1 vectors

            if(wcount_1 > (max_wire_1-min_wire_1+1) || count_1 > wcount_1)cout << "IND 0 1 W " << count_1 << " " << wcount_1 << " " << min_wire_1 << " " << max_wire_1 << endl;
			if(max_buco<10)hist_eff_0.push_back(count_1/(wcount_1));
			if(max_buco<10)hist_average_pitch_0.push_back(holder_1);

			cout << "WEST CRYO, PLANE IND1, LOGIC TPC 1 ";
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE IND1, LOGIC TPC 1
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 1;
			_plane = 0;
			_min_wire = min_wire_1;
			_max_wire = max_wire_1;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_1;
			_tot_wires = wcount_1;
			_avg_pitch = holder_1;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
			_which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			if(PRINT_WEST_IND1_TPC1)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();
		}
				
		// TPC2 IND1	
		if(count_2 > wire_threshold) 
		{
			//average pitch
			holder_2 = holder_2/count_2;
			//find min wire
			for(int k = 0; k < 1056; ++k) 
			{
				if(wire_array_0_2[k] == 1) 
				{
					min_wire_2 = k;
					break;		
				}
			}

			//find max wire
			//for(int k = 1056; k > 0; --k) 
			for(int k = 1055; k > 0; --k)
			{
				if(wire_array_0_2[k] == 1) 
				{
					max_wire_2 = k;
					break;		
				}
			}
					
			//find total amount of wires in the trace
			count_2 = 0; wcount_2 = 0; buco=0;max_buco=0;

			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_0.GetSize(); hit++)
			{
				if(tpcTPCW_0[hit] == 2 /*&& hasSPTPCW_0[hit]*/)
				{
					if(PRINT_WEST_IND1_TPC2)
					{
						cout << wireTPCW_0[hit] << " " << hit_time_TPCW_0[hit] << " " << hasSPTPCW_0[hit] << " " << (wireTPCW_0[hit] >= min_wire_2 && wireTPCW_0[hit] <= max_wire_2) << " " << w_0_2_w[wireTPCW_0[hit]] << " ";
						if(!isnan(hitx_TPCW_0[hit])){ cout << hitx_TPCW_0[hit] << " " << hity_TPCW_0[hit] << " " << hitz_TPCW_0[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_0[hit] >= min_wire_2 && wireTPCW_0[hit] <= max_wire_2)
					{
						wire_hits.push_back({wireTPCW_0[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_0[hit],hity_TPCW_0[hit],hitz_TPCW_0[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}
					
			double start_buco = -1;
			for(int k = min_wire_2; k <= max_wire_2; ++k) 
			{
                wcount_2+=w_0_2_w[k];
				if(w_0_2_w[k]>0 && wire_array_0_2[k] == 1 && buco>0) 
				{
					++count_2;

					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_0_2[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_0[hh.first]);
					this_hole.push_back(hity_TPCW_0[hh.first]);
					this_hole.push_back(hitz_TPCW_0[hh.first]);
					this_hole.push_back(hitx_TPCW_0[hh.second]);
					this_hole.push_back(hity_TPCW_0[hh.second]);
					this_hole.push_back(hitz_TPCW_0[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_0[hh.first]);
					this_hole.push_back(hit_diry_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_0[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_0[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.second]);
					this_hole.push_back(hit_time_TPCW_0[hh.first]);
					this_hole.push_back(hit_time_TPCW_0[hh.second]);
					this_hole.push_back(hit_mult_TPCW_0[hh.first]);
					this_hole.push_back(hit_mult_TPCW_0[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hind1buchi->Fill(holder_2,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;
				}
				else if(w_0_2_w[k]>0 && wire_array_0_2[k] == 1 && buco==0) ++count_2;
				else if(w_0_2_w[k]>0 && wire_array_0_2[k] == 0)
				{
					if(buco == 0)start_buco = k;
					++buco;
				}
			}
			wire_hits.clear();

			//count_2 = 0;
			//for(int k = min_wire_2; k <= max_wire_2; ++k) {
			//	if(wire_array_0_2[k] == 1) ++count_2;
			//}
			//store efficiency and average pitch in the ind1 vectors
            if(wcount_2 > (max_wire_2-min_wire_2+1) || count_2 > wcount_2)cout << "IND 0 2 W " << count_2 << " " << wcount_2 << " " << min_wire_2 << " " << max_wire_2 << endl;
			if(max_buco<10)hist_eff_0.push_back(count_2/(wcount_2));
			if(max_buco<10)hist_average_pitch_0.push_back(holder_2);

			cout << "WEST CRYO, PLANE IND1, LOGIC TPC 2 ";
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE IND1, LOGIC TPC 2
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 2;
			_plane = 0;
			_min_wire = min_wire_2;
			_max_wire = max_wire_2;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_2;
			_tot_wires = wcount_2;
			_avg_pitch = holder_2;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
			_which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			if(PRINT_WEST_IND1_TPC2)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();
			
		}
				
		// TPC3 IND1	
		if(count_3 > wire_threshold) 
		{
			//average pitch
			holder_3 = holder_3/count_3;

			//find min wire
			for(int k = 0; k < 1056; ++k) 
			{
				if(wire_array_0_3[k] == 1) 
				{
					min_wire_3 = k;
					break;		
				}
			}
					
			//find max wire
			//for(int k = 1056; k > 0; --k) 
			for(int k = 1055; k > 0; --k)
			{
				if(wire_array_0_3[k] == 1) 
				{
					max_wire_3 = k;
					break;		
				}
			}

			//find total amount of wires in the trace
			count_3 = 0; buco=0; wcount_3 = 0;max_buco=0;

			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_0.GetSize(); hit++)
			{
				if(tpcTPCW_0[hit] == 3 /*&& hasSPTPCW_0[hit]*/)
				{
					if(PRINT_WEST_IND1_TPC3)
					{
						cout << wireTPCW_0[hit] << " " << hit_time_TPCW_0[hit] << " " << hasSPTPCW_0[hit] << " " << (wireTPCW_0[hit] >= min_wire_3 && wireTPCW_0[hit] <= max_wire_3) << " " << w_0_3_w[wireTPCW_0[hit]] << " ";
						if(!isnan(hitx_TPCW_0[hit])){ cout << hitx_TPCW_0[hit] << " " << hity_TPCW_0[hit] << " " << hitz_TPCW_0[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_0[hit] >= min_wire_3 && wireTPCW_0[hit] <= max_wire_3)
					{
						wire_hits.push_back({wireTPCW_0[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_0[hit],hity_TPCW_0[hit],hitz_TPCW_0[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}

			double start_buco = -1;
			for(int k = min_wire_3; k <= max_wire_3; ++k) 
			{
                wcount_3+=w_0_3_w[k];
                if(w_0_3_w[k]>0 && wire_array_0_3[k] == 1 && buco>0) 
				{
					++count_3;

					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_0_3[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_0[hh.first]);
					this_hole.push_back(hity_TPCW_0[hh.first]);
					this_hole.push_back(hitz_TPCW_0[hh.first]);
					this_hole.push_back(hitx_TPCW_0[hh.second]);
					this_hole.push_back(hity_TPCW_0[hh.second]);
					this_hole.push_back(hitz_TPCW_0[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_0[hh.first]);
					this_hole.push_back(hit_diry_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_0[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_0[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_0[hh.second]);
					this_hole.push_back(hit_time_TPCW_0[hh.first]);
					this_hole.push_back(hit_time_TPCW_0[hh.second]);
					this_hole.push_back(hit_mult_TPCW_0[hh.first]);
					this_hole.push_back(hit_mult_TPCW_0[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hind1buchi->Fill(holder_3,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;
				}
			  	else if(w_0_3_w[k]>0 && wire_array_0_3[k] == 1 && buco==0) ++count_3;
			  	else if(w_0_3_w[k]>0 && wire_array_0_3[k] == 0)
				{
					if(buco == 0)start_buco = k;
					++buco;
				}
			}
			wire_hits.clear();

			//count_3 = 0;
			//for(int k = min_wire_3; k <= max_wire_3; ++k) {
			//	if(wire_array_0_3[k] == 1) ++count_3;
			//}
			//store efficiency and average pitch in the ind1 vectors
            if(wcount_3 > (max_wire_3-min_wire_3+1) || count_3 > wcount_3)cout << "IND 0 3 W " << count_3 << " " << wcount_3 << " " << min_wire_3 << " " << max_wire_3 << endl;
			if(max_buco<10)hist_eff_0.push_back(count_3/(wcount_3));
			if(max_buco<10)hist_average_pitch_0.push_back(holder_3);

			cout << "WEST CRYO, PLANE IND1, LOGIC TPC 3 "; 
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE IND1, LOGIC TPC 3
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 3;
			_plane = 0;
			_min_wire = min_wire_3;
			_max_wire = max_wire_3;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_3;
			_tot_wires = wcount_3;
			_avg_pitch = holder_3;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
			_which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			if(PRINT_WEST_IND1_TPC3)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();
		}
				
//
//	RESET ARRAYS AND COUNTERS
//

		for(int j = 0; j < 1056; ++j) 
		{
			wire_array_0_0[j] = 0;
			wire_array_0_1[j] = 0;
			wire_array_0_2[j] = 0;
			wire_array_0_3[j] = 0;
		}

		count_0 = 0;
		count_1 = 0;
		count_2 = 0;
		count_3 = 0;	

		holder_0 = 0;	
		holder_1 = 0;
		holder_2 = 0;
		holder_3 = 0;
				
//
//	PLANE IND2 WEST
//

		//loop over all hits of the trace
		for(int j = 0; j < wireTPCW_1.GetSize(); ++j) 
		{
			//conditions to consider a valid hit, can change
			if(dqdxTPCW_1[j] > 0 && pitchTPCW_1[j] < 4 && pitchTPCW_1[j] > 0 && !isnan(pitchTPCW_1[j])) 
			{
				if(tpcTPCW_1[j] == 0) 
				{
					holder_01 += pitchTPCW_1[j];
					++count_01;
				}
				else if(tpcTPCW_1[j] == 1) 
				{
					holder_01 += pitchTPCW_1[j];
					++count_01;
				}
				else if(tpcTPCW_1[j] == 2) 
				{
					holder_23 += pitchTPCW_1[j];
					++count_23;
				}
				else 
				{
					holder_23 += pitchTPCW_1[j];
					++count_23;
				}
			}
					
			//write down valid wire hits, less restrictions
            if(hasSPTPCW_1[j])
			{
				if(tpcTPCW_1[j] == 0) wire_array_1_01[wireTPCW_1[j]] = 1;
				else if(tpcTPCW_1[j] == 1) wire_array_1_01[wireTPCW_1[j] + 2536] = 1;
				else if(tpcTPCW_1[j] == 2) wire_array_1_23[wireTPCW_1[j]] = 1;
				else wire_array_1_23[wireTPCW_1[j] + 2536] = 1;
			}				
		}

				
		//check if there are at least the minimum amount of wires in the TPCs	
		/*
			for(int k = 0; k < 5600; ++k) 
			{

				if(wire_array_1_01[k] == 1)h1_01_w->Fill(k);
				if(wire_array_1_23[k] == 1)h1_23_w->Fill(k);
	
			}
		*/	
				
		// TPC01 IND2
		if(count_01 > wire_threshold) 
		{
			//average pitch
			holder_01 = holder_01/count_01;
					
			//find min wire
			for(int k = 0; k < 5600; ++k) 
			{
				if(wire_array_1_01[k] == 1) 
				{
					min_wire_01 = k;
					break;		
				}
			}
					
			//find max wire
			//for(int k = 5600; k > 0; --k) 
			for(int k = 5599; k > 0; --k) 
			{
				if(wire_array_1_01[k] == 1) 
				{
					max_wire_01 = k;
					break;		
				}
			}
					
			//find total amount of wires in the trace
			count_01 = 0; wcount_01 = 0; buco=0;max_buco=0;

			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_1.GetSize(); hit++)
			{

				if(tpcTPCW_1[hit] == 0 /*&& hasSPTPCW_1[hit]*/)
				{
					if(PRINT_WEST_IND2_TPC01)
					{
						cout << wireTPCW_1[hit] << " " << hit_time_TPCW_1[hit] << " " << hasSPTPCW_1[hit] << " " << (wireTPCW_1[hit] >= min_wire_01 && wireTPCW_1[hit] <= max_wire_01) << " " << w_1_01_w[wireTPCW_1[hit]] << " ";
						if(!isnan(hitx_TPCW_1[hit])){ cout << hitx_TPCW_1[hit] << " " << hity_TPCW_1[hit] << " " << hitz_TPCW_1[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_1[hit] >= min_wire_01 && wireTPCW_1[hit] <= max_wire_01)
					{
						wire_hits.push_back({wireTPCW_1[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_1[hit],hity_TPCW_1[hit],hitz_TPCW_1[hit]);
						coordinates.push_back(temp_coord);
					}
				}
				else if(tpcTPCW_1[hit] == 1)
				{
					if(PRINT_WEST_IND2_TPC01)
					{
						cout << wireTPCW_1[hit] + 2536 << " " << hit_time_TPCW_1[hit] << " " << hasSPTPCW_1[hit] << " " << (wireTPCW_1[hit] + 2536 >= min_wire_01 && wireTPCW_1[hit] + 2536 <= max_wire_01) << " " << w_1_01_w[wireTPCW_1[hit] + 2536] << " ";
						if(!isnan(hitx_TPCW_1[hit])){ cout << hitx_TPCW_1[hit] << " " << hity_TPCW_1[hit] << " " << hitz_TPCW_1[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_1[hit] + 2536 >= min_wire_01 && wireTPCW_1[hit] + 2536 <= max_wire_01)
					{
						wire_hits.push_back({wireTPCW_1[hit] + 2536,hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_1[hit],hity_TPCW_1[hit],hitz_TPCW_1[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}

			double start_buco = -1;
			for(int k = min_wire_01; k <= max_wire_01; ++k) 
			{
                wcount_01+=w_1_01_w[k];

                if(w_1_01_w[k]>0 && wire_array_1_01[k] == 1 && buco>0) 
				{
					++count_01;

					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_1_01[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_1[hh.first]);
					this_hole.push_back(hity_TPCW_1[hh.first]);
					this_hole.push_back(hitz_TPCW_1[hh.first]);
					this_hole.push_back(hitx_TPCW_1[hh.second]);
					this_hole.push_back(hity_TPCW_1[hh.second]);
					this_hole.push_back(hitz_TPCW_1[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_1[hh.first]);
					this_hole.push_back(hit_diry_TPCW_1[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_1[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_1[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_1[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_1[hh.second]);
					this_hole.push_back(hit_time_TPCW_1[hh.first]);
					this_hole.push_back(hit_time_TPCW_1[hh.second]);
					this_hole.push_back(hit_mult_TPCW_1[hh.first]);
					this_hole.push_back(hit_mult_TPCW_1[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hind2buchi->Fill(holder_01,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;
				}
				else if(w_1_01_w[k]>0 && wire_array_1_01[k] == 1 && buco==0) ++count_01;
				else if(w_1_01_w[k]>0 && wire_array_1_01[k] == 0)
				{
					if(buco == 0)start_buco = k;
					++buco;
				}
			}

			wire_hits.clear();

			//count_01 = 0;
			//for(int k = min_wire_01; k <= max_wire_01; ++k) {
			//	if(wire_array_1_01[k] == 1) ++count_01;
			//}
			//store efficiency and average pitch in the ind2 vectors

            if(wcount_01 > (max_wire_01-min_wire_01+1) || count_01 > wcount_01)cout << "IND 1 01 W " << count_01 << " " << wcount_01 << " " << min_wire_01 << " " << max_wire_01 << endl;
			if(max_buco<10)hist_eff_1.push_back(count_01/(wcount_01));
			if(max_buco<10)hist_average_pitch_1.push_back(holder_01);

			cout << "WEST CRYO, PLANE IND2, LOGIC TPC 01 ";
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE IND2, LOGIC TPC 01
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 0;
			_plane = 1;
			_min_wire = min_wire_01;
			_max_wire = max_wire_01;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_01;
			_tot_wires = wcount_01;
			_avg_pitch = holder_01;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
		    _which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			if(PRINT_WEST_IND2_TPC01)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();

		}
				
		// TPC23 IND2
		if(count_23 > wire_threshold) 
		{
			//average pitch
			holder_23 = holder_23/count_23;

			//find min wire
			for(int k = 0; k < 5600; ++k) 
			{
				if(wire_array_1_23[k] == 1) 
				{
					min_wire_23 = k;
					break;		
				}
			}
					
			//find max wire
			//for(int k = 5600; k > 0; --k) 
			for(int k = 5599; k > 0; --k) 
			{
				if(wire_array_1_23[k] == 1) 
				{
					max_wire_23 = k;
					break;		
				}
			}
					
			//find total amount of wires in the trace
			count_23 = 0; wcount_23 = 0; buco=0;max_buco=0;

			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_1.GetSize(); hit++)
			{
				if(tpcTPCW_1[hit] == 2 /*&& hasSPTPCW_1[hit]*/)
				{
					if(PRINT_WEST_IND2_TPC23)
					{
						cout << wireTPCW_1[hit] << " " << hit_time_TPCW_1[hit] << " " << hasSPTPCW_1[hit] << " " << (wireTPCW_1[hit] >= min_wire_23 && wireTPCW_1[hit] <= max_wire_23) << " " << w_1_23_w[wireTPCW_1[hit]] << " ";
						if(!isnan(hitx_TPCW_1[hit])){ cout << hitx_TPCW_1[hit] << " " << hity_TPCW_1[hit] << " " << hitz_TPCW_1[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_1[hit] >= min_wire_23 && wireTPCW_1[hit] <= max_wire_23)
					{
						wire_hits.push_back({wireTPCW_1[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_1[hit],hity_TPCW_1[hit],hitz_TPCW_1[hit]);
						coordinates.push_back(temp_coord);
					}
				}
				else if(tpcTPCW_1[hit] == 3)
				{
					if(PRINT_WEST_IND2_TPC23)
					{
						cout << wireTPCW_1[hit] + 2536 << " " << hit_time_TPCW_1[hit] << " " << hasSPTPCW_1[hit] << " " << (wireTPCW_1[hit] + 2536 >= min_wire_23 && wireTPCW_1[hit] + 2536 <= max_wire_23) << " " << w_1_23_w[wireTPCW_1[hit] + 2536] << " ";
						if(!isnan(hitx_TPCW_1[hit])){ cout << hitx_TPCW_1[hit] << " " << hity_TPCW_1[hit] << " " << hitz_TPCW_1[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_1[hit] + 2536 >= min_wire_23 && wireTPCW_1[hit] + 2536 <= max_wire_23)
					{
						wire_hits.push_back({wireTPCW_1[hit] + 2536,hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_1[hit],hity_TPCW_1[hit],hitz_TPCW_1[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}

			double start_buco = -1;
			for(int k = min_wire_23; k <= max_wire_23; ++k) 
			{
                wcount_23+=w_1_23_w[k];

                if(w_1_23_w[k]>0 && wire_array_1_23[k] == 1 && buco>0) 
				{
					++count_23;

					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_1_23[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_1[hh.first]);
					this_hole.push_back(hity_TPCW_1[hh.first]);
					this_hole.push_back(hitz_TPCW_1[hh.first]);
					this_hole.push_back(hitx_TPCW_1[hh.second]);
					this_hole.push_back(hity_TPCW_1[hh.second]);
					this_hole.push_back(hitz_TPCW_1[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_1[hh.first]);
					this_hole.push_back(hit_diry_TPCW_1[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_1[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_1[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_1[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_1[hh.second]);
					this_hole.push_back(hit_time_TPCW_1[hh.first]);
					this_hole.push_back(hit_time_TPCW_1[hh.second]);
					this_hole.push_back(hit_mult_TPCW_1[hh.first]);
					this_hole.push_back(hit_mult_TPCW_1[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hind2buchi->Fill(holder_23,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;
				}
				else if(w_1_23_w[k]>0 && wire_array_1_23[k] == 1 && buco==0) ++count_23;
				else if(w_1_23_w[k]>0 && wire_array_1_23[k] == 0)
				{
					if(buco == 0)start_buco = k;
					++buco;
				}
			}

			wire_hits.clear();

			//count_23 = 0;
			//for(int k = min_wire_23; k <= max_wire_23; ++k) {
			//	if(wire_array_1_23[k] == 1) ++count_23;
			//}
			//store efficiency and average pitch in the ind1 vectors

            if(wcount_23 > (max_wire_23-min_wire_23+1) || count_23 > wcount_23)cout << "IND 1 23 W " << count_23 << " " << wcount_23 << " " << min_wire_23 << " " << max_wire_23 << endl;
			if(max_buco<10)hist_eff_1.push_back(count_23/(wcount_23));
			if(max_buco<10)hist_average_pitch_1.push_back(holder_23);

			cout << "WEST CRYO, PLANE IND2, LOGIC TPC 23 ";
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE IND2, LOGIC TPC 23
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 2;
			_plane = 1;
			_min_wire = min_wire_23;
			_max_wire = max_wire_23;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_23;
			_tot_wires = wcount_23;
			_avg_pitch = holder_23;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
			_which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			if(PRINT_WEST_IND2_TPC23)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();

		}
				
//
//	RESET ARRAYS AND COUNTERS
//				
				
		for(int j = 0; j < 5600; ++j) 
		{
			wire_array_1_01[j] = 0;
			wire_array_1_23[j] = 0;
		}

		count_01 = 0;
		count_23 = 0;	

		holder_01 = 0;	
		holder_23 = 0;

//
//	PLANE COLL WEST
//

		//loop over all hits of the trace
				
		for(int j = 0; j < wireTPCW_2.GetSize(); ++j) 
		{
			//conditions to consider a valid hit, can change

			if(dqdxTPCW_2[j] > 0 && pitchTPCW_2[j] < 4 && pitchTPCW_2[j] > 0 && !isnan(pitchTPCW_2[j])) 
			{
				if(tpcTPCW_2[j] == 0) 
				{
					holder_01 += pitchTPCW_2[j];
					++count_01;
				}
				else if(tpcTPCW_2[j] == 1) 
				{
					holder_01 += pitchTPCW_2[j];
					++count_01;
				}
				else if(tpcTPCW_2[j] == 2) 
				{
					holder_23 += pitchTPCW_2[j];
					++count_23;
				}
				else 
				{
					holder_23 += pitchTPCW_2[j];
					++count_23;
				}
			}
					
			//write down valid wire hits, less restrictions
            if(hasSPTPCW_2[j])
			{
				if(tpcTPCW_2[j] == 0) wire_array_2_01[wireTPCW_2[j]] = 1;
				else if(tpcTPCW_2[j] == 1) wire_array_2_01[wireTPCW_2[j] + 2536] = 1;
				else if(tpcTPCW_2[j] == 2) wire_array_2_23[wireTPCW_2[j]] = 1;
				else wire_array_2_23[wireTPCW_2[j] + 2536] = 1;
			}
		}

		//check if there are at least the minimum amount of wires in the TPCs	
		/*				                          
			for(int k = 0; k < 5600; ++k) 
			{
				if(wire_array_2_01[k] == 1)h2_01_w->Fill(k);
				if(wire_array_2_23[k] == 1)h2_23_w->Fill(k);
 			}
		*/

		// TPC01 COLL

		if(count_01 > wire_threshold) 
		{
			//average pitch
			holder_01 = holder_01/count_01;

			//find min wire
			for(int k = 0; k < 5600; ++k) 
			{
				if(wire_array_2_01[k] == 1) 
				{
					min_wire_01 = k;
					break;		
				}
			}
					
			//find max wire
			//for(int k = 5600; k > 0; --k) 
			for(int k = 5599; k > 0; --k) 
			{
				if(wire_array_2_01[k] == 1) 
				{
					max_wire_01 = k;
					break;		
				}
			}

			//find total amount of wires in the trace
            count_01 = 0; buco=0; wcount_01=0;max_buco=0;

			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_2.GetSize(); hit++)
			{
				if(tpcTPCW_2[hit] == 0 /*&& hasSPTPCW_2[hit]*/)
				{
					if(PRINT_WEST_COLL_TPC01)
					{
						cout << wireTPCW_2[hit] << " " << hit_time_TPCW_2[hit] << " " << hasSPTPCW_2[hit] << " " << (wireTPCW_2[hit] >= min_wire_01 && wireTPCW_2[hit] <= max_wire_01) << " " << w_2_01_w[wireTPCW_2[hit]] << " ";
						if(!isnan(hitx_TPCW_2[hit])){ cout << hitx_TPCW_2[hit] << " " << hity_TPCW_2[hit] << " " << hitz_TPCW_2[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_2[hit] >= min_wire_01 && wireTPCW_2[hit] <= max_wire_01)
					{
						wire_hits.push_back({wireTPCW_2[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_2[hit],hity_TPCW_2[hit],hitz_TPCW_2[hit]);
						coordinates.push_back(temp_coord);
					}
				}
				else if(tpcTPCW_2[hit] == 1)
				{
					if(PRINT_WEST_COLL_TPC01)
					{
						cout << wireTPCW_2[hit] + 2536 << " " << hit_time_TPCW_2[hit] << " " << hasSPTPCW_2[hit] << " " << (wireTPCW_2[hit] + 2536 >= min_wire_01 && wireTPCW_2[hit] + 2536 <= max_wire_01) << " " << w_2_01_w[wireTPCW_2[hit] + 2536] << " ";
						if(!isnan(hitx_TPCW_2[hit])){ cout << hitx_TPCW_2[hit] << " " << hity_TPCW_2[hit] << " " << hitz_TPCW_2[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_2[hit] + 2536 >= min_wire_01 && wireTPCW_2[hit] + 2536 <= max_wire_01)
					{
						wire_hits.push_back({wireTPCW_2[hit] + 2536,hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_2[hit],hity_TPCW_2[hit],hitz_TPCW_2[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}

			double start_buco = -1;
			for(int k = min_wire_01; k <= max_wire_01; ++k) 
			{
                wcount_01+=w_2_01_w[k];

                if(w_2_01_w[k]>0 && wire_array_2_01[k] == 1 && buco>0) 
				{
					++count_01;

					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_2_01[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_2[hh.first]);
					this_hole.push_back(hity_TPCW_2[hh.first]);
					this_hole.push_back(hitz_TPCW_2[hh.first]);
					this_hole.push_back(hitx_TPCW_2[hh.second]);
					this_hole.push_back(hity_TPCW_2[hh.second]);
					this_hole.push_back(hitz_TPCW_2[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_2[hh.first]);
					this_hole.push_back(hit_diry_TPCW_2[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_2[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_2[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_2[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_2[hh.second]);
					this_hole.push_back(hit_time_TPCW_2[hh.first]);
					this_hole.push_back(hit_time_TPCW_2[hh.second]);
					this_hole.push_back(hit_mult_TPCW_2[hh.first]);
					this_hole.push_back(hit_mult_TPCW_2[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hcollbuchi->Fill(holder_01,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;
				}
				else if(w_2_01_w[k]>0 && wire_array_2_01[k] == 1 && buco==0) ++count_01;
				else if(w_2_01_w[k]>0 && wire_array_2_01[k] == 0)
				{
					if(buco == 0)start_buco = k;
					++buco;
				}
			}
			wire_hits.clear();
					
			//count_01 = 0;
			//for(int k = min_wire_01; k <= max_wire_01; ++k) 
			//{
				//if(wire_array_2_01[k] == 1) ++count_01;
			//}
			//store efficiency and average pitch in the ind2 vectors

            if(wcount_01 > (max_wire_01-min_wire_01+1) || count_01 > wcount_01)cout << "Coll 2 01 W " << count_01 << " " << wcount_01 << " " << min_wire_01 << " " << max_wire_01 << endl;
			if(max_buco<10)hist_eff_2.push_back(count_01/(wcount_01));
			if(max_buco<10)hist_average_pitch_2.push_back(holder_01);

			//print details of low efficiency traces
			//if(count_01/(max_wire_01 - min_wire_01 + 1) < upper_eff_print) 
			//{
			//	oFile << "WEST 01 COLL " << *runTPCW << " " << *evtTPCW;
			//	oFile << endl << "\t";
			//	oFile << "Efficiency: " << count_01/(max_wire_01 - min_wire_01 + 1);
			//	oFile << " over " << *lengthTPCW << " cm.";
			//	oFile << endl << "\t";
			//	oFile << *startxTPCW << " " << *startyTPCW << " " << *startzTPCW << " ";
			//	oFile << endl << "\t";
			//	oFile << *endxTPCW << " " << *endyTPCW << " " << *endzTPCW << " ";
			//	oFile << endl << endl;
			//}		

			cout << "WEST CRYO, PLANE COLL, LOGIC TPC 01 ";
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE COLL, LOGIC TPC 01
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 0;
			_plane = 2;
			_min_wire = min_wire_01;
			_max_wire = max_wire_01;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_01;
			_tot_wires = wcount_01;
			_avg_pitch = holder_01;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
			_which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			if(PRINT_WEST_COLL_TPC01)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();

		}
				
		// TPC23 COLL
		if(count_23 > wire_threshold) 
		{
			//average pitch
			holder_23 = holder_23/count_23;

			//find min wire
			for(int k = 0; k < 5600; ++k) 
			{
				if(wire_array_2_23[k] == 1) 
				{
					min_wire_23 = k;
					break;		
				}
			}
					
			//find max wire
			//for(int k = 5600; k > 0; --k) 
			for(int k = 5599; k > 0; --k) 
			{
				if(wire_array_2_23[k] == 1) 
				{
					max_wire_23 = k;
					break;		
				}
			}

			//find total amount of wires in the trace
            count_23 = 0; buco=0; wcount_23=0;max_buco=0;

			std::vector<std::vector<double>> holes_temp; //--> store all holes in this logic TPC
			std::vector<std::pair<int,int>> wire_hits;
			std::vector<TVector3> coordinates;


			for(int hit = 0; hit<wireTPCW_2.GetSize(); hit++)
			{
				if(tpcTPCW_2[hit] == 2 /*&& hasSPTPCW_2[hit]*/)
				{
					if(PRINT_WEST_COLL_TPC23)
					{
						cout << wireTPCW_2[hit] << " " << hit_time_TPCW_2[hit] << " " << hasSPTPCW_2[hit] << " " << (wireTPCW_2[hit] >= min_wire_23 && wireTPCW_2[hit] <= max_wire_23) << " " << w_2_23_w[wireTPCW_2[hit]] << " ";
						if(!isnan(hitx_TPCW_2[hit])){ cout << hitx_TPCW_2[hit] << " " << hity_TPCW_2[hit] << " " << hitz_TPCW_2[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_2[hit] >= min_wire_23 && wireTPCW_2[hit] <= max_wire_23)
					{
						wire_hits.push_back({wireTPCW_2[hit],hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_2[hit],hity_TPCW_2[hit],hitz_TPCW_2[hit]);
						coordinates.push_back(temp_coord);
					}
				}
				else if(tpcTPCW_2[hit] == 3)
				{
					if(PRINT_WEST_COLL_TPC23)
					{
						cout << wireTPCW_2[hit] + 2536 << " " << hit_time_TPCW_2[hit] << " " << hasSPTPCW_2[hit] << " " << (wireTPCW_2[hit] + 2536 >= min_wire_23 && wireTPCW_2[hit] + 2536 <= max_wire_23) << " " << w_2_23_w[wireTPCW_2[hit] + 2536] << " ";
						if(!isnan(hitx_TPCW_2[hit])){ cout << hitx_TPCW_2[hit] << " " << hity_TPCW_2[hit] << " " << hitz_TPCW_2[hit] << endl;}
						else {cout << -1 << " " << -1 << " " << -1 << endl;} 
					}

					if(wireTPCW_2[hit] + 2536 >= min_wire_23 && wireTPCW_2[hit] + 2536 <= max_wire_23)
					{
						wire_hits.push_back({wireTPCW_2[hit] + 2536,hit});
						TVector3 temp_coord;
						temp_coord.SetXYZ(hitx_TPCW_2[hit],hity_TPCW_2[hit],hitz_TPCW_2[hit]);
						coordinates.push_back(temp_coord);
					}
				}
			}

			double start_buco = -1;
			for(int k = min_wire_23; k <= max_wire_23; ++k) 
			{
                wcount_23+=w_2_23_w[k];

				if(w_2_23_w[k]>0 && wire_array_2_23[k] == 1 && buco>0) 
				{
					++count_23;

					std::vector<double> this_hole;

					int wire_before_hole = start_buco - 1;
					while(wire_array_2_23[wire_before_hole] == 0){wire_before_hole = wire_before_hole -1;}
					
					std::pair<int,int> hh = wireToHit(wire_hits,coordinates,wire_before_hole,k);

					// --> get more info about holes
					this_hole.push_back(buco);
					this_hole.push_back(hitx_TPCW_2[hh.first]);
					this_hole.push_back(hity_TPCW_2[hh.first]);
					this_hole.push_back(hitz_TPCW_2[hh.first]);
					this_hole.push_back(hitx_TPCW_2[hh.second]);
					this_hole.push_back(hity_TPCW_2[hh.second]);
					this_hole.push_back(hitz_TPCW_2[hh.second]);
					this_hole.push_back(wire_before_hole);
					this_hole.push_back(k);
					this_hole.push_back(hit_dirx_TPCW_2[hh.first]);
					this_hole.push_back(hit_diry_TPCW_2[hh.first]);
					this_hole.push_back(hit_dirz_TPCW_2[hh.first]);
					this_hole.push_back(hit_dirx_TPCW_2[hh.second]); 
					this_hole.push_back(hit_diry_TPCW_2[hh.second]);
					this_hole.push_back(hit_dirz_TPCW_2[hh.second]);
					this_hole.push_back(hit_time_TPCW_2[hh.first]);
					this_hole.push_back(hit_time_TPCW_2[hh.second]);
					this_hole.push_back(hit_mult_TPCW_2[hh.first]);
					this_hole.push_back(hit_mult_TPCW_2[hh.second]);

					holes_temp.push_back(this_hole);
					this_hole.clear();
					// <--

					hcollbuchi->Fill(holder_23,buco);
					if(buco>max_buco)max_buco=buco;
					buco=0;
					start_buco = -1;

				}
				else if(w_2_23_w[k]>0 && wire_array_2_23[k] == 1 && buco==0) ++count_23;
				else if(w_2_23_w[k]>0 && wire_array_2_23[k] == 0)
				{
					if(buco == 0)start_buco = k;
					++buco;
				}
			}
			wire_hits.clear();

			//count_23 = 0;
			//for(int k = min_wire_23; k <= max_wire_23; ++k) 
			//{
			//	if(wire_array_2_23[k] == 1) ++count_23;
			//}
			//store efficiency and average pitch in the ind1 vectors
                                        
			if(wcount_23 > (max_wire_23-min_wire_23+1) || count_23 > wcount_23)cout << "IND 2 23 W " << count_23 << " " << wcount_23 << " " << min_wire_23 << " " << max_wire_23 << endl;

			if(max_buco<10)hist_eff_2.push_back(count_23/(wcount_23));
			if(max_buco<10)hist_average_pitch_2.push_back(holder_23);
			//print details of low efficiency traces
			
			//if(count_23/(max_wire_23 - min_wire_23 + 1) < upper_eff_print) {
			//	oFile << "WEST 23 COLL " << *runTPCW << " " << *evtTPCW;
			//	oFile << endl << "\t";
			//	oFile << "Efficiency: " << count_23/(max_wire_23 - min_wire_23 + 1);
			//	oFile << " over " << *lengthTPCW << " cm.";
			//	oFile << endl << "\t";
			//	oFile << *startxTPCW << " " << *startyTPCW << " " << *startzTPCW << " ";
			//	oFile << endl << "\t";
			//	oFile << *endxTPCW << " " << *endyTPCW << " " << *endzTPCW << " ";
			//	oFile << endl << endl;
			//}			

			cout << "WEST CRYO, PLANE COLL, LOGIC TPC 23 ";
			// WRITE VARIABLES IN TREE --> WEST CRYO, PLANE COLL, LOGIC TPC 23
			_run = *runTPCW;
			_evt = *evtTPCW; 
			_cryo = *cryostatTPCW;
			_tpc = 2;
			_plane = 2;
			_min_wire = min_wire_23;
			_max_wire = max_wire_23;
			_start.SetXYZ(*startxTPCW,*startyTPCW,*startzTPCW);
			_end.SetXYZ(*endxTPCW,*endyTPCW,*endzTPCW);
			_hit_wires = count_23;
			_tot_wires = wcount_23;
			_avg_pitch = holder_23;
			_trk_length = *lengthTPCW;
			_max_buco = max_buco;
			_which_t0 = *whicht0TPCW;
			_t0PFP = *t0_PFP_TPCW;
			_t0CRTTrack = *t0_CRT_Track_TPCW;
			_t0CRTHit = *t0_CRT_Hit_TPCW;
			_nholes = holes_temp.size();
			_wire_holes = holes_temp;
			_G4ID = *truthG4W;
			_truth_fraction = truth_fraction;

			if(PRINT_WEST_COLL_TPC23)
			{
				cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
				for(const auto &hole : _wire_holes)
        		{
            		cout << "hole ";
            		for(const auto &hole_feature : hole)
            		{
                		cout << hole_feature << " ";
            		}
        		}
			}

			outree -> Fill();
			cout << "filled" << endl;

			holes_temp.clear();

		}

//
// RESET ARRAYS AND COUNTERS
//

	for(int j = 0; j < 5600; ++j) 
	{
		wire_array_2_01[j] = 0;
		wire_array_2_23[j] = 0;
	}

	count_01 = 0;
	count_23 = 0;	

	holder_01 = 0;	
	holder_23 = 0;
				
				
}
}			

//
// EAST
//
 
 
myReaderTPCE.Restart();
while(myReaderTPCE.Next()) 
{	

	//if(*truthG4E == -1)continue; //MC in Overlays
	//if(*truthG4E != -1)continue; // DATA in Overlays

	//if(*truthG4E == -1)continue; // --> MC in Overlays // --> HERE

	int bestplane = -1;
	std::array<float,3> nhit = {0,0,0};
	std::array<float,3> n_valid_hit = {0,0,0};
	float truth_fraction = 0;

	for(int hit = 0; hit < tpcTPCE_0.GetSize(); hit++)
	{
		nhit[0] ++;
		if(h0_truth_ne_E[hit] > 0 && h0_truth_e_E[hit] > 0){n_valid_hit[0]++;}
	}
	for(int hit = 0; hit < tpcTPCE_1.GetSize(); hit++)
	{
		nhit[1] ++;
		if(h1_truth_ne_E[hit] > 0 && h1_truth_e_E[hit] > 0){n_valid_hit[1]++;}
	}
	for(int hit = 0; hit < tpcTPCE_2.GetSize(); hit++)
	{
		nhit[2] ++;
		if(h2_truth_ne_E[hit] > 0 && h2_truth_e_E[hit] > 0){n_valid_hit[2]++;}
	}

	bestplane = 0;
	for (int i = 1; i < 3; ++i) 
	{
    	if (nhit[i] > nhit[bestplane])
        bestplane = i;
	}

	if (nhit[bestplane] > 0){truth_fraction = n_valid_hit[bestplane] / nhit[bestplane];}

	//if(truth_fraction < 0.8)continue; // --> HERE

	if(false)
	{
		cout << "EAST *** NEW TRACK, G4ID: " << *truthG4E << ", TRUTH FRACTION: " << truth_fraction << endl << endl;
		cout << "plane0 " << tpcTPCE_0.GetSize() << " " << nhit[0] << " | plane1 " << tpcTPCE_1.GetSize() << " " << nhit[1] << " | plane2 " << tpcTPCE_2.GetSize() << " " << nhit[2] << " " << endl << endl;
		for(int hit = 0; hit < h2_truth_ne_E.GetSize(); hit++)
		{
			cout << h2_truth_ne_E[hit] << " ";
		}
		cout << endl;
		for(int hit = 0; hit < h2_truth_e_E.GetSize(); hit++)
		{
			cout << h2_truth_e_E[hit] << " ";
		}
		cout << endl << endl;
	}

//muon length must be greater than 1 meter
if(*lengthTPCE > 100  && ( (*whicht0TPCE)==0 || (*whicht0TPCE)==2)) 
{
                        
	//outfiletracksE << (*runTPCE) << " " << (*evtTPCE) <<  " 0 " << (*startxTPCE) << " " << (*startyTPCE) << " " << (*startzTPCE) << " " << (*endxTPCE) << " " << (*endyTPCE) << " " << (*endzTPCE) << " " << (*lengthTPCE) << endl;
 
//
//	PLANE IND1 EAST
//
				
	//loop over all hits of the trace
	for(int j = 0; j < wireTPCE_0.GetSize(); ++j) 
	{
					
		//conditions to consider a valid pitch hit
		if(dqdxTPCE_0[j] > 0 && pitchTPCE_0[j] < 4 && pitchTPCE_0[j] > 0 && !isnan(pitchTPCE_0[j])) 
		{
						
			if(tpcTPCE_0[j] == 0) 
			{
				holder_0 += pitchTPCE_0[j];
				++count_0;
			}
			else if(tpcTPCE_0[j] == 1) 
			{
				holder_1 += pitchTPCE_0[j];
				++count_1;
			}
			else if(tpcTPCE_0[j] == 2) 
			{
				holder_2 += pitchTPCE_0[j];
				++count_2;
			}
			else 
			{
				holder_3 += pitchTPCE_0[j];
				++count_3;
			}
		}
					
		//write down valid wire hits, less restrictions
        if(hasSPTPCE_0[j])
		{
			if(tpcTPCE_0[j] == 0) wire_array_0_0[wireTPCE_0[j]] = 1;
			else if(tpcTPCE_0[j] == 1) wire_array_0_1[wireTPCE_0[j]] = 1;
			else if(tpcTPCE_0[j] == 2) wire_array_0_2[wireTPCE_0[j]] = 1;
			else wire_array_0_3[wireTPCE_0[j]] = 1;
		}
	}
 
	// TPC0 IND1
	if(count_0 > wire_threshold) 
	{
		//average pitch
		holder_0 = holder_0/count_0;
 
		//find min wire
		for(int k = 0; k < 1056; ++k) 
		{
			if(wire_array_0_0[k] == 1) 
			{
				min_wire_0 = k;
				break;		
			}
		}
					
		//find max wire
		//for(int k = 1056; k > 0; --k) 
		for(int k = 1055; k > 0; --k)
		{
			if(wire_array_0_0[k] == 1) 
			{
				max_wire_0 = k;
				break;		
			}
		}
				
		//find total amount of wires in the trace
        count_0 = 0; buco=0; wcount_0=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_0.GetSize(); hit++)
		{
			if(tpcTPCE_0[hit] == 0 /*&& hasSPTPCE_0[hit]*/)
			{
				if(PRINT_EAST_IND1_TPC0)
				{
					cout << wireTPCE_0[hit] << " " << hit_time_TPCE_0[hit] << " " << hasSPTPCE_0[hit] << " " << (wireTPCE_0[hit] >= min_wire_0 && wireTPCE_0[hit] <= max_wire_0) << " " << w_0_0_e[wireTPCE_0[hit]] << " ";
					if(!isnan(hitx_TPCE_0[hit])){ cout << hitx_TPCE_0[hit] << " " << hity_TPCE_0[hit] << " " << hitz_TPCE_0[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_0[hit] >= min_wire_0 && wireTPCE_0[hit] <= max_wire_0)
				{
					wire_hits.push_back({wireTPCE_0[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_0[hit], hity_TPCE_0[hit], hitz_TPCE_0[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_0; k <= max_wire_0; ++k) 
		{ 
            wcount_0 += w_0_0_e[k];
            if(w_0_0_e[k]>0 && wire_array_0_0[k] == 1 && buco>0) 
			{
				++count_0; 
				
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_0_0[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_0[hh.first]);
				this_hole.push_back(hity_TPCE_0[hh.first]);
				this_hole.push_back(hitz_TPCE_0[hh.first]);
				this_hole.push_back(hitx_TPCE_0[hh.second]);
				this_hole.push_back(hity_TPCE_0[hh.second]);
				this_hole.push_back(hitz_TPCE_0[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_0[hh.first]);
				this_hole.push_back(hit_diry_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_0[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_0[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.second]);
				this_hole.push_back(hit_time_TPCE_0[hh.first]);
				this_hole.push_back(hit_time_TPCE_0[hh.second]);
				this_hole.push_back(hit_mult_TPCE_0[hh.first]);
				this_hole.push_back(hit_mult_TPCE_0[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hind1buchi->Fill(holder_0, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
			else if(w_0_0_e[k]>0 && wire_array_0_0[k] == 1 && buco==0) ++count_0;
			else if(w_0_0_e[k]>0 && wire_array_0_0[k] == 0)
			{ 
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_0 > (max_wire_0-min_wire_0+1) || count_0 > wcount_0) cout << "IND 0 0 E " << count_0 << " " << wcount_0 << " " << min_wire_0 << " " << max_wire_0 << endl;
		if(max_buco<10) hist_eff_0.push_back(count_0/(wcount_0));
		if(max_buco<10) hist_average_pitch_0.push_back(holder_0);
 
		cout << "EAST CRYO, PLANE IND1, LOGIC TPC 0 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE IND1, LOGIC TPC 0
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 0;
		_plane = 0;
		_min_wire = min_wire_0;
		_max_wire = max_wire_0;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_0;
		_tot_wires = wcount_0;
		_avg_pitch = holder_0;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_IND1_TPC0)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
 
	}
				
	// TPC1 IND1
	if(count_1 > wire_threshold) 
	{
		//average pitch
		holder_1 = holder_1/count_1;
 
		//find min wire
		for(int k = 0; k < 1056; ++k) 
		{
			if(wire_array_0_1[k] == 1) 
			{
				min_wire_1 = k;
				break;		
			}
		}
 
		//find max wire
		//for(int k = 1056; k > 0; --k) 
		for(int k = 1055; k > 0; --k) 
		{
			if(wire_array_0_1[k] == 1) 
			{
				max_wire_1 = k;
				break;		
			}
		}
 
		//find total amount of wires in the trace
        count_1 = 0; buco=0; wcount_1=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_0.GetSize(); hit++)
		{
			if(tpcTPCE_0[hit] == 1 /*&& hasSPTPCE_0[hit]*/)
			{
				if(PRINT_EAST_IND1_TPC1)
				{
					cout << wireTPCE_0[hit] << " " << hit_time_TPCE_0[hit] << " " << hasSPTPCE_0[hit] << " " << (wireTPCE_0[hit] >= min_wire_1 && wireTPCE_0[hit] <= max_wire_1) << " " << w_0_1_e[wireTPCE_0[hit]] << " ";
					if(!isnan(hitx_TPCE_0[hit])){ cout << hitx_TPCE_0[hit] << " " << hity_TPCE_0[hit] << " " << hitz_TPCE_0[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_0[hit] >= min_wire_1 && wireTPCE_0[hit] <= max_wire_1)
				{
					wire_hits.push_back({wireTPCE_0[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_0[hit], hity_TPCE_0[hit], hitz_TPCE_0[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_1; k <= max_wire_1; ++k) 
		{
            wcount_1 += w_0_1_e[k];
            if(w_0_1_e[k]>0 && wire_array_0_1[k] == 1 && buco>0) 
			{
				++count_1;
 
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_0_1[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_0[hh.first]);
				this_hole.push_back(hity_TPCE_0[hh.first]);
				this_hole.push_back(hitz_TPCE_0[hh.first]);
				this_hole.push_back(hitx_TPCE_0[hh.second]);
				this_hole.push_back(hity_TPCE_0[hh.second]);
				this_hole.push_back(hitz_TPCE_0[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_0[hh.first]);
				this_hole.push_back(hit_diry_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_0[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_0[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.second]);
				this_hole.push_back(hit_time_TPCE_0[hh.first]);
				this_hole.push_back(hit_time_TPCE_0[hh.second]);
				this_hole.push_back(hit_mult_TPCE_0[hh.first]);
				this_hole.push_back(hit_mult_TPCE_0[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hind1buchi->Fill(holder_1, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
		  	else if(w_0_1_e[k]>0 && wire_array_0_1[k] == 1 && buco==0) ++count_1;
		  	else if(w_0_1_e[k]>0 && wire_array_0_1[k] == 0)
			{
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_1 > (max_wire_1-min_wire_1+1) || count_1 > wcount_1) cout << "IND 0 1 E " << count_1 << " " << wcount_1 << " " << min_wire_1 << " " << max_wire_1 << endl;
		if(max_buco<10) hist_eff_0.push_back(count_1/(wcount_1));
		if(max_buco<10) hist_average_pitch_0.push_back(holder_1);
 
		cout << "EAST CRYO, PLANE IND1, LOGIC TPC 1 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE IND1, LOGIC TPC 1
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 1;
		_plane = 0;
		_min_wire = min_wire_1;
		_max_wire = max_wire_1;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_1;
		_tot_wires = wcount_1;
		_avg_pitch = holder_1;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_IND1_TPC1)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
 
	}
				
	// TPC2 IND1	
	if(count_2 > wire_threshold) 
	{
		//average pitch
		holder_2 = holder_2/count_2;
 
		//find min wire
		for(int k = 0; k < 1056; ++k) 
		{
			if(wire_array_0_2[k] == 1) 
			{
				min_wire_2 = k;
				break;		
			}
		}
				
		//find max wire
		//for(int k = 1056; k > 0; --k) 
		for(int k = 1055; k > 0; --k) 
		{
			if(wire_array_0_2[k] == 1) 
			{
				max_wire_2 = k;
				break;		
			}
		}
 
		//find total amount of wires in the trace
        count_2 = 0; buco=0; wcount_2=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_0.GetSize(); hit++)
		{
			if(tpcTPCE_0[hit] == 2 /*&& hasSPTPCE_0[hit]*/)
			{
				if(PRINT_EAST_IND1_TPC2)
				{
					cout << wireTPCE_0[hit] << " " << hit_time_TPCE_0[hit] << " " << hasSPTPCE_0[hit] << " " << (wireTPCE_0[hit] >= min_wire_2 && wireTPCE_0[hit] <= max_wire_2) << " " << w_0_2_e[wireTPCE_0[hit]] << " ";
					if(!isnan(hitx_TPCE_0[hit])){ cout << hitx_TPCE_0[hit] << " " << hity_TPCE_0[hit] << " " << hitz_TPCE_0[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_0[hit] >= min_wire_2 && wireTPCE_0[hit] <= max_wire_2)
				{
					wire_hits.push_back({wireTPCE_0[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_0[hit], hity_TPCE_0[hit], hitz_TPCE_0[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_2; k <= max_wire_2; ++k) 
		{
            wcount_2 += w_0_2_e[k];
			if(w_0_2_e[k]>0 && wire_array_0_2[k] == 1 && buco>0) 
			{
				++count_2;
 
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_0_2[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_0[hh.first]);
				this_hole.push_back(hity_TPCE_0[hh.first]);
				this_hole.push_back(hitz_TPCE_0[hh.first]);
				this_hole.push_back(hitx_TPCE_0[hh.second]);
				this_hole.push_back(hity_TPCE_0[hh.second]);
				this_hole.push_back(hitz_TPCE_0[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_0[hh.first]);
				this_hole.push_back(hit_diry_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_0[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_0[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.second]);
				this_hole.push_back(hit_time_TPCE_0[hh.first]);
				this_hole.push_back(hit_time_TPCE_0[hh.second]);
				this_hole.push_back(hit_mult_TPCE_0[hh.first]);
				this_hole.push_back(hit_mult_TPCE_0[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hind1buchi->Fill(holder_2, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
			else if(w_0_2_e[k]>0 && wire_array_0_2[k] == 1 && buco==0) ++count_2;
			else if(w_0_2_e[k]>0 && wire_array_0_2[k] == 0)
			{
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_2 > (max_wire_2-min_wire_2+1) || count_2 > wcount_2) cout << "IND 0 2 E " << count_2 << " " << wcount_2 << " " << min_wire_2 << " " << max_wire_2 << endl;
		if(max_buco<10) hist_eff_0.push_back(count_2/(wcount_2));
		if(max_buco<10) hist_average_pitch_0.push_back(holder_2);
 
		cout << "EAST CRYO, PLANE IND1, LOGIC TPC 2 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE IND1, LOGIC TPC 2
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 2;
		_plane = 0;
		_min_wire = min_wire_2;
		_max_wire = max_wire_2;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_2;
		_tot_wires = wcount_2;
		_avg_pitch = holder_2;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_IND1_TPC2)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
		
	}
				
	// TPC3 IND1	
	if(count_3 > wire_threshold) 
	{
		//average pitch
		holder_3 = holder_3/count_3;
 
		//find min wire
		for(int k = 0; k < 1056; ++k) 
		{
			if(wire_array_0_3[k] == 1) 
			{
				min_wire_3 = k;
				break;		
			}
		}
 
		//find max wire
		//for(int k = 1056; k > 0; --k) 
		for(int k = 1055; k > 0; --k) 
		{
			if(wire_array_0_3[k] == 1) 
			{
				max_wire_3 = k;
				break;		
			}
		}
 
		//find total amount of wires in the trace
        count_3 = 0; buco=0; wcount_3=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_0.GetSize(); hit++)
		{
			if(tpcTPCE_0[hit] == 3 /*&& hasSPTPCE_0[hit]*/)
			{
				if(PRINT_EAST_IND1_TPC3)
				{
					cout << wireTPCE_0[hit] << " " << hit_time_TPCE_0[hit] << " " << hasSPTPCE_0[hit] << " " << (wireTPCE_0[hit] >= min_wire_3 && wireTPCE_0[hit] <= max_wire_3) << " " << w_0_3_e[wireTPCE_0[hit]] << " ";
					if(!isnan(hitx_TPCE_0[hit])){ cout << hitx_TPCE_0[hit] << " " << hity_TPCE_0[hit] << " " << hitz_TPCE_0[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_0[hit] >= min_wire_3 && wireTPCE_0[hit] <= max_wire_3)
				{
					wire_hits.push_back({wireTPCE_0[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_0[hit], hity_TPCE_0[hit], hitz_TPCE_0[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_3; k <= max_wire_3; ++k) 
		{
        	wcount_3 += w_0_3_e[k];
		  	if(w_0_3_e[k]>0 && wire_array_0_3[k] == 1 && buco>0) 
			{
				++count_3;
 
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_0_3[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_0[hh.first]);
				this_hole.push_back(hity_TPCE_0[hh.first]);
				this_hole.push_back(hitz_TPCE_0[hh.first]);
				this_hole.push_back(hitx_TPCE_0[hh.second]);
				this_hole.push_back(hity_TPCE_0[hh.second]);
				this_hole.push_back(hitz_TPCE_0[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_0[hh.first]);
				this_hole.push_back(hit_diry_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_0[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_0[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_0[hh.second]);
				this_hole.push_back(hit_time_TPCE_0[hh.first]);
				this_hole.push_back(hit_time_TPCE_0[hh.second]);
				this_hole.push_back(hit_mult_TPCE_0[hh.first]);
				this_hole.push_back(hit_mult_TPCE_0[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hind1buchi->Fill(holder_3, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
		  	else if(w_0_3_e[k]>0 && wire_array_0_3[k] == 1 && buco==0) ++count_3;
		  	else if(w_0_3_e[k]>0 && wire_array_0_3[k] == 0)
			{
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_3 > (max_wire_3-min_wire_3+1) || count_3 > wcount_3) cout << "IND 0 3 E " << count_3 << " " << wcount_3 << " " << min_wire_3 << " " << max_wire_3 << endl;
		if(max_buco<10) hist_eff_0.push_back(count_3/(wcount_3));
		if(max_buco<10) hist_average_pitch_0.push_back(holder_3);
 
		cout << "EAST CRYO, PLANE IND1, LOGIC TPC 3 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE IND1, LOGIC TPC 3
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 3;
		_plane = 0;
		_min_wire = min_wire_3;
		_max_wire = max_wire_3;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_3;
		_tot_wires = wcount_3;
		_avg_pitch = holder_3;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_IND1_TPC3)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
 
	}
 
//	
// RESET ARRAYS AND COUNTERS
//	
 
	for(int j = 0; j < 1056; ++j) 
	{
		wire_array_0_0[j] = 0;
		wire_array_0_1[j] = 0;
		wire_array_0_2[j] = 0;
		wire_array_0_3[j] = 0;
	}
 
	count_0 = 0;
	count_1 = 0;
	count_2 = 0;
	count_3 = 0;	
 
	holder_0 = 0;	
	holder_1 = 0;
	holder_2 = 0;
	holder_3 = 0;
				
//
//	PLANE IND2 EAST
//
 
	//loop over all hits of the trace
	for(int j = 0; j < wireTPCE_1.GetSize(); ++j) 
	{
		//conditions to consider a valid hit, can change
		if(dqdxTPCE_1[j] > 0 && pitchTPCE_1[j] < 4 && pitchTPCE_1[j] > 0 && !isnan(pitchTPCE_1[j]))
		{
			if(tpcTPCE_1[j] == 0) 
			{
				holder_01 += pitchTPCE_1[j];
				++count_01;
			}
			else if(tpcTPCE_1[j] == 1) 
			{
				holder_01 += pitchTPCE_1[j];				
				++count_01;
			}
			else if(tpcTPCE_1[j] == 2) 
			{
				holder_23 += pitchTPCE_1[j];
				++count_23;
			}
			else 
			{
				holder_23 += pitchTPCE_1[j];
				++count_23;
			}
		}
					
		//write down valid wire hits, less restrictions
		if(hasSPTPCE_1[j])
		{
			if(tpcTPCE_1[j] == 0) wire_array_1_01[wireTPCE_1[j]] = 1;
			else if(tpcTPCE_1[j] == 1) wire_array_1_01[wireTPCE_1[j] + 2536] = 1;
			else if(tpcTPCE_1[j] == 2) wire_array_1_23[wireTPCE_1[j]] = 1;
			else wire_array_1_23[wireTPCE_1[j] + 2536] = 1;
		}
	}
 
	// TPC01 IND2
	if(count_01 > wire_threshold) 
	{
		//average pitch
		holder_01 = holder_01/count_01;
 
		//find min wire
		for(int k = 0; k < 5600; ++k) 
		{
			if(wire_array_1_01[k] == 1) 
			{
				min_wire_01 = k;
				break;		
			}
		}
		//find max wire
		//for(int k = 5600; k > 0; --k) 
		for(int k = 5599; k > 0; --k) 
		{
			if(wire_array_1_01[k] == 1) 
			{
				max_wire_01 = k;
				break;		
			}
		}
 
		//find total amount of wires in the trace
        count_01 = 0; buco=0; wcount_01=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_1.GetSize(); hit++)
		{
			if(tpcTPCE_1[hit] == 0 /*&& hasSPTPCE_1[hit]*/)
			{
				if(PRINT_EAST_IND2_TPC01)
				{
					cout << wireTPCE_1[hit] << " " << hit_time_TPCE_1[hit] << " " << hasSPTPCE_1[hit] << " " << (wireTPCE_1[hit] >= min_wire_01 && wireTPCE_1[hit] <= max_wire_01) << " " << w_1_01_e[wireTPCE_1[hit]] << " ";
					if(!isnan(hitx_TPCE_1[hit])){ cout << hitx_TPCE_1[hit] << " " << hity_TPCE_1[hit] << " " << hitz_TPCE_1[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_1[hit] >= min_wire_01 && wireTPCE_1[hit] <= max_wire_01)
				{
					wire_hits.push_back({wireTPCE_1[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_1[hit], hity_TPCE_1[hit], hitz_TPCE_1[hit]);
					coordinates.push_back(temp_coord);
				}
			}
			else if(tpcTPCE_1[hit] == 1)
			{
				if(PRINT_EAST_IND2_TPC01)
				{
					cout << wireTPCE_1[hit] + 2536 << " " << hit_time_TPCE_1[hit] << " " << hasSPTPCE_1[hit] << " " << (wireTPCE_1[hit] + 2536 >= min_wire_01 && wireTPCE_1[hit] + 2536 <= max_wire_01) << " " << w_1_01_e[wireTPCE_1[hit] + 2536] << " ";
					if(!isnan(hitx_TPCE_1[hit])){ cout << hitx_TPCE_1[hit] << " " << hity_TPCE_1[hit] << " " << hitz_TPCE_1[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_1[hit] + 2536 >= min_wire_01 && wireTPCE_1[hit] + 2536 <= max_wire_01)
				{
					wire_hits.push_back({wireTPCE_1[hit] + 2536, hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_1[hit], hity_TPCE_1[hit], hitz_TPCE_1[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_01; k <= max_wire_01; ++k) 
		{
            wcount_01 += w_1_01_e[k];
		  	if(w_1_01_e[k]>0 && wire_array_1_01[k] == 1 && buco>0) 
			{
				++count_01;
 
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_1_01[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_1[hh.first]);
				this_hole.push_back(hity_TPCE_1[hh.first]);
				this_hole.push_back(hitz_TPCE_1[hh.first]);
				this_hole.push_back(hitx_TPCE_1[hh.second]);
				this_hole.push_back(hity_TPCE_1[hh.second]);
				this_hole.push_back(hitz_TPCE_1[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_1[hh.first]);
				this_hole.push_back(hit_diry_TPCE_1[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_1[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_1[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_1[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_1[hh.second]);
				this_hole.push_back(hit_time_TPCE_1[hh.first]);
				this_hole.push_back(hit_time_TPCE_1[hh.second]);
				this_hole.push_back(hit_mult_TPCE_1[hh.first]);
				this_hole.push_back(hit_mult_TPCE_1[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hind2buchi->Fill(holder_01, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
		  	else if(w_1_01_e[k]>0 && wire_array_1_01[k] == 1 && buco==0) ++count_01;
		  	else if(w_1_01_e[k]>0 && wire_array_1_01[k] == 0)
			{
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_01 > (max_wire_01-min_wire_01+1) || count_01 > wcount_01) cout << "IND 1 01 E " << count_01 << " " << wcount_01 << " " << min_wire_01 << " " << max_wire_01 << endl;
		if(max_buco<10) hist_eff_1.push_back(count_01/(wcount_01));
		if(max_buco<10) hist_average_pitch_1.push_back(holder_01);
 
		cout << "EAST CRYO, PLANE IND2, LOGIC TPC 01 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE IND2, LOGIC TPC 01
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 0;
		_plane = 1;
		_min_wire = min_wire_01;
		_max_wire = max_wire_01;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_01;
		_tot_wires = wcount_01;
		_avg_pitch = holder_01;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_IND2_TPC01)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
 
	}
				
	// TPC23 IND2
	if(count_23 > wire_threshold) 
	{
		//average pitch
		holder_23 = holder_23/count_23;
 
		//find min wire
		for(int k = 0; k < 5600; ++k) 
		{
			if(wire_array_1_23[k] == 1) 
			{
				min_wire_23 = k;
				break;		
			}
		}
 
		//find max wire
		//for(int k = 5600; k > 0; --k) 
		for(int k = 5599; k > 0; --k) 
		{
			if(wire_array_1_23[k] == 1) 
			{
				max_wire_23 = k;
				break;		
			}
		}
 
		//find total amount of wires in the trace
        count_23 = 0; buco=0; wcount_23=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_1.GetSize(); hit++)
		{
			if(tpcTPCE_1[hit] == 2 /*&& hasSPTPCE_1[hit]*/)
			{
				if(PRINT_EAST_IND2_TPC23)
				{
					cout << wireTPCE_1[hit] << " " << hit_time_TPCE_1[hit] << " " << hasSPTPCE_1[hit] << " " << (wireTPCE_1[hit] >= min_wire_23 && wireTPCE_1[hit] <= max_wire_23) << " " << w_1_23_e[wireTPCE_1[hit]] << " ";
					if(!isnan(hitx_TPCE_1[hit])){ cout << hitx_TPCE_1[hit] << " " << hity_TPCE_1[hit] << " " << hitz_TPCE_1[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_1[hit] >= min_wire_23 && wireTPCE_1[hit] <= max_wire_23)
				{
					wire_hits.push_back({wireTPCE_1[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_1[hit], hity_TPCE_1[hit], hitz_TPCE_1[hit]);
					coordinates.push_back(temp_coord);
				}
			}
			else if(tpcTPCE_1[hit] == 3)
			{
				if(PRINT_EAST_IND2_TPC23)
				{
					cout << wireTPCE_1[hit] + 2536 << " " << hit_time_TPCE_1[hit] << " " << hasSPTPCE_1[hit] << " " << (wireTPCE_1[hit] + 2536 >= min_wire_23 && wireTPCE_1[hit] + 2536 <= max_wire_23) << " " << w_1_23_e[wireTPCE_1[hit] + 2536] << " ";
					if(!isnan(hitx_TPCE_1[hit])){ cout << hitx_TPCE_1[hit] << " " << hity_TPCE_1[hit] << " " << hitz_TPCE_1[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_1[hit] + 2536 >= min_wire_23 && wireTPCE_1[hit] + 2536 <= max_wire_23)
				{
					wire_hits.push_back({wireTPCE_1[hit] + 2536, hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_1[hit], hity_TPCE_1[hit], hitz_TPCE_1[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_23; k <= max_wire_23; ++k) 
		{
            wcount_23 += w_1_23_e[k];
		  	if(w_1_23_e[k]>0 && wire_array_1_23[k] == 1 && buco>0) 
			{
				++count_23;
 
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_1_23[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_1[hh.first]);
				this_hole.push_back(hity_TPCE_1[hh.first]);
				this_hole.push_back(hitz_TPCE_1[hh.first]);
				this_hole.push_back(hitx_TPCE_1[hh.second]);
				this_hole.push_back(hity_TPCE_1[hh.second]);
				this_hole.push_back(hitz_TPCE_1[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_1[hh.first]);
				this_hole.push_back(hit_diry_TPCE_1[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_1[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_1[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_1[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_1[hh.second]);
				this_hole.push_back(hit_time_TPCE_1[hh.first]);
				this_hole.push_back(hit_time_TPCE_1[hh.second]);
				this_hole.push_back(hit_mult_TPCE_1[hh.first]);
				this_hole.push_back(hit_mult_TPCE_1[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hind2buchi->Fill(holder_23, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
		  	else if(w_1_23_e[k]>0 && wire_array_1_23[k] == 1 && buco==0) ++count_23;
		  	else if(w_1_23_e[k]>0 && wire_array_1_23[k] == 0)
			{
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_23 > (max_wire_23-min_wire_23+1) || count_23 > wcount_23) cout << "IND 1 23 E " << count_23 << " " << wcount_23 << " " << min_wire_23 << " " << max_wire_23 << endl;
		if(max_buco<10) hist_eff_1.push_back(count_23/(wcount_23));
		if(max_buco<10) hist_average_pitch_1.push_back(holder_23);
 
		cout << "EAST CRYO, PLANE IND2, LOGIC TPC 23 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE IND2, LOGIC TPC 23
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 2;
		_plane = 1;
		_min_wire = min_wire_23;
		_max_wire = max_wire_23;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_23;
		_tot_wires = wcount_23;
		_avg_pitch = holder_23;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_IND2_TPC23)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
 
	}
				
//
//RESET ARRAYS AND COUNTERS
//
 
	for(int j = 0; j < 5600; ++j) 
	{
		wire_array_1_01[j] = 0;
		wire_array_1_23[j] = 0;
	}
 
	count_01 = 0;
	count_23 = 0;	
 
	holder_01 = 0;	
	holder_23 = 0;
 
//
//	PLANE COLL EAST
//
 
	//loop over all hits of the trace
	for(int j = 0; j < wireTPCE_2.GetSize(); ++j) 
	{
		//conditions to consider a valid hit, can change
		if(dqdxTPCE_2[j] > 0 && pitchTPCE_2[j] < 4 && pitchTPCE_2[j] > 0 && !isnan(pitchTPCE_2[j])) 
		{
			if(tpcTPCE_2[j] == 0) 
			{
				holder_01 += pitchTPCE_2[j];
				++count_01;
			}
			else if(tpcTPCE_2[j] == 1) 
			{
				holder_01 += pitchTPCE_2[j];
				++count_01;
			}
			else if(tpcTPCE_2[j] == 2) 
			{
				holder_23 += pitchTPCE_2[j];
				++count_23;
			}
			else 
			{
				holder_23 += pitchTPCE_2[j];	
				++count_23;
			}
		}
		
		//write down valid wire hits, less restrictions
        if(hasSPTPCE_2[j])
		{
			if(tpcTPCE_2[j] == 0) wire_array_2_01[wireTPCE_2[j]] = 1;
			else if(tpcTPCE_2[j] == 1) wire_array_2_01[wireTPCE_2[j] + 2536] = 1;
			else if(tpcTPCE_2[j] == 2) wire_array_2_23[wireTPCE_2[j]] = 1;
			else wire_array_2_23[wireTPCE_2[j] + 2536] = 1;
		}
	}
 
	// TPC01 COLL
	if(count_01 > wire_threshold) 
	{
		//average pitch
		holder_01 = holder_01/count_01;
		
		//find min wire
		for(int k = 0; k < 5600; ++k) 
		{
			if(wire_array_2_01[k] == 1) 
			{
				min_wire_01 = k;
				break;		
			}
		}
 
		//find max wire
		//for(int k = 5600; k > 0; --k) 
		for(int k = 5599; k > 0; --k) 
		{
			if(wire_array_2_01[k] == 1) 
			{
				max_wire_01 = k;
				break;		
			}
		}
 
		//find total amount of wires in the trace
        count_01 = 0; buco=0; wcount_01=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_2.GetSize(); hit++)
		{
			if(tpcTPCE_2[hit] == 0 /*&& hasSPTPCE_2[hit]*/)
			{
				if(PRINT_EAST_COLL_TPC01)
				{
					cout << wireTPCE_2[hit] << " " << hit_time_TPCE_2[hit] << " " << hasSPTPCE_2[hit] << " " << (wireTPCE_2[hit] >= min_wire_01 && wireTPCE_2[hit] <= max_wire_01) << " " << w_2_01_e[wireTPCE_2[hit]] << " ";
					if(!isnan(hitx_TPCE_2[hit])){ cout << hitx_TPCE_2[hit] << " " << hity_TPCE_2[hit] << " " << hitz_TPCE_2[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_2[hit] >= min_wire_01 && wireTPCE_2[hit] <= max_wire_01)
				{
					wire_hits.push_back({wireTPCE_2[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_2[hit], hity_TPCE_2[hit], hitz_TPCE_2[hit]);
					coordinates.push_back(temp_coord);
				}
			}
			else if(tpcTPCE_2[hit] == 1)
			{
				if(PRINT_EAST_COLL_TPC01)
				{
					cout << wireTPCE_2[hit] + 2536 << " " << hit_time_TPCE_2[hit] << " " << hasSPTPCE_2[hit] << " " << (wireTPCE_2[hit] + 2536 >= min_wire_01 && wireTPCE_2[hit] + 2536 <= max_wire_01) << " " << w_2_01_e[wireTPCE_2[hit] + 2536] << " ";
					if(!isnan(hitx_TPCE_2[hit])){ cout << hitx_TPCE_2[hit] << " " << hity_TPCE_2[hit] << " " << hitz_TPCE_2[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_2[hit] + 2536 >= min_wire_01 && wireTPCE_2[hit] + 2536 <= max_wire_01)
				{
					wire_hits.push_back({wireTPCE_2[hit] + 2536, hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_2[hit], hity_TPCE_2[hit], hitz_TPCE_2[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_01; k <= max_wire_01; ++k) 
		{
        	wcount_01 += w_2_01_e[k];
		  	if(w_2_01_e[k]>0 && wire_array_2_01[k] == 1 && buco>0)
			{
				++count_01;
 
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_2_01[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_2[hh.first]);
				this_hole.push_back(hity_TPCE_2[hh.first]);
				this_hole.push_back(hitz_TPCE_2[hh.first]);
				this_hole.push_back(hitx_TPCE_2[hh.second]);
				this_hole.push_back(hity_TPCE_2[hh.second]);
				this_hole.push_back(hitz_TPCE_2[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_2[hh.first]);
				this_hole.push_back(hit_diry_TPCE_2[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_2[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_2[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_2[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_2[hh.second]);
				this_hole.push_back(hit_time_TPCE_2[hh.first]);
				this_hole.push_back(hit_time_TPCE_2[hh.second]);
				this_hole.push_back(hit_mult_TPCE_2[hh.first]);
				this_hole.push_back(hit_mult_TPCE_2[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hcollbuchi->Fill(holder_01, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
		  	else if(w_2_01_e[k]>0 && wire_array_2_01[k] == 1 && buco==0) ++count_01;
		  	else if(w_2_01_e[k]>0 && wire_array_2_01[k] == 0)
			{
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_01 > (max_wire_01-min_wire_01+1) || count_01 > wcount_01) cout << "Coll 2 01 E " << count_01 << " " << wcount_01 << " " << min_wire_01 << " " << max_wire_01 << endl;
		if(max_buco<10) hist_eff_2.push_back(count_01/(wcount_01));
		if(max_buco<10) hist_average_pitch_2.push_back(holder_01);
		
		//print details of low efficiency traces
		//if(count_01/(max_wire_01 - min_wire_01 + 1) < upper_eff_print) 
		//{
		//	oFile << "EAST 01 COLL " << *runTPCE << " " << *evtTPCE;
		//	oFile << endl << "\t";
		//	oFile << "Efficiency: " << count_01/(max_wire_01 - min_wire_01 + 1);
		//	oFile << " over " << *lengthTPCE << " cm.";
		//	oFile << endl << "\t";
		//	oFile << *startxTPCE << " " << *startyTPCE << " " << *startzTPCE << " ";
		//	oFile << endl << "\t";
		//	oFile << *endxTPCE << " " << *endyTPCE << " " << *endzTPCE << " ";
		//	oFile << endl << endl;
		//}		
 
		cout << "EAST CRYO, PLANE COLL, LOGIC TPC 01 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE COLL, LOGIC TPC 01
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 0;
		_plane = 2;
		_min_wire = min_wire_01;
		_max_wire = max_wire_01;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_01;
		_tot_wires = wcount_01;
		_avg_pitch = holder_01;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_COLL_TPC01)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
	}
				
	// TPC23 COLL
	if(count_23 > wire_threshold) 
	{
		//average pitch
		holder_23 = holder_23/count_23;
		
		//find min wire
		for(int k = 0; k < 5600; ++k) 
		{
			if(wire_array_2_23[k] == 1) 
			{
				min_wire_23 = k;
				break;		
			}
		}
 
		//find max wire
		//for(int k = 5600; k > 0; --k) 
		for(int k = 5599; k > 0; --k) 
		{
			if(wire_array_2_23[k] == 1) 
			{
				max_wire_23 = k;
				break;		
			}
		}
 
		//find total amount of wires in the trace
        count_23 = 0; buco=0; wcount_23=0; max_buco=0;
 
		std::vector<std::vector<double>> holes_temp;
		std::vector<std::pair<int,int>> wire_hits;
		std::vector<TVector3> coordinates;
 
		for(int hit = 0; hit < wireTPCE_2.GetSize(); hit++)
		{
			if(tpcTPCE_2[hit] == 2 /*&& hasSPTPCE_2[hit]*/)
			{
				if(PRINT_EAST_COLL_TPC23)
				{
					cout << wireTPCE_2[hit] << " " << hit_time_TPCE_2[hit] << " " << hasSPTPCE_2[hit] << " " << (wireTPCE_2[hit] >= min_wire_23 && wireTPCE_2[hit] <= max_wire_23) << " " << w_2_23_e[wireTPCE_2[hit]] << " ";
					if(!isnan(hitx_TPCE_2[hit])){ cout << hitx_TPCE_2[hit] << " " << hity_TPCE_2[hit] << " " << hitz_TPCE_2[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_2[hit] >= min_wire_23 && wireTPCE_2[hit] <= max_wire_23)
				{
					wire_hits.push_back({wireTPCE_2[hit], hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_2[hit], hity_TPCE_2[hit], hitz_TPCE_2[hit]);
					coordinates.push_back(temp_coord);
				}
			}
			else if(tpcTPCE_2[hit] == 3)
			{
				if(PRINT_EAST_COLL_TPC23)
				{
					cout << wireTPCE_2[hit] + 2536 << " " << hit_time_TPCE_2[hit] << " " << hasSPTPCE_2[hit] << " " << (wireTPCE_2[hit] + 2536 >= min_wire_23 && wireTPCE_2[hit] + 2536 <= max_wire_23) << " " << w_2_23_e[wireTPCE_2[hit] + 2536] << " ";
					if(!isnan(hitx_TPCE_2[hit])){ cout << hitx_TPCE_2[hit] << " " << hity_TPCE_2[hit] << " " << hitz_TPCE_2[hit] << endl;}
					else {cout << -1 << " " << -1 << " " << -1 << endl;} 
				}
 
				if(wireTPCE_2[hit] + 2536 >= min_wire_23 && wireTPCE_2[hit] + 2536 <= max_wire_23)
				{
					wire_hits.push_back({wireTPCE_2[hit] + 2536, hit});
					TVector3 temp_coord;
					temp_coord.SetXYZ(hitx_TPCE_2[hit], hity_TPCE_2[hit], hitz_TPCE_2[hit]);
					coordinates.push_back(temp_coord);
				}
			}
		}
 
		double start_buco = -1;
		for(int k = min_wire_23; k <= max_wire_23; ++k) 
		{
            wcount_23 += w_2_23_e[k];
		  	if(w_2_23_e[k]>0 && wire_array_2_23[k] == 1 && buco>0) 
			{
				++count_23;
				
				std::vector<double> this_hole;
 
				int wire_before_hole = start_buco - 1;
				while(wire_array_2_23[wire_before_hole] == 0){wire_before_hole = wire_before_hole - 1;}
				
				std::pair<int,int> hh = wireToHit(wire_hits, coordinates, wire_before_hole, k);
 
				this_hole.push_back(buco);
				this_hole.push_back(hitx_TPCE_2[hh.first]);
				this_hole.push_back(hity_TPCE_2[hh.first]);
				this_hole.push_back(hitz_TPCE_2[hh.first]);
				this_hole.push_back(hitx_TPCE_2[hh.second]);
				this_hole.push_back(hity_TPCE_2[hh.second]);
				this_hole.push_back(hitz_TPCE_2[hh.second]);
				this_hole.push_back(wire_before_hole);
				this_hole.push_back(k);
				this_hole.push_back(hit_dirx_TPCE_2[hh.first]);
				this_hole.push_back(hit_diry_TPCE_2[hh.first]);
				this_hole.push_back(hit_dirz_TPCE_2[hh.first]);
				this_hole.push_back(hit_dirx_TPCE_2[hh.second]); 
				this_hole.push_back(hit_diry_TPCE_2[hh.second]);
				this_hole.push_back(hit_dirz_TPCE_2[hh.second]);
				this_hole.push_back(hit_time_TPCE_2[hh.first]);
				this_hole.push_back(hit_time_TPCE_2[hh.second]);
				this_hole.push_back(hit_mult_TPCE_2[hh.first]);
				this_hole.push_back(hit_mult_TPCE_2[hh.second]);
 
				holes_temp.push_back(this_hole);
				this_hole.clear();
 
				hcollbuchi->Fill(holder_23, buco);
				if(buco > max_buco) max_buco = buco;
				buco = 0;
				start_buco = -1;
			}
		  	else if(w_2_23_e[k]>0 && wire_array_2_23[k] == 1 && buco==0) ++count_23;
		  	else if(w_2_23_e[k]>0 && wire_array_2_23[k] == 0)
			{
				if(buco == 0) start_buco = k;
				++buco;
			}
		}
 
		wire_hits.clear();
 
		if(wcount_23 > (max_wire_23-min_wire_23+1) || count_23 > wcount_23) cout << "Coll 2 23 E " << count_23 << " " << wcount_23 << " " << min_wire_23 << " " << max_wire_23 << endl;
		if(max_buco<10) hist_eff_2.push_back(count_23/(wcount_23));
		if(max_buco<10) hist_average_pitch_2.push_back(holder_23);
		
		//print details of low efficiency traces
		//if(count_23/(max_wire_23 - min_wire_23 + 1) < upper_eff_print) 
		//{
		//	oFile << "EAST 23 COLL " << *runTPCE << " " << *evtTPCE;
		//	oFile << endl << "\t";
		//	oFile << "Efficiency: " << count_23/(max_wire_23 - min_wire_23 + 1);
		//	oFile << " over " << *lengthTPCE << " cm.";
		//	oFile << endl << "\t";
		//	oFile << *startxTPCE << " " << *startyTPCE << " " << *startzTPCE << " ";
		//	oFile << endl << "\t";
		//	oFile << *endxTPCE << " " << *endyTPCE << " " << *endzTPCE << " ";
		//	oFile << endl << endl;
		//}		
 
		cout << "EAST CRYO, PLANE COLL, LOGIC TPC 23 ";
		// WRITE VARIABLES IN TREE --> EAST CRYO, PLANE COLL, LOGIC TPC 23
		_run = *runTPCE;
		_evt = *evtTPCE; 
		_cryo = *cryostatTPCE;
		_tpc = 2;
		_plane = 2;
		_min_wire = min_wire_23;
		_max_wire = max_wire_23;
		_start.SetXYZ(*startxTPCE, *startyTPCE, *startzTPCE);
		_end.SetXYZ(*endxTPCE, *endyTPCE, *endzTPCE);
		_hit_wires = count_23;
		_tot_wires = wcount_23;
		_avg_pitch = holder_23;
		_trk_length = *lengthTPCE;
		_max_buco = max_buco;
		_which_t0 = *whicht0TPCE;
		_t0PFP = *t0_PFP_TPCE;
		_t0CRTTrack = *t0_CRT_Track_TPCE;
		_t0CRTHit = *t0_CRT_Hit_TPCE;
		_nholes = holes_temp.size();
		_wire_holes = holes_temp;
		_G4ID = *truthG4E;
		_truth_fraction = truth_fraction;

		if(PRINT_EAST_COLL_TPC23)
		{
			cout << _min_wire << " " << _max_wire << " | " << _nholes << " holes: ";
			for(const auto &hole : _wire_holes)
        	{
        		cout << "hole ";
        		for(const auto &hole_feature : hole)
        		{
            		cout << hole_feature << " ";
        		}
        	}
		}
 
		outree -> Fill();
		cout << "filled" << endl;
 
		holes_temp.clear();
	}
				
	//reset arrays and counters
	for(int j = 0; j < 5600; ++j) 
	{
		wire_array_2_01[j] = 0;
		wire_array_2_23[j] = 0;
	}
 
	count_01 = 0;
	count_23 = 0;	
 
	holder_01 = 0;	
	holder_23 = 0;
				
				
}
}		

	
	myFile -> Close();	
}
		
	




//		//
//		//	EFFICIENCY vs AVERAGE PITCH HISTOGRAMS AND TProfile
//		// 
	
	bin = 400;

	//initialize histograms
	TH2F *h_eff_vs_pitch_0 = new TH2F("h_eff_vs_pitch_0","",bin,0,4,bin,0,1.05);
	TH2F *h_eff_vs_pitch_1 = new TH2F("h_eff_vs_pitch_1","",bin,0,4,bin,0,1.05);
	TH2F *h_eff_vs_pitch_2 = new TH2F("h_eff_vs_pitch_2","",bin,0,4,bin,0,1.05);

	//fill histograms
	for(int i = 0; i < hist_eff_0.size(); ++i) h_eff_vs_pitch_0->Fill(hist_average_pitch_0[i], hist_eff_0[i]);
	for(int i = 0; i < hist_eff_1.size(); ++i) h_eff_vs_pitch_1->Fill(hist_average_pitch_1[i], hist_eff_1[i]);
	for(int i = 0; i < hist_eff_2.size(); ++i) h_eff_vs_pitch_2->Fill(hist_average_pitch_2[i], hist_eff_2[i]);


	/*
	auto *eff_vs_pitch_0_c = new TCanvas();
	h_eff_vs_pitch_0->SetTitle("Efficiency vs average pitch ind1, ALL;Average pitch;Efficiency");
	h_eff_vs_pitch_0->Draw("COLZ");
	eff_vs_pitch_0_c->SaveAs("h_eff_vs_pitch_all_0.root");

	auto *eff_vs_pitch_1_c = new TCanvas();
	h_eff_vs_pitch_1->SetTitle("Efficiency vs average pitch ind2, ALL;Average pitch;Efficiency");
	h_eff_vs_pitch_1->Draw("COLZ");
	eff_vs_pitch_1_c->SaveAs("h_eff_vs_pitch_all_1.root");

	auto *eff_vs_pitch_2_c = new TCanvas();
	h_eff_vs_pitch_2->SetTitle("Efficiency vs average pitch coll, ALL;Average pitch;Efficiency");
	h_eff_vs_pitch_2->Draw("COLZ");
	eff_vs_pitch_2_c->SaveAs("h_eff_vs_pitch_all_2.root"); 
	*/

	//create TProfile graphs
	auto *prof_eff_vs_pitch_0 = new TProfile("prof_eff_vs_pitch_0","",bin,0,4,0,1.05,"");
	auto *prof_eff_vs_pitch_1 = new TProfile("prof_eff_vs_pitch_1","",bin,0,4,0,1.05,"");
	auto *prof_eff_vs_pitch_2 = new TProfile("prof_eff_vs_pitch_2","",bin,0,4,0,1.05,"");

	for(int i = 0; i < hist_eff_0.size(); ++i) prof_eff_vs_pitch_0->Fill(hist_average_pitch_0[i], hist_eff_0[i], 1);
	for(int i = 0; i < hist_eff_1.size(); ++i) prof_eff_vs_pitch_1->Fill(hist_average_pitch_1[i], hist_eff_1[i], 1);
	for(int i = 0; i < hist_eff_2.size(); ++i) prof_eff_vs_pitch_2->Fill(hist_average_pitch_2[i], hist_eff_2[i], 1);

	/*
	auto *prof_eff_vs_pitch_0_c = new TCanvas();
	prof_eff_vs_pitch_0->SetTitle("TProfile Efficiency vs average pitch ind1, ALL;Average pitch;Efficiency");
	prof_eff_vs_pitch_0->Draw();
	prof_eff_vs_pitch_0_c->SaveAs("prof_eff_vs_pitch_0.root");

	auto *prof_eff_vs_pitch_1_c = new TCanvas();
	prof_eff_vs_pitch_1->SetTitle("TProfile Efficiency vs average pitch ind2, ALL;Average pitch;Efficiency");
	prof_eff_vs_pitch_1->Draw();
	prof_eff_vs_pitch_1_c->SaveAs("prof_eff_vs_pitch_1.root");

	auto *prof_eff_vs_pitch_2_c = new TCanvas();
	prof_eff_vs_pitch_2->SetTitle("TProfile Efficiency vs average pitch coll, ALL;Average pitch;Efficiency");
	prof_eff_vs_pitch_2->Draw();
	prof_eff_vs_pitch_2_c->SaveAs("prof_eff_vs_pitch_2.root");	
	*/

	TFile *outoutfilefile=new TFile(Form("all_histo_%s.root",sample_name.c_str()),"RECREATE");
	outoutfilefile -> cd();

	h_eff_vs_pitch_0->Write();
	h_eff_vs_pitch_1->Write();
	h_eff_vs_pitch_2->Write();

    prof_eff_vs_pitch_0->Write();        
    prof_eff_vs_pitch_1->Write(); 
	prof_eff_vs_pitch_2->Write();

	hind1buchi->Write();
	hind2buchi->Write();
	hcollbuchi->Write();


//
// OUTFILE CLOSING
//
    outoutfilefile->Close();


//
// TREE WRITING 
//	
	tree_outfile -> cd();
	outree->Write(0,TObject::kOverwrite);
	tree_outfile -> Close();

	/*
		TFile *outoutfilefile2=new TFile("wire_histo.root","RECREATE");	
		h0_0_w->Write();h0_1_w->Write();h0_2_w->Write();h0_3_w->Write();
        h0_0_e->Write();h0_1_e->Write();h0_2_e->Write();h0_3_e->Write();
        h1_01_w->Write();h2_01_w->Write();h1_23_w->Write();h2_23_w->Write();
        h1_01_e->Write();h2_01_e->Write();h1_23_e->Write();h2_23_e->Write();
        outoutfilefile2->Close();
	*/

}

