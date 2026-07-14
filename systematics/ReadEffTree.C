void ReadEffTree()
{
    //
    // READING THE TREE
    //

    TFile * infile = TFile::Open("tree_outfile_complete_tree_overlays.root","READ");

    TTree * intree = (TTree*)infile -> Get("outree");

    int run = -1;
	int evt = -1;
	short cryo = -1; //0 --> EAST		1 --> WEST
	short tpc = -1; //0 e 2 --> IND2/COLL, 0,1,2,3 --> IND1 
	short plane = -1;
	int min_wire = -1;
	int max_wire = -1;
	TVector3 *start = nullptr;
	TVector3 *end = nullptr;
	float hit_wires = -1;
	float tot_wires = -1;
	float avg_pitch = -1;
	float trk_length = -1;
    int max_buco = -1;
    short which_t0 = -1;
    float t0PFP = -1;
    float t0CRTTrack = -1;
    float t0CRTHit = -1;
	int nholes = -1;
	std::vector<std::vector<double>> *wire_holes = nullptr; // --> it contains:
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

    int g4id = -999;
    float truth_fraction = -1;

    intree -> SetBranchAddress("run",&run);
	intree -> SetBranchAddress("evt",&evt);
	intree -> SetBranchAddress("cryo",&cryo);
	intree -> SetBranchAddress("tpc",&tpc);
	intree -> SetBranchAddress("plane",&plane);
	intree -> SetBranchAddress("min_wire",&min_wire);
	intree -> SetBranchAddress("max_wire",&max_wire);
	intree -> SetBranchAddress("start",&start);
	intree -> SetBranchAddress("end",&end);
	intree -> SetBranchAddress("hit_wires",&hit_wires);
	intree -> SetBranchAddress("tot_wires",&tot_wires);
	intree -> SetBranchAddress("avg_pitch",&avg_pitch);
	intree -> SetBranchAddress("trk_length",&trk_length);
    intree -> SetBranchAddress("max_buco",&max_buco);
    intree -> SetBranchAddress("whicht0",&which_t0);
    intree -> SetBranchAddress("t0PFP",&t0PFP);
    intree -> SetBranchAddress("t0CRTTrack",&t0CRTTrack);
    intree -> SetBranchAddress("t0CRTHit",&t0CRTHit);
	intree -> SetBranchAddress("nholes",&nholes);
	intree -> SetBranchAddress("wire_holes",&wire_holes);
    intree -> SetBranchAddress("G4ID",&g4id);
    intree -> SetBranchAddress("truth_fraction",&truth_fraction); 

    //
    // DEFINING THE HISTOS
    //

    int bin = 400;
    auto *prof_eff_vs_pitch_0 = new TProfile("prof_eff_vs_pitch_0_fromTree","",bin,0,4,0,1.05,"");
	auto *prof_eff_vs_pitch_1 = new TProfile("prof_eff_vs_pitch_1_fromTree","",bin,0,4,0,1.05,"");
	auto *prof_eff_vs_pitch_2 = new TProfile("prof_eff_vs_pitch_2_fromTree","",bin,0,4,0,1.05,"");

    //
    // LOOP OVER TREE ENTRIES
    //

    for(int i = 0; i < intree->GetEntries(); i++)
    {
        intree->GetEntry(i);

        if(false) 
        {
            cout << run << " " << evt << " " << nholes << " holes " << endl;
            for(const auto &hole : *wire_holes)
            {
                cout << "hole " << endl;
                for(const auto &hole_feature : hole)
                {
                    cout << hole_feature << " ";
                }
                cout << endl;
            }
            cout << endl;
        }

        //if(tot_wires < 0){cout << "TOT_WIRES: " << tot_wires << " TPC: " << tpc << " PLANE: " << plane << " CRYO: " << cryo << " MIN_WIRE: " << min_wire << " MAX_WIRE: " << max_wire << endl;}

        if(g4id == -1 || truth_fraction < 0.8)continue;

        if(plane == 0)
        {
            if(max_buco<10)prof_eff_vs_pitch_0 -> Fill(avg_pitch,hit_wires/tot_wires,1);
        }
        else if(plane == 1)
        {
            if(max_buco<10)prof_eff_vs_pitch_1 -> Fill(avg_pitch,hit_wires/tot_wires,1);
        }
        else if(plane == 2)
        {
            if(max_buco<10)prof_eff_vs_pitch_2 -> Fill(avg_pitch,hit_wires/tot_wires,1);
        }
    }

    TFile * histo_file = TFile::Open("all_histo.root","RECREATE");
    histo_file -> cd();
    prof_eff_vs_pitch_0->Write(0,TObject::kOverwrite);
    prof_eff_vs_pitch_1->Write(0,TObject::kOverwrite);
    prof_eff_vs_pitch_2->Write(0,TObject::kOverwrite);

    histo_file -> Close();
    infile -> Close();

}