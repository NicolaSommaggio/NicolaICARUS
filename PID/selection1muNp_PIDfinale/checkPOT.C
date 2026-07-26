//find /pnfs/sbn/scratch/users/twester/Run2_WM_ReCAF2026/flatcaf/ -mindepth 3 -maxdepth 3 -type f -name "*.root" > file_list.txt
void sum_pot() {

    ifstream file_list("file_list.txt");

    std::string line;

    double tot_pot = 0;

    TFile * f;

    while(std::getline(file_list,line))
    {
        std::stringstream ss(line);
        std::string path;
        ss >> path;

        cout << path << endl;

        f = TFile::Open(path.c_str(),"READ");

        f->Print();

        TH1D * poth = (TH1D*)f->Get("TotalPOT");

        tot_pot = tot_pot + poth->GetBinContent(1);

        cout << tot_pot << endl;

        f->Close();
        //delete f;
        //delete poth;

    }


}