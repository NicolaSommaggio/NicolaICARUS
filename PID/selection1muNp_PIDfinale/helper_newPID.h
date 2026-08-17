
struct Node { int left, right, feature; double threshold, value; };

struct singleTree {
    std::vector<Node> nodes;
    double predict(const std::vector<double>& x) const {
        int cur = 0;
        while (nodes[cur].feature != -2)
            cur = x[nodes[cur].feature] <= nodes[cur].threshold
                  ? nodes[cur].left : nodes[cur].right;
        return nodes[cur].value;
    }
};

struct GBDTModel {
    double learning_rate;
    int n_classes;
    std::vector<double> init_scores;  // uno per classe
    std::vector<singleTree> trees;          // n_stages * n_classes alberi

    // Ritorna il vettore di probabilità per ogni classe (softmax)
    std::vector<double> predict_proba(const std::vector<double>& x) const {
        int n_stages = trees.size() / n_classes;

        // Accumula score per ogni classe
        std::vector<double> scores = init_scores;
        for (int s = 0; s < n_stages; s++)
            for (int k = 0; k < n_classes; k++)
                scores[k] += learning_rate * trees[s * n_classes + k].predict(x);

        // Softmax
        double max_s = *std::max_element(scores.begin(), scores.end());
        std::vector<double> proba(n_classes);
        double sum = 0.0;
        for (int k = 0; k < n_classes; k++) {
            proba[k] = std::exp(scores[k] - max_s);
            sum += proba[k];
        }
        for (auto& p : proba) p /= sum;
        return proba;
    }

    int predict_class(const std::vector<double>& x) const {
        auto proba = predict_proba(x);
        return std::max_element(proba.begin(), proba.end()) - proba.begin();
    }
};

GBDTModel load_model(const char* filename) {
    std::ifstream f(filename);
    GBDTModel model;

    int n_stages;
    f >> model.learning_rate >> n_stages >> model.n_classes;

    model.init_scores.resize(model.n_classes);
    for (int k = 0; k < model.n_classes; k++)
        f >> model.init_scores[k];

    for (int s = 0; s < n_stages * model.n_classes; s++) {
        int n_nodes; f >> n_nodes;
        singleTree tree;
        for (int i = 0; i < n_nodes; i++) {
            Node n;
            f >> n.left >> n.right >> n.feature >> n.threshold >> n.value;
            tree.nodes.push_back(n);
        }
        model.trees.push_back(tree);
    }
    return model;
}


int find_idx(double rr)
{
    if(rr==25.)return 47;
    double a = std::round(rr);
    double max=-1;
    double min=-1;
    if(a>rr)
    {
        min=a-0.5;
        max=a;
    }
    else
    {
        min=a;
        max=a+0.5;
    }
    
    int n = (max-1.5)/0.5;

    return n;
}

int true_selection(
  const caf::SRSpillProxy* sr, 
  const caf::Proxy<caf::SRSlice>& islc, 
  std::size_t ipfp, 
  int plane = 999, 
  double dedx_min = 1., 
  double dedx_max = 30., 
  double rr_min = 1.,
  double rr_max = 25., 
  double mult = 1.,
  double wm_cut = true
)
{
    //-1 : unclassified
    //0 : muon rising
    //1 : muon mip
    //2 : proton rising
    //3 : proton interacting
    //4 : pion rising
    //5 : pion interacting 

    int bestplane = 2;
    if(plane == 0 || plane == 1 || plane == 2)
    {
      bestplane = plane; //--> IF YOU DON'T WANT A SPECIFIC PLANE DO NOT SPECIFY PLANE OR SET PLANE TO 999
    }
    
    bool hasValidHits=false;

    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()>0)
    {
      //check if it has valid hits for for computing the likelihood
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
      {
        if(
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<rr_max && 
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>rr_min && 
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>dedx_min && 
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<dedx_max && 
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].mult <= mult)
        {
          if(wm_cut)
          {
            if(wiremod::WireModHitCut(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].phi, islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch, islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].integral, plane)){continue;}
          }
          hasValidHits=true;
        }
      }
    }

    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()==0 || hasValidHits==false)
    {
      if(plane == 0 || plane == 1 || plane == 2)
      {
        return -1; //--> IF YOU DON'T WANT A SPECIFIC PLANE DO NOT SPECIFY PLANE OR SET PLANE TO 999
      }

      bestplane = islc.reco.pfp[ipfp].trk.bestplane;
      if(bestplane==-1)return -1;

      //check if it has valid hits for for computing the likelihood
      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
      {
        if(
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<rr_max && 
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>rr_min && 
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>dedx_min && 
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<dedx_max &&
          islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].mult <= mult)
        {
          if(wm_cut)
          {
            if(wiremod::WireModHitCut(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].phi, islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].pitch, islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].integral, plane)){continue;}
          }
          hasValidHits=true;
        }
      }

      if(hasValidHits==false)return -1;
    }

    //controllo che vertice reco e true siano abbastanza vicini per eliminare la componente di cosmici
    TVector3 vertex_reco;
    vertex_reco.SetXYZ(islc.vertex.x, islc.vertex.y, islc.vertex.z);
    TVector3 vertex_true;
    vertex_true.SetXYZ(islc.truth.position.x, islc.truth.position.y, islc.truth.position.z);
    if((vertex_true-vertex_reco).Mag()>100.)return -1;

    //compute distance between the endpoints
    TVector3 end_true;
    end_true.SetXYZ(islc.reco.pfp[ipfp].trk.truth.p.end.x, islc.reco.pfp[ipfp].trk.truth.p.end.y, islc.reco.pfp[ipfp].trk.truth.p.end.z);
    TVector3 end_hit;
    int nhits = (int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size();
    double endx = islc.reco.pfp[ipfp].trk.calo[bestplane].points[nhits-1].x;
    double endy = islc.reco.pfp[ipfp].trk.calo[bestplane].points[nhits-1].y;
    double endz = islc.reco.pfp[ipfp].trk.calo[bestplane].points[nhits-1].z;
    end_hit.SetXYZ(endx,endy,endz);
    double end_distance = (end_hit-end_true).Mag();

    //looking at the energy match
    double daughter_electron_energy_match = -1;
    double daughter_proton_energy_match = -1;
    double daughter_pion_energy_match = -1;

    for (const auto& true_p : sr->true_particles)
    {
      for (auto const& match: islc.reco.pfp[ipfp].trk.truth.matches)
      {
        if(true_p.G4ID==match.G4ID)
        {
          if(std::abs(true_p.pdg)==11 && (int)true_p.parent == islc.reco.pfp[ipfp].trk.truth.p.G4ID)
          {
            if(match.energy/3. > daughter_electron_energy_match)daughter_electron_energy_match=match.energy/3.;
          }
          if(std::abs(true_p.pdg)==2212 && (int)true_p.parent == islc.reco.pfp[ipfp].trk.truth.p.G4ID)
          {
            if(match.energy/3. > daughter_proton_energy_match)daughter_proton_energy_match=match.energy/3.;
          }
          if(std::abs(true_p.pdg)==211 && (int)true_p.parent == islc.reco.pfp[ipfp].trk.truth.p.G4ID)
          {
            if(match.energy/3. > daughter_pion_energy_match)daughter_pion_energy_match=match.energy/3.;
          }
        }
      }
    } 

    int true_class=-1;
    //TRUE CLASSIFICATION
    if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==13)
    {
        if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness>=0.5 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity>=0.5)
        {
            if(daughter_electron_energy_match<=0.01)
            {
                if(end_distance<=3.)true_class = 0;       
                else true_class = 1;
            }
            else true_class=6;
        } 
        else true_class=-1;
    }
    if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==2212)
    {
        if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness>=0.3 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity>=0.3)
        {
            if
            ( 
                (islc.reco.pfp[ipfp].trk.truth.p.end_process==54 && end_distance <= 5.) ||
                (islc.reco.pfp[ipfp].trk.truth.p.end_process!=54 && (daughter_proton_energy_match>0.055 || daughter_pion_energy_match>0.055))
            ) true_class = 2;
            else true_class = 3;
        }
        else true_class=-1;
    }
    if(std::abs(islc.reco.pfp[ipfp].trk.truth.p.pdg)==211)
    {
        if(islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_completeness>=0.3 && islc.reco.pfp[ipfp].trk.truth.bestmatch.hit_purity>=0.3)
        {
            if(( islc.reco.pfp[ipfp].trk.truth.p.end_process==3 || islc.reco.pfp[ipfp].trk.truth.p.end_process==45 ) && end_distance <= 1.5) true_class = 4;
            else true_class = 5;
        }
        else true_class=-1;
    }

    return true_class;
}

int find_best_plane(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
    //there should be at least 1 hit with rr between 1 cm and 25 cm and with the dedx<30 MeV/cm
    //if the condition is satified in collection, collection is the best plane regardless
    int bestplane = 2;
    bool hasValidHits=false;

    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()>0)
    {
        //check if it has valid hits for for computing the likelihood
        for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
        {
            if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
            {
                hasValidHits=true;
            }
        }   
    }
    
    if((int)islc.reco.pfp[ipfp].trk.calo[bestplane].points.size()==0 || hasValidHits==false)
    {
        bestplane = islc.reco.pfp[ipfp].trk.bestplane;
        if(bestplane==-1)return -1;

        //check if it has valid hits for for computing the likelihood
        for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[bestplane].points.size(); ++ihit)
        {
            if(islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[bestplane].points[ihit].dedx<30.)
            {
              hasValidHits=true;
            }
        }

        if(hasValidHits==false)return -1;
    }

    return bestplane;
}


std::array<std::vector<TH1D*>,6> load_prob_densities(std::string plane, TFile *inputf)
{
  std::vector<TH1D*> probDensitiesMU_class0;
    std::vector<TH1D*> probDensitiesMU_class1;
    std::vector<TH1D*> probDensitiesPRO_class2;
    std::vector<TH1D*> probDensitiesPRO_class3;
    std::vector<TH1D*> probDensitiesPI_class4;
    std::vector<TH1D*> probDensitiesPI_class5;

    std::array<std::vector<TH1D*>,6> probDensities={
        probDensitiesMU_class0,
        probDensitiesMU_class1,
        probDensitiesPRO_class2,
        probDensitiesPRO_class3,
        probDensitiesPI_class4,
        probDensitiesPI_class5
    };

    std::vector<std::pair<std::string,int>> classes = {{"muon_class0",0}, {"muon_class1",1}, {"proton_class2",2}, {"proton_class3",3}, {"pion_class4",4}, {"pion_class5",5}};

    for(const auto &clas : classes)
    {
        for(double i=1.5; i<=25.0; i+=0.5)
        {
            TDirectory *dclass = (TDirectory*)inputf->Get(clas.first.c_str());
            TDirectory *d = (TDirectory*)dclass->Get(Form("%.1f",i));

            TH1D *dEdx_coll = (TH1D*)d->Get(Form("dEdx_%s_rr_%.1f", plane.c_str(), i));
            probDensities[clas.second].push_back(dEdx_coll);
        }
    }

    return probDensities;
}

std::vector<double> likelihood(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp, std::array<std::vector<TH1D*>,6> &coll_pd, std::array<std::vector<TH1D*>,6> &ind1_pd, std::array<std::vector<TH1D*>,6> &ind2_pd)
{    
  auto probDensities = coll_pd;

  std::vector<double> likelihood_return;
  for(int i=0; i<15; i++){likelihood_return.push_back(-1);}
  
  std::vector<double> likelihood_ratios;
            
  std::array<double,6> lkl; //it contains the likelihoods in the 6 hypotheses for the current track

  bool hasValidHits=false;

  //ciclo sulla particle hypothesis
  for(int j=0; j<6; j++)
  {
    auto density = probDensities[j];
    double log_lkh=0;
        
    for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[2].points.size(); ++ihit)
    {
      if(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx<30.)
      {
        hasValidHits=true;
        int idx=find_idx(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].rr);
        int bin = density[idx]->FindBin(islc.reco.pfp[ipfp].trk.calo[2].points[ihit].dedx);
        log_lkh = log_lkh + (-1*std::log(density[idx]->GetBinContent(bin)/density[idx]->GetBinWidth(bin)));
        //if(std::isinf(log_lkh)){ cout << dat.track.rr->at(hit) << " " << dat.track.dE->at(hit) << " | " << density[idx]->GetBinContent(bin) << " " << std::log(density[idx]->GetBinContent(bin)) << " " << density[idx]->GetBinWidth(bin) << " | " << par << " " << subclass << " " << idx << "(" << dat.track.rr->at(hit) << ")" << " " << bin << endl; break;}
      }
    }
    lkl[j]=log_lkh;
  }//cycle over particle hypothesis

  if(!hasValidHits)
  {
    int use_plane = -1;
    if(islc.reco.pfp[ipfp].trk.bestplane==0 || islc.reco.pfp[ipfp].trk.bestplane==1 || islc.reco.pfp[ipfp].trk.bestplane==2){use_plane=islc.reco.pfp[ipfp].trk.bestplane;}
    else{use_plane=2;}

    if(use_plane == 0){probDensities = ind1_pd;}
    else if(use_plane == 1){probDensities = ind2_pd;}
    else if(use_plane == 2){probDensities = coll_pd;}

    for(int j=0; j<6; j++)
    { 
      auto density = probDensities[j];
      double log_lkh=0;

      for(std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit)
      {
        if(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr<25 && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr>1. && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx>1. && islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx<30.)
        {
          hasValidHits=true;
          int idx=find_idx(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
          int bin = density[idx]->FindBin(islc.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
          log_lkh = log_lkh + (-1*std::log(density[idx]->GetBinContent(bin)/density[idx]->GetBinWidth(bin)));
          //if(std::isinf(log_lkh)){ cout << dat.track.rr->at(hit) << " " << dat.track.dE->at(hit) << " | " << density[idx]->GetBinContent(bin) << " " << std::log(density[idx]->GetBinContent(bin)) << " " << density[idx]->GetBinWidth(bin) << " | " << par << " " << subclass << " " << idx << "(" << dat.track.rr->at(hit) << ")" << " " << bin << endl; break;}
        }
      }
      lkl[j]=log_lkh;
    }//cycle over particle hypothesis
  }

  if(!hasValidHits){return likelihood_return;}

  //GETTING LIKELIHOOD RATIOS FOR EACH TRACK
  for(int k=0; k<6; k++)
  {
    for(int t=k+1; t<6; t++)
    {
      //if(isnan(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90))cout << "nan value ecountered " << lkl[k] << " " << lkl[t] << endl;
      likelihood_ratios.push_back(std::atan((lkl[k]-lkl[t])/3.)*180/M_PI/90);
    }
  }  
  return likelihood_ratios;
}

double compute_depE(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp, int plane)
{
    double dep_E = 0;
    for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp].trk.calo[plane].points.size(); ++ihit )
    {
        if(islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].rr <= 5.)
        {   
            dep_E = dep_E + islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].dedx*islc.reco.pfp[ipfp].trk.calo[plane].points[ihit].pitch;
        }
    }
    return dep_E;
}

std::vector<double> compute_ke(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp, int plane)
{
  double ke_calo = islc.reco.pfp[ipfp].trk.calo[plane].ke;

  TVector3 p_from_range_mu;
  p_from_range_mu.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_muon)*islc.reco.pfp[ipfp].trk.dir.z);
  double ke_range_mu = sqrt(pow(105.658,2)+pow(p_from_range_mu.Mag()*1000,2))-105.658;

  TVector3 p_from_range_pi;
  p_from_range_pi.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.x, (islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.y, (islc.reco.pfp[ipfp].trk.rangeP.p_pion)*islc.reco.pfp[ipfp].trk.dir.z);
  double ke_range_pi = sqrt(pow(139.570,2)+pow(p_from_range_pi.Mag()*1000,2))-139.570;

  TVector3 p_from_range_pro;
  p_from_range_pro.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
  double ke_range_pro = sqrt(pow(938.3,2)+pow(p_from_range_pro.Mag()*1000,2))-938.3;

  double ke_difference_mu_hyp = (ke_range_mu - ke_calo)/ke_calo;
  double ke_difference_pi_hyp = (ke_range_pi - ke_calo)/ke_calo;
  double ke_difference_pro_hyp = (ke_range_pro - ke_calo)/ke_calo;

  //return {ke_difference_mu_hyp,ke_difference_pi_hyp,ke_difference_pro_hyp};
  return {ke_calo, ke_range_mu, ke_range_pro, ke_range_pi};
}

std::vector<double> compute_daughter_vars(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
    int daughter_max_hits = 0;
    int daughter_id = -1;
    std::vector<double> daughter_rr;
    std::vector<double> daughter_dedx;
    double daughter_depE = 0;
    std::array<double,3> daughter_direction={-1,-1,-1};
    double angle_end = -1;

    for(int daughter=0; daughter<int(islc.reco.pfp[ipfp].ndaughters); daughter++)
    {  
        int d = islc.reco.pfp[ipfp].daughters[daughter];
    
        std::vector<double> temp_dedx;
        std::vector<double> temp_rr;
        double temp_depE = 0;
        std::array<double,3> temp_direction = {-1,-1,-1};
        for (std::size_t jpfp(0); jpfp < islc.reco.npfp; ++jpfp)
        {
            if(islc.reco.pfp[jpfp].id!=d)continue;

            int d_bestplane = islc.reco.pfp[jpfp].trk.bestplane;

            if(d_bestplane==-1)continue;

            for ( std::size_t ihit(0); ihit < islc.reco.pfp[jpfp].trk.calo[d_bestplane].points.size(); ++ihit )
            {
              if(islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].rr <= 5)
              {
                temp_rr.push_back(islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].rr);
                temp_dedx.push_back(islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].dedx);
                temp_depE = temp_depE + islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].dedx*islc.reco.pfp[jpfp].trk.calo[d_bestplane].points[ihit].pitch;
              }
            }
            temp_direction = {islc.reco.pfp[jpfp].trk.dir.x, islc.reco.pfp[jpfp].trk.dir.y, islc.reco.pfp[jpfp].trk.dir.z};

            break;
        }

        if(temp_rr.empty())continue;
        if((int)temp_rr.size()>daughter_max_hits)
        {
          daughter_max_hits=(int)temp_rr.size();
          daughter_id = daughter;
          daughter_rr = temp_rr;
          daughter_dedx = temp_dedx;
          daughter_depE = temp_depE;
          daughter_direction = temp_direction;
        }
    }

    std::array<double,3> track_direction_end = {islc.reco.pfp[ipfp].trk.dir_end.x, islc.reco.pfp[ipfp].trk.dir_end.y, islc.reco.pfp[ipfp].trk.dir_end.z};
    angle_end = std::acos(track_direction_end[0]*daughter_direction[0] + track_direction_end[1]*daughter_direction[1] + track_direction_end[2]*daughter_direction[2]);

    if(daughter_id==-1)return {-1,-1};

    return {daughter_depE,angle_end};
}


/*
//PID model OLD MC ------------------------------------------------------------------

TFile * f_prob_densities_coll = TFile::Open("/exp/icarus/app/users/nsommagg/NicolaICARUS/PID/selection_newPID/collection.root", "READ");
TFile * f_prob_densities_ind1 = TFile::Open("/exp/icarus/app/users/nsommagg/NicolaICARUS/PID/selection_newPID/induction1.root", "READ");
TFile * f_prob_densities_ind2 = TFile::Open("/exp/icarus/app/users/nsommagg/NicolaICARUS/PID/selection_newPID/induction2.root", "READ");

std::string path_BDT_model =  "/exp/icarus/app/users/nsommagg/NicolaICARUS/PID/selection_newPID/provaBDT_dedxmag1/GBDT_dedx_mag1.txt";
GBDTModel model = load_model(path_BDT_model.c_str());
//-----------------------------------------------------------------------------------
*/



//PID model NEW MC ------------------------------------------------------------------

//FNAL
TFile * f_prob_densities_coll = TFile::Open("/exp/icarus/data/users/nsommagg/PDF_COLL_NEW_MC.root", "READ");
TFile * f_prob_densities_ind1 = TFile::Open("/exp/icarus/data/users/nsommagg/PDF_IND1_NEW_MC.root", "READ");
TFile * f_prob_densities_ind2 = TFile::Open("/exp/icarus/data/users/nsommagg/PDF_IND2_NEW_MC.root", "READ");

//CNAF
//TFile * f_prob_densities_coll = TFile::Open("PDF_COLL_NEW_MC.root", "READ");
//TFile * f_prob_densities_ind1 = TFile::Open("PDF_IND1_NEW_MC.root", "READ");
//TFile * f_prob_densities_ind2 = TFile::Open("PDF_IND2_NEW_MC.root", "READ");

//FNAL
std::string path_BDT_model =  "/exp/icarus/app/users/nsommagg/NicolaICARUS/PID/selection_newPID/TEST_BDT_NEW_MC/GBDT_MODEL_NEW_MC_EXPORT.txt";

//CNAF
//std::string path_BDT_model =  "GBDT_MODEL_NEW_MC_EXPORT.txt";

GBDTModel model = load_model(path_BDT_model.c_str());
//-----------------------------------------------------------------------------------


std::array<std::vector<TH1D*>,6> prob_d_coll = load_prob_densities("coll",f_prob_densities_coll);
std::array<std::vector<TH1D*>,6> prob_d_ind1 = load_prob_densities("ind1",f_prob_densities_ind1);
std::array<std::vector<TH1D*>,6> prob_d_ind2 = load_prob_densities("ind2",f_prob_densities_ind2);

int PIDclass(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
  //-1 : unclassified
  //0 : muon rising
  //1 : muon mip
  //2 : proton rising
  //3 : proton interacting
  //4 : pion rising
  //5 : pion interacting 

  int bestplane = find_best_plane(islc,ipfp);

  if(bestplane == -1)return -1; //no valid hits (RR < 25 cm and RR > 1 cm and dEdx < 30 MeV/cm and dEdx > 1 MeV/cm) in neither of the 3 planes

  std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
  double depE = compute_depE(islc,ipfp,bestplane);
  std::vector<double> dvars = compute_daughter_vars(islc,ipfp);
  //std::vector<double> ke_vars = compute_ke(islc,ipfp,bestplane);

  std::vector<double> track_features;

  track_features.insert(track_features.end(), lr.begin(), lr.end());
  track_features.push_back(depE);
  track_features.insert(track_features.end(), dvars.begin(), dvars.end());
  //track_features.insert(track_features.end(), ke_vars.begin(), ke_vars.end());

  int prediction = model.predict_class(track_features); 
  
  return prediction;
}

std::vector<double> PIDproba(const caf::Proxy<caf::SRSlice>& islc, std::size_t ipfp)
{
  //it return a 6D vector containing the probabilities for the track to be each of the 6 classes 

  int bestplane = find_best_plane(islc,ipfp);

  if(bestplane == -1)return {-1.,-1.,-1.,-1.,-1.,-1.}; //no valid hits (RR < 25 cm and RR > 1 cm and dEdx < 30 MeV/cm) in neither of the 3 planes

  std::vector<double> lr = likelihood(islc,ipfp, prob_d_coll, prob_d_ind1, prob_d_ind2);
  double depE = compute_depE(islc,ipfp,bestplane);
  std::vector<double> dvars = compute_daughter_vars(islc,ipfp);
  //std::vector<double> ke_vars = compute_ke(islc,ipfp,bestplane);

  std::vector<double> track_features;

  track_features.insert(track_features.end(), lr.begin(), lr.end());
  track_features.push_back(depE);
  track_features.insert(track_features.end(), dvars.begin(), dvars.end());
  //track_features.insert(track_features.end(), ke_vars.begin(), ke_vars.end());

  std::vector<double> prediction_proba = model.predict_proba(track_features); 
  
  return prediction_proba;
}

