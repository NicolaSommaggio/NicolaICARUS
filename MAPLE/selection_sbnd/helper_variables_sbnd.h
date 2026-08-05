//#include "helper_sbnd.h"
#include "helper_sbnd_original.h"



double Neutrino_energy_visible_true_Np_SLICE(const caf::SRSpillProxy* sr,
                                      const caf::Proxy<caf::SRSlice>& islc)
{

    // ----------------- Counters -----------------

    double E_p_sum = 0.0;
    double E_mu = 0.0;

    const int use_plane = DEFAULT_PLANE;

    // ----------------- Primary loop -----------------
    for (const auto& ipart : islc.truth.prim) {

        if (ipart.G4ID < 0 || ipart.cryostat < 0) continue;

        TruthPID pid = classify_truth_pid(ipart.pdg); 
        double depE = ipart.plane[ipart.cryostat][use_plane].visE * 1000.0;

        // -------- Muons --------
        if (pid == TruthPID::Muon) {
            TVector3 p_mu(ipart.genp.x*1000.0,ipart.genp.y*1000.0,ipart.genp.z*1000.0); // MeV

            E_mu = std::sqrt(p_mu.Mag()*p_mu.Mag() + MUON_MASS*MUON_MASS); 
        }


        // -------- Protons --------
        if (pid == TruthPID::Proton) {
            for (const auto& d : sr->true_particles)
                if (same_g4id(d.parent, ipart.G4ID))
                    depE += d.plane[ipart.cryostat][use_plane].visE * 1000.0;


            if (depE > PROTON_KE_MIN) {

            TVector3 p_pr(ipart.genp.x,ipart.genp.y,ipart.genp.z); // GeV

            double T_p = kinetic_energy(PROTON_MASS, p_pr.Mag());  // momentum in GeV!! 

        E_p_sum += (T_p + PROTON_BINDING_ENERGY);                
            }
        }
    }

        //std::cout << "Energy " << (E_p_sum+E_mu)/1000.0 << endl;
    return (E_p_sum+E_mu)/1000.0;
}

struct LongestProtonResult {
    int index = -1;
    double length = 0.0;
};


LongestProtonResult find_longest_proton(const caf::Proxy<caf::SRSlice>& islc,
                                        int muon_index,
                                        double dist_cut)
{
    LongestProtonResult result;

    for (std::size_t ipfp = 0; ipfp < islc.reco.npfp; ++ipfp) {

        if (static_cast<int>(ipfp) == muon_index) {
            continue;
        }

        PFPID pid = static_cast<PFPID>(id_pfp(islc, ipfp, dist_cut));

        if (pid != PFPID::Proton) {
            continue;
        }

        double length = islc.reco.pfp[ipfp].trk.len;  // adjust if needed

        if (length > result.length) {
            result.length = length;
            result.index  = static_cast<int>(ipfp);
        }
    }

    return result;
}


double Transverse_angle ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    double p_p_tot = sqrt(p_p_x*p_p_x+p_p_y*p_p_y+p_p_z*p_p_z);

    //Total momentum
    p_tot_x=p_p_x+p_mu_x;
    p_tot_y=p_p_y+p_mu_y;
    p_tot_z=p_p_z+p_mu_z;

    //TRANSVERSE PLANE - ANGLE! 
    double norm_mu=sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y);
    double norm_pro=sqrt(p_p_x*p_p_x+p_p_y*p_p_y);    
       
     
    return (p_mu_x*p_p_x+p_mu_y*p_p_y)/(norm_mu*norm_pro); 
    }

double T3D_angle_mup ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    TVector3 mu_vector(p_mu_x,p_mu_y,p_mu_z);

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    TVector3 pro_vector(p_p_x,p_p_y,p_p_z);

     
    return cos(mu_vector.Angle(pro_vector)); 
    }


double Transverse_mom_reco ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) {
    
    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;
    float E_mu,E_p;    
    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV
    E_mu=sqrt(p_mu_tot*p_mu_tot+(105.658*105.658)/(1000*1000));

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.z;
    double p_p_tot = sqrt(p_p_x*p_p_x+p_p_y*p_p_y+p_p_z*p_p_z);
    E_p=sqrt(p_p_tot*p_p_tot+(938.272*938.272)/(1000*1000));

    //Total momentum
    p_tot_x=p_p_x+p_mu_x;
    p_tot_y=p_p_y+p_mu_y;
    p_tot_z=p_p_z+p_mu_z;

    double p_T=sqrt(p_tot_x*p_tot_x+p_tot_y*p_tot_y);


return p_T;
}
/*
double Transverse_mom_reco_Np ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu) {

        
        float p_p_x=0; float p_p_y =0; float p_p_z=0;
        float p_tot_x =0; float p_tot_y=0; float p_tot_z=0;
        int ipfp_pro = -1; 

        float p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
        float p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
        float p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;

        for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
            if(int(ipfp)==ipfp_mu)continue;
            if(id_pfp(islc, ipfp,10)==1){
                p_p_x +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x;
                p_p_y +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y;
                p_p_z +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z;
                ipfp_pro=ipfp;

            }
                    }

    if(ipfp_mu!=-1 && ipfp_pro!=-1){
            p_tot_x=p_p_x+p_mu_x;
            p_tot_y=p_p_y+p_mu_y;
            p_tot_z=p_p_z+p_mu_z;
                      
    }
        
       double p_T = (sqrt(p_tot_x*p_tot_x+p_tot_y*p_tot_y));


return p_T;
}

*/

double Transverse_mom_reco_mu ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) {
    
    float p_mu_x,p_mu_y;   
    //Muon momentum
    p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;


    double p_T=sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y);


return p_T;
}

double Transverse_mom_reco_pro ( const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, int ipfp_pro ) {
    
    float p_p_x,p_p_y;

    //Proton momentum
    p_p_x =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc.reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc.reco.pfp[ipfp_pro].trk.dir.y;


    double p_T=sqrt(p_p_x*p_p_x+p_p_y*p_p_y);


return p_T;
}



const SpillMultiVar kProton_length_leading([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;

    if (longest_proton_idx >= 0) {
        vector_active.push_back(islc.reco.pfp[longest_proton_idx].trk.len);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kProton_kinetic_leading([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
        if (longest_proton_idx >= 0) {
        
        TVector3 p_pr(islc.reco.pfp[longest_proton_idx].trk.dir.x * islc.reco.pfp[longest_proton_idx].trk.rangeP.p_proton,
              islc.reco.pfp[longest_proton_idx].trk.dir.y * islc.reco.pfp[longest_proton_idx].trk.rangeP.p_proton,
              islc.reco.pfp[longest_proton_idx].trk.dir.z * islc.reco.pfp[longest_proton_idx].trk.rangeP.p_proton);

        double ke = kinetic_energy(PROTON_MASS, p_pr.Mag());
        vector_active.push_back(ke/1000.0);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});



const SpillMultiVar kProton_chi2mu([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
        if (longest_proton_idx >= 0) {
    std::vector<double> chi2 = compute_chi2(islc, longest_proton_idx, DEFAULT_PLANE);
        vector_active.push_back(chi2[static_cast<int>(Chi2PID::Muon)]);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kProton_chi2pro([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
    std::vector<double> chi2 = compute_chi2(islc, longest_proton_idx, DEFAULT_PLANE);
        vector_active.push_back(chi2[static_cast<int>(Chi2PID::Proton)]);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kMuon_chi2mu([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
    std::vector<double> chi2 = compute_chi2(islc, ipfp_mu, DEFAULT_PLANE);
        vector_active.push_back(chi2[static_cast<int>(Chi2PID::Muon)]);
         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kMuon_chi2pro([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
    std::vector<double> chi2 = compute_chi2(islc, ipfp_mu, DEFAULT_PLANE);
        vector_active.push_back(chi2[static_cast<int>(Chi2PID::Proton)]);
         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kVertex_x([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

        vector_active.push_back(islc.vertex.x); }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kVertex_y([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

        vector_active.push_back(islc.vertex.y); }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kVertex_z([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

        vector_active.push_back(islc.vertex.z); }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});



const SpillMultiVar kMuon_length([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.len);
         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kMuon_endx([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.end.x);
         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kMuon_endy([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.end.y);
         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kMuon_endz([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.end.z);
         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kProton_endx([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(islc.reco.pfp[longest_proton_idx].trk.end.x);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kProton_endy([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(islc.reco.pfp[longest_proton_idx].trk.end.y);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kProton_endz([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(islc.reco.pfp[longest_proton_idx].trk.end.z);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kTransverse_angle([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(Transverse_angle(islc,ipfp_mu,longest_proton_idx));
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kT3D_angle_mup([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(T3D_angle_mup(islc,ipfp_mu,longest_proton_idx));
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kTransverse_mom_reco([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(Transverse_mom_reco(islc,ipfp_mu,longest_proton_idx));
        //  NEED TO USE THIS ONE Transverse_mom_reco_Np!!!!!!
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kTransverse_mom_reco_mu([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(Transverse_mom_reco_mu(islc,ipfp_mu,longest_proton_idx));
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kTransverse_mom_reco_pro([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(Transverse_mom_reco_pro(islc,ipfp_mu,longest_proton_idx));
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kNumber_protons([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
            ParticleCounts counts = count_particles(islc, ipfp_mu, 10);
        vector_active.push_back(counts.protons);
        }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});




const SpillMultiVar kTrue_visible_Enu_reco([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        TruthClass cls = classification_type_debug(sr, islc);
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

        switch (cls) {
            case TruthClass::kOneMuOneP: vector_active.push_back(Neutrino_energy_visible_true_Np_SLICE(sr,islc)); break;
            case TruthClass::kOneMuNp:   vector_active.push_back(Neutrino_energy_visible_true_Np_SLICE(sr,islc)); break;
            case TruthClass::kCosmic:    vector_active.push_back(-9999); break;
            default:
                vector_active.push_back(-9999);
                break;

        }
 }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});



const SpillMultiVar kMuon_trackScore([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        vector_active.push_back(islc.reco.pfp[ipfp_mu].trackScore);
         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kProton_trackScore([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        vector_active.push_back(islc.reco.pfp[longest_proton_idx].trackScore);
        }}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});







const SpillMultiVar kVertex_x_true([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        if(!std::isnan(islc.truth.position.x))vector_active.push_back(islc.truth.position.x);
        if(std::isnan(islc.truth.position.x))vector_active.push_back(-9999);

 }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kVertex_y_true([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        if(!std::isnan(islc.truth.position.y))vector_active.push_back(islc.truth.position.y);
        if(std::isnan(islc.truth.position.y))vector_active.push_back(-9999);
 }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kVertex_z_true([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        if(!std::isnan(islc.truth.position.z))vector_active.push_back(islc.truth.position.z);
        if(std::isnan(islc.truth.position.z))vector_active.push_back(-9999);

 }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kMuon_endx_true([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        if(!std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.end.x))vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.end.x);
        if(std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.end.x))vector_active.push_back(-9999);

         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kMuon_endy_true([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        if(!std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.end.y))vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.end.y);
        if(std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.end.y))vector_active.push_back(-9999);         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kMuon_endz_true([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        if(!std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.end.z))vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.truth.p.end.z);
        if(std::isnan(islc.reco.pfp[ipfp_mu].trk.truth.p.end.z))vector_active.push_back(-9999);         }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});






const SpillMultiVar kEff_raster_angle([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for ( auto const& nu : sr->mc.nu ){
    // ----------------- Sanity -----------------
    TruthClass cls_nu = classification_type_MC(sr,nu);
    if(cls_nu == TruthClass::kOneMuOneP || cls_nu == TruthClass::kOneMuNp){

        int maxCut = 1; // default: total only
        double angle_3D_mup = -9;

        for (auto const& islc : sr->slc){
            int ipfp_mu = -1; 

            if(islc.tmatch.eff<0.5)continue;
            if(islc.truth.index != nu.index) continue;
            
            CutResults res = automatic_selection_1muNp_new_SBND_cutflow(sr, islc, 10, 100);

            //if(maxCut<MaxCutPassed(res)){ maxCut = MaxCutPassed(res);}


                if(maxCut<MaxCutPassed(res)){
                    maxCut = MaxCutPassed(res);
                    ipfp_mu = find_muon(islc, 10);
                    if(ipfp_mu!=-1){
                    auto longest_proton = find_longest_proton(islc, ipfp_mu, 10);
                    int longest_proton_idx = longest_proton.index;
                    if (longest_proton_idx >= 0) {
                    angle_3D_mup =T3D_angle_mup(islc,ipfp_mu,longest_proton_idx);
                    }}                    
                    
                    }

        } //loop reco slice
        
        vector_active.push_back(angle_3D_mup); 

    }     
        
    }


 	return vector_active;
});



const SpillMultiVar kNu_score([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new_SBND(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

        vector_active.push_back(islc.nu_score); }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});
