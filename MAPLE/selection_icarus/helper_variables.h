//#include "helper/helper_eff_cf.h"
//#include "helper/helper_eff_cf_Lmu_EKpro_CHI2var_TRKSCOREvar_NOmup_only_CRTveto_Lmu_FV_Trigger.h"
#include "helper/helper_eff_cf_Lmu_EKpro_CHI2var_TRKSCOREvar_mup_only_CRTveto_Lmu_FV_Trigger.h"
//#include "helper/helper_eff_cf_Lmu_EKpro_CHI2var_TRKSCOREvar_mup_only_CRTveto_Lmu_FV.h"
//#include "helper/helper_eff_cf_Lmu_EKpro_CHI2var_TRKSCOREvar_mup_only_CRTveto_Lmu.h"
//#include "helper/helper_eff_cf_Lmu_EKpro_CHI2var_TRKSCOREvar_mup_only_CRTveto.h"
//#include "helper/helper_eff_cf_Lmu_EKpro_CHI2var_TRKSCOREvar_mup_only.h"
//#include "helper_eff_cf_Lmu_EKpro_CHI2var_TRKSCOREvar.h"
//#include "helper_eff_cf_Lmu_EKpro_CHI2var.h"
//#include "helper_eff_cf_Lmu_EKpro.h"
//#include "helper_variations_calo.h"

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


const SpillMultiVar kMuon_length([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
            vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.len);}
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kProton_length_leading([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

        vector_active.push_back(islc.vertex.z); }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});


const SpillMultiVar kMuon_endx([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
            ParticleCounts counts = count_particles(islc, ipfp_mu, 10);
        vector_active.push_back(counts.protons);
        }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar kBarycenter_delta([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

      double barycenter=bar_flash(sr,islc);  //This is z!!!
	  double barycenter_x=bar_flash_x(sr,islc);

	  double bar_charge=islc.charge_center.z;
	  double delta=fabs(barycenter-bar_charge);
      if(barycenter<-9999 && abs(barycenter_x)<1){vector_active.push_back(-10);}
      if(!(barycenter<-9999 && abs(barycenter_x)<1)){vector_active.push_back(delta);}
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
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
            
            CutResults res = automatic_selection_1muNp_new_cutflow(sr, islc, 10, 100);

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

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){

        vector_active.push_back(islc.nu_score); }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});






const SpillMultiVar kEff_raster_nuscore([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for ( auto const& nu : sr->mc.nu ){
    // ----------------- Sanity -----------------
    TruthClass cls_nu = classification_type_MC(sr,nu);
    if(cls_nu == TruthClass::kOneMuOneP || cls_nu == TruthClass::kOneMuNp){

        int maxCut = 1; // default: total only
        double nu_score = -9;

        for (auto const& islc : sr->slc){
            int ipfp_mu = -1; 

            if(islc.tmatch.eff<0.5)continue;
            if(islc.truth.index != nu.index) continue;
            
            CutResults res = automatic_selection_1muNp_new_cutflow(sr, islc, 10, 100);

            //if(maxCut<MaxCutPassed(res)){ maxCut = MaxCutPassed(res);}


                if(maxCut<MaxCutPassed(res)){
                    maxCut = MaxCutPassed(res);
                    ipfp_mu = find_muon(islc, 10);
                    if(ipfp_mu!=-1){
                       nu_score = islc.nu_score;
                       }                    
                    
                    }

        } //loop reco slice
        
        vector_active.push_back(nu_score); 

    }     
        
    }


 	return vector_active;
});





const SpillMultiVar kCryoSel([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;    

    if(automatic_selection_1muNp_new(sr, islc,10,100 )){
        ipfp_mu = find_muon(islc, 10);
        if(ipfp_mu!=-1){
        int cryo_light = cryo_selection_from_light(sr);
        if (cryo_light < 0)
        {
            vector_active.push_back(-1);
            continue;
        }

        bool light_in_east = (cryo_light == 2 || cryo_light == 0);
        bool light_in_west = (cryo_light == 2 || cryo_light == 1);
        // x > 0 WEST   x < 0 EAST !!!!

        if (islc.vertex.x > 0 && light_in_west){ vector_active.push_back(1);}
        else if (islc.vertex.x < 0 && light_in_east){vector_active.push_back(0);}
        else {vector_active.push_back(-1);}

       }
        
    }//selected by the automatic selection
    
    }//loop over slices

 	return vector_active;
});



const SpillMultiVar kEff_cryosel([](const caf::SRSpillProxy* sr) -> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& nu : sr->mc.nu)
    {
        // ----------------- Sanity -----------------
        TruthClass cls_nu = classification_type_MC(sr, nu);
        if (cls_nu != TruthClass::kOneMuOneP && cls_nu != TruthClass::kOneMuNp)
            continue;

        int maxCut = 1; // default: total only
        bool found_best_slice = false;
        double best_vertex_x = 0;

        for (auto const& islc : sr->slc)
        {
            if (islc.tmatch.eff < 0.5) continue;
            if (islc.truth.index != nu.index) continue;

            CutResults res = automatic_selection_1muNp_new_cutflow(sr, islc, 10, 100);

            if (maxCut < MaxCutPassed(res))
            {
                maxCut = MaxCutPassed(res);
                int ipfp_mu = find_muon(islc, 10);
                if (ipfp_mu != -1)
                {
                    best_vertex_x = islc.vertex.x;
                    found_best_slice = true;
                }
            }
        } // loop reco slice

        if (!found_best_slice)
        {
            vector_active.push_back(-1);
            continue;
        }

        int cryo_light = cryo_selection_from_light(sr);
        if (cryo_light < 0)
        {
            vector_active.push_back(-1);
            continue;
        }

        bool light_in_east = (cryo_light == 2 || cryo_light == 0);
        bool light_in_west = (cryo_light == 2 || cryo_light == 1);
        // x > 0 WEST   x < 0 EAST !!!!

        if (best_vertex_x > 0 && light_in_west)
        {
            vector_active.push_back(1);
        }
        else if (best_vertex_x < 0 && light_in_east)
        {
            vector_active.push_back(0);
        }
        else
        {
            vector_active.push_back(-1);
        }

    } // loop over nu

    return vector_active;
});




const SpillMultiVar k_E_nu_true([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    for (auto const& islc : sr->slc)
    {   
        if(automatic_selection_1muNp_new(sr, islc,10,100 ))
        {    
            vector_active.push_back(islc.truth.E);
        }//selected by the automatic selection
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar k_pe([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;
    constexpr double kMinTime = -0.6;
    constexpr double kMaxTime = 1.8;

    for (auto const& islc : sr->slc)
    {   
        if(automatic_selection_1muNp_new(sr, islc,10,100 ))
        {    
            double max_pe = 0;
            for (const auto& match : sr->crtpmt_matches)
            {
                if (match.flashGateTime <= kMinTime || match.flashGateTime >= kMaxTime)continue;

                if(match.flashPE > max_pe){max_pe = match.flashPE;}
            }
            vector_active.push_back(max_pe);
        }//selected by the automatic selection
    }//loop over slices

 	return vector_active;
});

const SpillMultiVar k_slice_indedx([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;

    double slice_index = -1;
    for (auto const& islc : sr->slc)
    {   
        slice_index = slice_index + 1;
        if(automatic_selection_1muNp_new(sr, islc,10,100 ))
        {    
            vector_active.push_back(slice_index);
        }//selected by the automatic selection
    }//loop over slices

 	return vector_active;
});
