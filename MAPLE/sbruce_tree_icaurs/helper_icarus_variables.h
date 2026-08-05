#include "helper_ICARUS.h"
//#include "helper_ICARUS_calosys.h"

struct LongestProtonResult {
    int index = -1;
    double length = 0.0;
};


LongestProtonResult find_longest_proton(const caf::SRSliceProxy* slc,
                                        int muon_index,
                                        double dist_cut)
{
    LongestProtonResult result;

    for (std::size_t ipfp = 0; ipfp < slc->reco.npfp; ++ipfp) {

        if (static_cast<int>(ipfp) == muon_index) {
            continue;
        }

        PFPID pid = static_cast<PFPID>(id_pfp(slc, ipfp, dist_cut));

        if (pid != PFPID::Proton) {
            continue;
        }

        double length = slc->reco.pfp[ipfp].trk.len;  // adjust if needed

        if (length > result.length) {
            result.length = length;
            result.index  = static_cast<int>(ipfp);
        }
    }

    return result;
}


double Transverse_angle ( const caf::SRSliceProxy* slc, int muon_index, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.x; //GeV
    p_mu_y=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.y;
    p_mu_z=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV

    //Proton momentum
    p_p_x =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.z;
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

double T3D_angle_mup ( const caf::SRSliceProxy* slc, int muon_index, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;


    //Muon momentum
    p_mu_x=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.x; //GeV
    p_mu_y=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.y;
    p_mu_z=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.z;
    TVector3 mu_vector(p_mu_x,p_mu_y,p_mu_z);

    //Proton momentum
    p_p_x =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.z;
    TVector3 pro_vector(p_p_x,p_p_y,p_p_z);

     
    return cos(mu_vector.Angle(pro_vector)); 
    }


double Transverse_mom_reco ( const caf::SRSliceProxy* slc, int muon_index, int ipfp_pro ) {
    
    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    float p_tot_x,p_tot_y,p_tot_z;
    float E_mu,E_p;    
    //Muon momentum
    p_mu_x=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.x; //GeV
    p_mu_y=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.y;
    p_mu_z=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.z;
    double p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV
    E_mu=sqrt(p_mu_tot*p_mu_tot+(105.658*105.658)/(1000*1000));

    //Proton momentum
    p_p_x =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.z;
    double p_p_tot = sqrt(p_p_x*p_p_x+p_p_y*p_p_y+p_p_z*p_p_z);
    E_p=sqrt(p_p_tot*p_p_tot+(938.272*938.272)/(1000*1000));

    //Total momentum
    p_tot_x=p_p_x+p_mu_x;
    p_tot_y=p_p_y+p_mu_y;
    p_tot_z=p_p_z+p_mu_z;

    double p_T=sqrt(p_tot_x*p_tot_x+p_tot_y*p_tot_y);


return p_T;
}


double Transverse_mom_reco_mu ( const caf::SRSliceProxy* slc, int muon_index, int ipfp_pro ) {
    
    float p_mu_x,p_mu_y;   
    //Muon momentum
    p_mu_x=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.x; //GeV
    p_mu_y=(slc->reco.pfp[muon_index].trk.rangeP.p_muon)*slc->reco.pfp[muon_index].trk.dir.y;


    double p_T=sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y);


return p_T;
}

double Transverse_mom_reco_pro ( const caf::SRSliceProxy* slc, int muon_index, int ipfp_pro ) {
    
    float p_p_x,p_p_y;

    //Proton momentum
    p_p_x =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(slc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*slc->reco.pfp[ipfp_pro].trk.dir.y;


    double p_T=sqrt(p_p_x*p_p_x+p_p_y*p_p_y);


return p_T;
}




const Var kTransverse_mom_reco([](const caf::SRSliceProxy* slc){
    double Pt = 0;
    
    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        Pt = Transverse_mom_reco_Np(slc,muon_index);
        
    }//selected by the automatic selection
    

 	return Pt;
});


const Var kProton_length_leading([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;

    if (longest_proton_idx >= 0) {
        res = slc->reco.pfp[longest_proton_idx].trk.len;
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kMuon_length([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){ 
        const int muon_index = find_muon(slc, 10);       
        res = slc->reco.pfp[muon_index].trk.len;
         
        
    }//selected by the automatic selection
    
    

 	return res;
});
const Var kProton_kinetic_leading([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
        if (longest_proton_idx >= 0) {
        
        TVector3 p_pr(slc->reco.pfp[longest_proton_idx].trk.dir.x * slc->reco.pfp[longest_proton_idx].trk.rangeP.p_proton,
              slc->reco.pfp[longest_proton_idx].trk.dir.y * slc->reco.pfp[longest_proton_idx].trk.rangeP.p_proton,
              slc->reco.pfp[longest_proton_idx].trk.dir.z * slc->reco.pfp[longest_proton_idx].trk.rangeP.p_proton);

        double ke = kinetic_energy(PROTON_MASS, p_pr.Mag());
        res = ke/1000.0;
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});



const Var kProton_chi2mu([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
        if (longest_proton_idx >= 0) {
    std::vector<double> chi2 = compute_chi2(slc, longest_proton_idx, DEFAULT_PLANE);
        res = chi2[static_cast<int>(Chi2PID::Muon)];
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});

const Var kProton_chi2pro([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
    std::vector<double> chi2 = compute_chi2(slc, longest_proton_idx, DEFAULT_PLANE);
        res = chi2[static_cast<int>(Chi2PID::Proton)];
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kMuon_chi2mu([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
    std::vector<double> chi2 = compute_chi2(slc, muon_index, DEFAULT_PLANE);
        res = chi2[static_cast<int>(Chi2PID::Muon)];
         
        
    }//selected by the automatic selection
    
    

 	return res;
});

const Var kMuon_chi2pro([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
    std::vector<double> chi2 = compute_chi2(slc, muon_index, DEFAULT_PLANE);
        res = chi2[static_cast<int>(Chi2PID::Proton)];
         
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kVertex_x([](const caf::SRSliceProxy* slc)
{
    double res = 0;


    if(Selection_1muNp_slice(slc)){        
        res = slc->vertex.x; 
        
    }//selected by the automatic selection
    
    

 	return res;
});

const Var kVertex_y([](const caf::SRSliceProxy* slc)
{
    double res = 0;


    if(Selection_1muNp_slice(slc)){
        res = slc->vertex.y; 
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kVertex_z([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
    if(Selection_1muNp_slice(slc)){        
        res = slc->vertex.z; 
        
    }//selected by the automatic selection
    
    

 	return res;
});




const Var kMuon_endx([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){ 
        const int muon_index = find_muon(slc, 10);       
        res = slc->reco.pfp[muon_index].trk.end.x;
         
        
    }//selected by the automatic selection
    
    

 	return res;
});

const Var kMuon_endy([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        res = slc->reco.pfp[muon_index].trk.end.y;
         
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kMuon_endz([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        res = slc->reco.pfp[muon_index].trk.end.z;
         
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kProton_endx([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = slc->reco.pfp[longest_proton_idx].trk.end.x;
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});

const Var kProton_endy([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = slc->reco.pfp[longest_proton_idx].trk.end.y;
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kProton_endz([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = slc->reco.pfp[longest_proton_idx].trk.end.z;
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kTransverse_angle([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = Transverse_angle(slc,muon_index,longest_proton_idx);
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kT3D_angle_mup([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = T3D_angle_mup(slc,muon_index,longest_proton_idx);
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kTransverse_mom_reco_mu([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = Transverse_mom_reco_mu(slc,muon_index,longest_proton_idx);
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});

const Var kTransverse_mom_reco_pro([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = Transverse_mom_reco_pro(slc,muon_index,longest_proton_idx);
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kNumber_protons([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
            ParticleCounts counts = count_particles(slc, muon_index, 10);
        res = counts.protons;
        
        
    }//selected by the automatic selection
    
    

 	return res;
});



const Var kMuon_trackScore([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        res = slc->reco.pfp[muon_index].trackScore;
         
        
    }//selected by the automatic selection
    
    

 	return res;
});


const Var kProton_trackScore([](const caf::SRSliceProxy* slc)
{
    double res = 0;

    
        int muon_index = -1;    

    if(Selection_1muNp_slice(slc)){
        const int muon_index = find_muon(slc, 10);
        
        auto longest_proton = find_longest_proton(slc, muon_index, 10);
        int longest_proton_idx = longest_proton.index;
    if (longest_proton_idx >= 0) {
        res = slc->reco.pfp[longest_proton_idx].trackScore;
        }
        
    }//selected by the automatic selection
    
    

 	return res;
});


