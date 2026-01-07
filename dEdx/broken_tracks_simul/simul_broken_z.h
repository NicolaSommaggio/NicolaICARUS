const SpillMultiVar kMuon_endz_simul_MC_Np([](const caf::SRSpillProxy* sr)-> std::vector<double>
{
    std::vector<double> vector_active;
    float p_mu_x = 0; float p_mu_y=0; float p_mu_z=0; double p_mu_tot=0;
         float p_mu_x_a = 0; float p_mu_y_a=0; float p_mu_z_a=0;double p_mu_tot_a=0;
        float p_p_x=0; float p_p_y =0; float p_p_z=0;
        float p_tot_x =0; float p_tot_y=0; float p_tot_z=0;
        float p_tot_x_a =0; float p_tot_y_a=0; float p_tot_z_a=0;
         double E_mu =0 ; double E_p=0; double E_mu_a =0 ;double E_p_s=0;
         double E_nu =0 ; double E_nu_a =0 ;
	 double p_T=0; double p_T_a=0;
         float max_p_e=0;
          float max_p_l=0;
          float max_p_x=0;
          float max_p_y=0;
          float max_p_z=0;
          float max_p_dirx=0;
          float max_p_diry=0;
          float max_p_dirz=0;
          float mu_dirx=0;
          float mu_diry=0;
          float mu_dirz=0;

    float start_muon, end_muon, l_muon;
    for (auto const& islc : sr->slc){
        int ipfp_mu = -1;
        int ipfp_pro = -1;

    if(automatic_selection_1muNp(sr, islc,10,100)){

        ipfp_mu = automatic_selection_mu_index(sr,islc,10,100);
        //std::cout << " SCAN Run: " << sr->hdr.run << " Event: " << sr->hdr.evt << std::endl;
	if(ipfp_mu!=-1){

	  p_mu_x=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
	  p_mu_y=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.y;
	  p_mu_z=(islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc.reco.pfp[ipfp_mu].trk.dir.z;
	  mu_dirx=islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
	  mu_diry=islc.reco.pfp[ipfp_mu].trk.dir.y;
	  mu_dirz=islc.reco.pfp[ipfp_mu].trk.dir.z;
	  p_mu_tot = sqrt(p_mu_x*p_mu_x+p_mu_y*p_mu_y+p_mu_z*p_mu_z);             //GeV
	  E_mu=1000*sqrt(p_mu_tot*p_mu_tot+(105.658*105.658)/(1000*1000));
        
	  for ( std::size_t ipfp(0); ipfp < islc.reco.npfp ; ++ipfp ){
            if(int(ipfp)==ipfp_mu)continue;
            if(id_pfp(islc, ipfp, 10)==1){
	      TVector3 Start_mom_v2;
	      Start_mom_v2.SetXYZ((islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y,(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z);
	      E_p += (sqrt(pow(938.3,2)+pow(Start_mom_v2.Mag()*1000,2))-938.3);
	      E_p_s= (sqrt(pow(938.3,2)+pow(Start_mom_v2.Mag()*1000,2))-938.3);
	      p_p_x +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.x;
	      p_p_y +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.y;
	      p_p_z +=(islc.reco.pfp[ipfp].trk.rangeP.p_proton)*islc.reco.pfp[ipfp].trk.dir.z;
	      ipfp_pro=ipfp;
              if(E_p>max_p_e){max_p_e=E_p;max_p_l=islc.reco.pfp[ipfp].trk.len;max_p_x=islc.reco.pfp[ipfp].trk.end.x;max_p_y=islc.reco.pfp[ipfp].trk.end.y;max_p_z=islc.reco.pfp[ipfp].trk.end.z;max_p_dirx=islc.reco.pfp[ipfp].trk.dir.x;max_p_diry=islc.reco.pfp[ipfp].trk.dir.y;max_p_dirz=islc.reco.pfp[ipfp].trk.dir.z;}
            }
	  }

	  E_nu=(E_mu+E_p)/1000.;
            p_tot_x=p_p_x+p_mu_x;
            p_tot_y=p_p_y+p_mu_y;
            p_tot_z=p_p_z+p_mu_z;
            p_T = (sqrt(p_tot_x*p_tot_x+p_tot_y*p_tot_y));
	    E_nu_a=E_nu;
	    p_T_a=p_T;

      start_muon=islc.reco.pfp[ipfp_mu].trk.start.z;
      end_muon=islc.reco.pfp[ipfp_mu].trk.end.z;
      l_muon=islc.reco.pfp[ipfp_mu].trk.len;

      if(start_muon<0. && end_muon<0.) {vector_active.push_back(end_muon);MyFile<<E_nu<<" "<<E_nu_a<<" "<<p_T<<" "<<p_T_a<<" "<<end_muon<<" "<<end_muon<<" "<<start_muon<<" "<<l_muon<<" "<<l_muon<<" "<<p_mu_tot<<" "<<p_mu_tot<<" "<<mu_dirx<<" "<<mu_diry<<" "<<mu_dirz<<" "<<max_p_e<<" "<<max_p_l<<" "<<max_p_dirx<<" "<<max_p_diry<<" "<<max_p_dirz<<std::endl;}
      else if(start_muon>0. && end_muon>0.) {vector_active.push_back(end_muon);MyFile<<E_nu<<" "<<E_nu_a<<" "<<p_T<<" "<<p_T_a<<" "<<end_muon<<" "<<end_muon<<" "<<start_muon<<" "<<l_muon<<" "<<l_muon<<" "<<p_mu_tot<<" "<<p_mu_tot<<" "<<mu_dirx<<" "<<mu_diry<<" "<<mu_dirz<<" "<<max_p_e<<" "<<max_p_l<<" "<<max_p_dirx<<" "<<max_p_diry<<" "<<max_p_dirz<<std::endl;}
      else{
	float random=r3->Rndm();
	//std::cout << random << std::endl;
	//if(random>308./1349.){vector_active.push_back(end_muon);MyFile<<E_nu<<" "<<E_nu_a<<" "<<p_T<<" "<<p_T_a<<" "<<end_muon<<" "<<end_muon<<std::endl;}
	if(random>268./1397.){vector_active.push_back(end_muon);MyFile<<E_nu<<" "<<E_nu_a<<" "<<p_T<<" "<<p_T_a<<" "<<end_muon<<" "<<end_muon<<" "<<start_muon<<" "<<l_muon<<" "<<l_muon<<" "<<p_mu_tot<<" "<<p_mu_tot<<" "<<mu_dirx<<" "<<mu_diry<<" "<<mu_dirz<<" "<<max_p_e<<" "<<max_p_l<<" "<<max_p_dirx<<" "<<max_p_diry<<" "<<max_p_dirz<<std::endl;}
	else
	  {
	    float new_l_mu=l_muon*fabs(start_muon)/fabs(end_muon-start_muon);
	    if(new_l_mu>50.)
	      {
		vector_active.push_back(0.);
		float mom_new_mu=GetTrackMomentum(new_l_mu, 13);
		p_mu_x_a=(mom_new_mu)*islc.reco.pfp[ipfp_mu].trk.dir.x; //GeV
		p_mu_y_a=(mom_new_mu)*islc.reco.pfp[ipfp_mu].trk.dir.y;
		p_mu_z_a=(mom_new_mu)*islc.reco.pfp[ipfp_mu].trk.dir.z;
		p_mu_tot_a = sqrt(p_mu_x_a*p_mu_x_a+p_mu_y_a*p_mu_y_a+p_mu_z_a*p_mu_z_a);             //GeV
		E_mu_a=1000*sqrt(p_mu_tot_a*p_mu_tot_a+(105.658*105.658)/(1000*1000));
		E_nu_a=(E_mu_a+E_p)/1000.;
		p_tot_x_a=p_p_x+p_mu_x_a;
		p_tot_y_a=p_p_y+p_mu_y_a;
		p_tot_z_a=p_p_z+p_mu_z_a;
		p_T_a = (sqrt(p_tot_x_a*p_tot_x_a+p_tot_y_a*p_tot_y_a));
		MyFile<<E_nu<<" "<<E_nu_a<<" "<<p_T<<" "<<p_T_a<<" "<<end_muon<<" 0."<<" "<<start_muon<<" "<<l_muon<<" "<<new_l_mu<<" "<<p_mu_tot<<" "<<p_mu_tot_a<<" "<<mu_dirx<<" "<<mu_diry<<" "<<mu_dirz<<" "<<max_p_e<<" "<<max_p_l<<" "<<max_p_dirx<<" "<<max_p_diry<<" "<<max_p_dirz<<std::endl;
	      }
	    else
	      {
		MyFile<<E_nu<<" -1. "<<p_T<<" -1. "<<end_muon<<" -1001."<<" "<<start_muon<<" "<<l_muon<<" -1."<<" "<<p_mu_tot<<" -1. "<<mu_dirx<<" "<<mu_diry<<" "<<mu_dirz<<" "<<max_p_e<<" "<<max_p_l<<" "<<max_p_dirx<<" "<<max_p_diry<<" "<<max_p_dirz<<std::endl;
	      }
	  }
      }
    }

    }//selected by the automatic selection
    
    }//loop over slices



 	return vector_active;
});
