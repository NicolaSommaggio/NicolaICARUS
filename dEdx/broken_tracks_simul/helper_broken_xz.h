

std::vector<double> HITde;
std::vector<double> HITr;
std::vector<double> HITpitch;
std::vector<int> HITmult;
std::vector<double> HITx;
std::vector<double> HITy;
std::vector<double> HITz;

std::vector<std::vector<double>> HITdeMU;
std::vector<std::vector<double>> HITrrMU;
std::vector<std::vector<int>> HITmultMU;
std::vector<std::vector<double>> HITpitchMU;
std::vector<std::vector<double>> HITxMU;
std::vector<std::vector<double>> HITyMU;
std::vector<std::vector<double>> HITzMU;

std::ofstream dump("BrokenSimulDump_YZvariations.txt");

void FillCaloInfox(const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, bool broken)
{ 
        std::vector<double> vector_active;
        double offset_rr=0;
		double x_broken_position=10e6;

        if(broken==false){offset_rr=0;}
        else if(broken==true)
        {
            double best_distance = 10e6;
            for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit )
            {
                    if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx<500.)
                        {
                            if(fabs(islc.reco.pfp[ipfp_mu].trk.start.x)>210.07)
                            {
                                if(fabs(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x)>210.07)
                                {
                                    double distance = fabs(fabs(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x)-210.07);
                                    if(distance <= best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr;
										x_broken_position = islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x;
                                    }
                                    
                                }
                            }
                            else if(fabs(islc.reco.pfp[ipfp_mu].trk.start.x)<210.07)
                            {
                                if(fabs(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x)<210.07)
                                {
                                    double distance = fabs(fabs(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x)-210.07);
                                    if(distance < best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr;
										x_broken_position = islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x;
                                    }
                                    
                                }
                            }
                        }
            }
			dump << x_broken_position << " " << offset_rr << " " << best_distance << " ";
        }

        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit )
            {
                    if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx<500.)
                        {    

                            HITde.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx);           

                            HITpitch.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].pitch);

                            HITr.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr-offset_rr);

                            HITmult.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].mult);

                            HITx.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x);
                            vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x);

                            HITy.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].y);

                            HITz.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z);
                        }

            }//loop over all hits of that particle

            HITdeMU.push_back(HITde);
            HITrrMU.push_back(HITr);
            HITpitchMU.push_back(HITpitch);
            HITmultMU.push_back(HITmult);
            HITxMU.push_back(HITx);
            HITyMU.push_back(HITy);
            HITzMU.push_back(HITz);

            
            HITx.clear();
            HITy.clear();
            HITz.clear();
            HITmult.clear();
            HITpitch.clear();
            HITde.clear();
            HITr.clear();

};

void FillCaloInfoz(const caf::Proxy<caf::SRSlice>& islc, int ipfp_mu, bool broken)
{ 
        std::vector<double> vector_active;
        double offset_rr=0;
		double z_broken_position=10e6;

        if(broken==false){offset_rr=0;}
        else if(broken==true)
        {
            double best_distance = 10e6;
            for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit )
            {
                    if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx<500.)
                        {
                            if(islc.reco.pfp[ipfp_mu].trk.start.z>0.)
                            {
                                if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z>0.)
                                {
                                    double distance = fabs(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z);
                                    if(distance < best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr;
										z_broken_position = islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z;
                                    }
                                    
                                }
                            }
                            else if(islc.reco.pfp[ipfp_mu].trk.start.z<0.)
                            {
                                if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z<0.)
                                {
                                    double distance = fabs(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z);
                                    if(distance < best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr;
										z_broken_position = islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z;
                                    }
                                    
                                }
                            }
                        }
            }
			dump << z_broken_position << " " << offset_rr << " " << best_distance << " ";
        }

        for ( std::size_t ihit(0); ihit < islc.reco.pfp[ipfp_mu].trk.calo[2].points.size(); ++ihit )
            {
                    if(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx>0.5 && islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx<500.)
                        {    

                            HITde.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].dedx);           

                            HITpitch.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].pitch);

                            HITr.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].rr-offset_rr);

                            HITmult.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].mult);

                            HITx.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].x);

                            HITy.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].y);

                            HITz.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z);
                            vector_active.push_back(islc.reco.pfp[ipfp_mu].trk.calo[2].points[ihit].z);
                        }

            }//loop over all hits of that particle

            HITdeMU.push_back(HITde);
            HITrrMU.push_back(HITr);
            HITpitchMU.push_back(HITpitch);
            HITmultMU.push_back(HITmult);
            HITxMU.push_back(HITx);
            HITyMU.push_back(HITy);
            HITzMU.push_back(HITz);

            
            HITx.clear();
            HITy.clear();
            HITz.clear();
            HITmult.clear();
            HITpitch.clear();
            HITde.clear();
            HITr.clear();

};




std::vector<double> endx;
std::vector<double> endz;
const SpillMultiVar simul_broken_xz([](const caf::SRSpillProxy* sr)-> std::vector<double>
{

    double start_z_muon, end_z_muon, start_x_muon, end_x_muon, l_muon, start_x_muon_t, end_x_muon_t;
    int flag=-1;
    for (auto const& islc : sr->slc)
    {     
      int ipfp_mu=-1;
      int ipfp_pro=-1;
      double distanza;

      
        if (automatic_selection_1muNp(sr,islc,10,100))//1muNp reco   _mod per 1mu1p
	        {
            ipfp_mu=find_muon(islc,10); //trovo l'indice del muone tra tutte le particelle nella slice

                //if(!(islc.truth.index>=0 && (classification_type(sr,islc)==2 || classification_type(sr,islc)==5 ))) {continue;} //1muNp E>50 true

                if(ipfp_mu!=-1)
                {  
                    flag=-1;

                    start_z_muon=islc.reco.pfp[ipfp_mu].trk.start.z;
                    end_z_muon=islc.reco.pfp[ipfp_mu].trk.end.z;
                    l_muon=islc.reco.pfp[ipfp_mu].trk.len;
                    start_x_muon=fabs(islc.reco.pfp[ipfp_mu].trk.start.x);
                    end_x_muon=fabs(islc.reco.pfp[ipfp_mu].trk.end.x);
                    start_x_muon_t=islc.reco.pfp[ipfp_mu].trk.start.x;
                    end_x_muon_t=islc.reco.pfp[ipfp_mu].trk.end.x;
					
					dump << start_x_muon_t << " " << end_x_muon_t << " " << start_z_muon << " " << end_z_muon << " "; 

                    
                    //se start e end sono nella stessa TPC
                    if(start_x_muon<210.07 && end_x_muon<210.07) 
                    {
                        endx.push_back(islc.reco.pfp[ipfp_mu].trk.end.x);
                        FillCaloInfox(islc,ipfp_mu,false);
						dump << 0 << " -1 -1 -1 " ;
                    }
                    else if(start_x_muon>210.07 && end_x_muon>210.07) 
                    {
                        endx.push_back(islc.reco.pfp[ipfp_mu].trk.end.x);
                        FillCaloInfox(islc,ipfp_mu,false);
						dump << 0 << " -1 -1 -1 ";
                    }
                    //se sono in TPC diverse
                    else
                    {
	                    float random=r3->Rndm();
                        //condizione diversa per EAST E WEST
                        if((random>1161.*0.33/6133. && islc.reco.pfp[ipfp_mu].trk.end.x<0.) || (random>1161.*0.11/6133. && islc.reco.pfp[ipfp_mu].trk.end.x>0.))
                        {
                            endx.push_back(islc.reco.pfp[ipfp_mu].trk.end.x);
                            FillCaloInfox(islc,ipfp_mu,false);
							dump << 0 << " -1 -1 -1 ";
                        }
	                    else
	                    {
                            flag=1;
	                        float new_l_mu=l_muon*fabs(start_x_muon-210.07)/fabs(end_x_muon-start_x_muon);
                            //se viene considerato segnale
	                        if(new_l_mu>50.)
	                        {
		                        endx.push_back(210.*islc.reco.pfp[ipfp_mu].trk.end.x/end_x_muon);
                                dump << 1 << " ";
								FillCaloInfox(islc,ipfp_mu,true);
		                    }
							else if(new_l_mu<=50.){dump << 2 << " -1 -1 -1 ";}
	                    }
                    }//cathode
                    

                    if(flag==-1){

                    //se start e end sono dalla stessa parte
                    if(start_z_muon<0. && end_z_muon<0.) 
                    {
                        endz.push_back(end_z_muon);
                        FillCaloInfoz(islc,ipfp_mu,false);
						dump << 0 << " -1 -1 -1 ";
                    }
                    else if(start_z_muon>0. && end_z_muon>0.) 
                    {
                        endz.push_back(end_z_muon);
                        FillCaloInfoz(islc,ipfp_mu,false);
						dump << 0 << " -1 -1 -1 ";
					}
                    //se sono in parti diverse
                    else
                    {
	                    float random=r3->Rndm();
	                    if(random>268./1397.)
                        {
                            endz.push_back(end_z_muon);
                            FillCaloInfoz(islc,ipfp_mu,false);
							dump << 0 << " -1 -1 -1 ";
                        }
	                    else
	                    {
	                        float new_l_mu=l_muon*fabs(start_z_muon)/fabs(end_z_muon-start_z_muon);
                            //se viene considerato segnale
	                        if(new_l_mu>50.)
	                        {
		                        endz.push_back(0.);
                                dump << 1 << " ";
								FillCaloInfoz(islc,ipfp_mu,true);
		                    }
							else if(new_l_mu<=50.){dump << 2 << " -1 -1 -1 ";}

	                    }
                    }//z0
                    }// se non è rotta al catodo e la lunghezza rimasta > 50 cm
                    else if(flag==1){dump << "3 -1 -1 -1 ";}



                }
				dump << endl;
            }
    }
    

    return endx;

});