//#include "ReadTree.C"
std::ofstream dump("BrokenSimulDump_YZvariations.txt");


//std::vector<std::vector<int>> HITmultMU;
//std::vector<std::vector<double>> HITpitchMU;
//std::vector<std::vector<double>> HITxMU;
//std::vector<std::vector<double>> HITyMU;
//std::vector<std::vector<double>> HITzMU;
TRandom *r3 = new TRandom3();
void FillCaloInfox(EventsData &dat, int track, bool broken, std::vector<std::vector<double>> &HITdeMU, std::vector<std::vector<double>> &HITrrMU )
{ 
        std::vector<double> HITde;
        std::vector<double> HITr;
        //std::vector<double> HITpitch;
        //std::vector<int> HITmult;
        //std::vector<double> HITx;
        //std::vector<double> HITy;
        //std::vector<double> HITz;

        double offset_rr=0;
		double x_broken_position=10e6;

        dat.tree->GetEntry(track);

        if(broken==false){offset_rr=0;}
        else if(broken==true)
        {
            double best_distance = 10e6;
            for (int hit=0; hit<int(dat.track.rr->size()); hit++)
            {
                            if(fabs(dat.track.start_reco->at(0))>210.07)
                            {
                                if(fabs(dat.track.hitx->at(hit))>210.07)
                                {
                                    double distance = fabs(fabs(dat.track.hitx->at(hit))-210.07);
                                    if(distance <= best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  dat.track.rr->at(hit);
										x_broken_position = dat.track.hitx->at(hit);
                                    }
                                    
                                }
                            }
                            else if(fabs(dat.track.start_reco->at(0))<210.07)
                            {
                                if(fabs(dat.track.hitx->at(hit))<210.07)
                                {
                                    double distance = fabs(fabs(dat.track.hitx->at(hit))-210.07);
                                    if(distance < best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  dat.track.rr->at(hit);
										x_broken_position = dat.track.hitx->at(hit);
                                    }
                                    
                                }
                            }
                        
            }
			dump << x_broken_position << " " << offset_rr << " " << best_distance << " ";
        }

        for ( int hit=0; hit<int(dat.track.rr->size()); hit++ )
            {   

                            HITde.push_back(dat.track.dE->at(hit));           

                            //HITpitch.push_back(dat.track.pitch->at(hit));

                            HITr.push_back(dat.track.rr->at(hit)-offset_rr);

                            //HITmult.push_back(dat.track.mult->at(hit));

                            //HITx.push_back(dat.track.hitx->at(hit));

                            //HITy.push_back(dat.track.hity->at(hit));

                            //HITz.push_back(dat.track.hitz->at(hit));

            }//loop over all hits of that particle

            HITdeMU.push_back(HITde);
            HITrrMU.push_back(HITr);
            //HITpitchMU.push_back(HITpitch);
            //HITmultMU.push_back(HITmult);
            //HITxMU.push_back(HITx);
            //HITyMU.push_back(HITy);
            //HITzMU.push_back(HITz);


};

void FillCaloInfoz(EventsData &dat, int track, bool broken, std::vector<std::vector<double>> &HITdeMU, std::vector<std::vector<double>> &HITrrMU)
{ 
        std::vector<double> HITde;
        std::vector<double> HITr;
        //std::vector<double> HITpitch;
        //std::vector<int> HITmult;
        //std::vector<double> HITx;
        //std::vector<double> HITy;
        //std::vector<double> HITz;

        dat.tree->GetEntry(track);
        double offset_rr=0;
		double z_broken_position=10e6;

        if(broken==false){offset_rr=0;}
        else if(broken==true)
        {
            double best_distance = 10e6;
            for ( int hit=0; hit<int(dat.track.rr->size()); hit++ )
            {
                            if(dat.track.start_reco->at(2)>0.)
                            {
                                if(dat.track.hitz->at(hit)>0.)
                                {
                                    double distance = fabs(dat.track.hitz->at(hit));
                                    if(distance < best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  dat.track.rr->at(hit);
										z_broken_position = dat.track.hitz->at(hit);
                                    }
                                    
                                }
                            }
                            else if(dat.track.start_reco->at(2)<0.)
                            {
                                if(dat.track.hitz->at(hit)<0.)
                                {
                                    double distance = fabs(dat.track.hitz->at(hit));
                                    if(distance < best_distance)
                                    {
                                        best_distance=distance;
                                        offset_rr =  dat.track.rr->at(hit);
										z_broken_position = dat.track.hitz->at(hit);
                                    }
                                    
                                }
                            }
            }
			dump << z_broken_position << " " << offset_rr << " " << best_distance << " ";
        }

        for ( int hit=0; hit<int(dat.track.rr->size()); hit++ )
            {   

                            HITde.push_back(dat.track.dE->at(hit));           

                            //HITpitch.push_back(dat.track.pitch->at(hit));

                            HITr.push_back(dat.track.rr->at(hit)-offset_rr);

                            //HITmult.push_back(dat.track.mult->at(hit));

                            //HITx.push_back(dat.track.hitx->at(hit));

                            //HITy.push_back(dat.track.hity->at(hit));

                            //HITz.push_back(dat.track.hitz->at(hit));

            }//loop over all hits of that particle

            HITdeMU.push_back(HITde);
            HITrrMU.push_back(HITr);
            //HITpitchMU.push_back(HITpitch);
            //HITmultMU.push_back(HITmult);
            //HITxMU.push_back(HITx);
            //HITyMU.push_back(HITy);
            //HITzMU.push_back(HITz);

};





std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> simul_broken_xz(EventsData dat)
{

    std::vector<double> endx;
    std::vector<double> endz;
    std::vector<std::vector<double>> HITdeMU;
    std::vector<std::vector<double>> HITrrMU;

    double start_z_muon, end_z_muon, start_x_muon, end_x_muon, l_muon, start_x_muon_t, end_x_muon_t;
    int flag=-1;
    for (int track=0; track<int(dat.tree->GetEntries()); track++)
    {     

      dat.tree->GetEntry(track);  

      flag=-1;

            start_z_muon=dat.track.start_reco->at(2);
            end_z_muon=dat.track.end_reco->at(2);
            l_muon=dat.track.len_reco;
            start_x_muon=fabs(dat.track.start_reco->at(0));
            end_x_muon=fabs(dat.track.end_reco->at(0));
            start_x_muon_t=dat.track.start_reco->at(0);
            end_x_muon_t=dat.track.end_reco->at(0);
					
			dump << start_x_muon_t << " " << end_x_muon_t << " " << start_z_muon << " " << end_z_muon << " "; 
                    
            //se start e end sono nella stessa TPC
            if(start_x_muon<210.07 && end_x_muon<210.07) 
            {
                endx.push_back(dat.track.end_reco->at(0));
                FillCaloInfox(dat,track,false,HITdeMU,HITrrMU);
				dump << 0 << " -1 -1 -1 " ;
            }
            else if(start_x_muon>210.07 && end_x_muon>210.07) 
            {
                endx.push_back(dat.track.end_reco->at(0));
                FillCaloInfox(dat,track,false,HITdeMU,HITrrMU);
				dump << 0 << " -1 -1 -1 ";
            }
            //se sono in TPC diverse
            else
            {
	            float random=r3->Rndm();
                //condizione diversa per EAST E WEST
                if((random>1161.*0.33/6133. && dat.track.end_reco->at(0)<0.) || (random>1161.*0.11/6133. && dat.track.end_reco->at(0)>0.))
                {
                    endx.push_back(dat.track.end_reco->at(0));
                    FillCaloInfox(dat,track,false,HITdeMU,HITrrMU);
					dump << 0 << " -1 -1 -1 ";
                }
	            else
	            {
                    flag=1;
	                float new_l_mu=l_muon*fabs(start_x_muon-210.07)/fabs(end_x_muon-start_x_muon);
                    //se viene considerato segnale
	                if(new_l_mu>50.)
	                {
		                endx.push_back(210.*dat.track.end_reco->at(0)/end_x_muon);
                        dump << 1 << " ";
						FillCaloInfox(dat,track,true,HITdeMU,HITrrMU);
		            }
					else if(new_l_mu<=50.){dump << 2 << " -1 -1 -1 ";}
	            }
        }//cathode
                    

        if(flag==-1){

            //se start e end sono dalla stessa parte
            if(start_z_muon<0. && end_z_muon<0.) 
            {
                endz.push_back(end_z_muon);
                FillCaloInfoz(dat,track,false,HITdeMU,HITrrMU);
				dump << 0 << " -1 -1 -1 ";
            }
            else if(start_z_muon>0. && end_z_muon>0.) 
            {
                endz.push_back(end_z_muon);
                FillCaloInfoz(dat,track,false,HITdeMU,HITrrMU);
				dump << 0 << " -1 -1 -1 ";
			}
            //se sono in parti diverse
            else
            {
	            float random=r3->Rndm();
	            if(random>268./1397.)
                {
                    endz.push_back(end_z_muon);
                    FillCaloInfoz(dat,track,false,HITdeMU,HITrrMU);
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
					    FillCaloInfoz(dat,track,true,HITdeMU,HITrrMU);
		            }
					else if(new_l_mu<=50.){dump << 2 << " -1 -1 -1 ";}

	            }
            }//z0
        }// se non è rotta al catodo e la lunghezza rimasta > 50 cm
        else if(flag==1){dump << "3 -1 -1 -1 ";}
  
		dump << endl;
            
    }
    
    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> broken_calo;
    broken_calo.first=HITdeMU;
    broken_calo.second=HITrrMU;

    return broken_calo;

};