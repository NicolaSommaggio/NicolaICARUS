void Leggo()
{
  int runm[200];
  int evtm[200];
  float vtxm[200];
  float vtym[200];
  float vtzm[200];

  int aa,bb;
  float xx,yy,zz;

  ifstream ifile("selevts.txt");
  for(int ijk=0;ijk<168;ijk++)
    {
      ifile>>aa>>bb>>xx>>yy>>zz;
      runm[ijk]=aa;
      evtm[ijk]=bb;
      vtxm[ijk]=xx;
      vtym[ijk]=yy;
      vtzm[ijk]=zz;
    }

    int eventi_inpiu = 0;
    ifstream ifile2("tot_selected.txt");
    for(int ijk=0;ijk<168;ijk++)
    {
      ifile2>>aa>>bb>>xx>>yy>>zz;
      bool found = false;
      for(int ijk2=0;ijk2<168;ijk2++)
	    {
	      if(runm[ijk2]==aa && evtm[ijk2]==bb)
	      {
	        if(fabs(vtxm[ijk2]-xx)<1. && fabs(vtym[ijk2]-yy)<1. && fabs(vtzm[ijk2]-zz)<1.){found = true; continue;}
	      }
	    }
      if(found == false){cout << aa << " " << bb << " " << xx << " " << yy << " " << zz << endl; eventi_inpiu++;}
    }  
    cout << eventi_inpiu << " more events" << endl; 
}
