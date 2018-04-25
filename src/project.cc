//////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
  
#include <include/common.h>
#include <src/functions.h>
#include <src/histogram.h>

///////////////////////////////////////////////////////////////////////////////
//			MAIN
///////////////////////////////////////////////////////////////////////////////
int main(int argc,char *argv[]){
	int ind_arg;
	char *prm_name;

	time_t time_start,time_end;

	char filein[201];
	ifstream file_in;
	char fileout[201],fileout_enemat[201],fileout_mep[201];
	char fileout_DXene[201],fileout_DXhist[201];

	double kT,temperature;
	vector<double> str_pnt;
	vector<double> end_pnt;

	bool debug,flag,collect_samples;
	int prjdim,nummin;
	double distmin;
	char cmt;
	histogram HIST;

//---------------------------------------------------------
//		START-UP	BEGIN
//Floatint point exception
#ifdef HAVE_FEENABLEEXCEPT
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_OVERFLOW);
	feenableexcept(FE_INVALID);
#endif

//Initialize random number generator
	srand((unsigned)time(0));

//To exclude the path from the program name
	if(strrchr(argv[0],'/'))prm_name=(strrchr(argv[0],'/')+1);
	else prm_name=argv[0];

//Timer
	time(&time_start);

//Input-Output file names
	strcpy(filein,prm_name);
	strcat(filein,".in.dat");
	strcpy(fileout,"NULL");
	strcpy(fileout_enemat,"NULL");
	strcpy(fileout_mep,"NULL");
	strcpy(fileout_DXene,"NULL");
	strcpy(fileout_DXhist,"NULL");

//Initialization
	debug = false;
	temperature = 300;	//Temperature [k]
	prjdim = 1;	//Dimensionality of the projection
	nummin = 100; // Number of minima searched
	distmin = 2; // Distance between the minima
	cmt = '#'; // First character of comment lines # = gnuplot % = matlab
			// also the output format changes in accordance
	collect_samples = false; // Write also the number of samples in the projections
	str_pnt.clear(); // Start point for MEP
	end_pnt.clear(); // End point for MEP

//Parameters reading
	for(ind_arg=1;ind_arg<(argc);ind_arg++){
	        if(argv[ind_arg][0] == '-')
	        {
	      	  if(strcmp(argv[ind_arg]+1,"h")==0){//Help
	      		  printhelp(prm_name);
	      		  return 0;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"i")==0){//Input File
	      		  strcpy(filein,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"o")==0){//Output File
	      		  strcpy(fileout,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"enemat")==0){//Output energy in matlab format
	      		  strcpy(fileout_enemat,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"DXene")==0){//Output energy in DX format
	      		  strcpy(fileout_DXene,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"DXhist")==0){//Output histograms in DX format
	      		  strcpy(fileout_DXhist,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"mep")==0){//Output Minimum Energy Path
	      		  strcpy(fileout_mep,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"matlab")==0){//Use matlab format for mep ?
			  cmt = '%';
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"prjdim")==0){//Dimensionality of the projection
	      		  prjdim = atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"prjsamples")==0){//Project also the number of samples
			  collect_samples = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"nummin")==0){//Number of searched minima
	      		  nummin = atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"distmin")==0){//Distance between minima
	      		  distmin = atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"temperature")==0){//Temperature
	      		  temperature=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"mepstr")==0){//Start point for Minimum Energy Path
			  if ( (ind_arg + 1) == argc) flag = false;
			  else flag = true;
			  while (flag) {
			  	if((ind_arg + 1) == argc)flag = false;
			  	else {
					if (strspn(argv[ind_arg+1],".-0123456789") == strlen(argv[ind_arg+1])) {
	      		  			str_pnt.push_back(atof(argv[ind_arg+1]));
			  			ind_arg++;
					} else {
						flag = false;
					}
				}
			  }
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"mepend")==0){//End point for Minimum Energy Path
			  if ( (ind_arg + 1) == argc) flag = false;
			  else flag = true;
			  while (flag) {
			  	if((ind_arg + 1) == argc)flag = false;
			  	else {
					if (strspn(argv[ind_arg+1],".-0123456789") == strlen(argv[ind_arg+1])) {
	      		  			end_pnt.push_back(atof(argv[ind_arg + 1]));
			  			ind_arg++;
					} else {
						flag = false;
					}
				}
			  }
	      	  }
	      	  else{//Wrong Option
	      		  cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
			  for(ind_arg=0;ind_arg<(argc);ind_arg++)cerr<<argv[ind_arg]<<" ";
	      		  cerr<<endl;
	      		  printhelp(prm_name);
	      		  return 1;
	      	  }
	        }
	}

//Derived parameters
	kT = 1.9845e-3 * temperature;	//[Kcal/mol]
	
//Output control
	file_in.open(filein);
	if(!file_in){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE INPUT FILE"<<endl;
		return 1;
	}

	cout<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;
//		START-UP	END
//---------------------------------------------------------

  
//---------------------------------------------------------
//		PROJECTION	START
	HIST.readrestart(file_in);//Construct the HIST object reading the restart file
	if(strncmp(fileout,"NULL",4))HIST.WriteProjection(fileout,prjdim,collect_samples,kT);//Write projections in fileout
	if(strncmp(fileout_enemat,"NULL",4))HIST.WriteEneMat(fileout_enemat,kT);//Write the energy for matlab
	if(strncmp(fileout_DXhist,"NULL",4))HIST.WriteDxHist(fileout_DXhist);//Write the histogram in DX format
	if(strncmp(fileout_DXene,"NULL",4))HIST.WriteDxEne(fileout_DXene);//Write the energy in DX format
	if(strncmp(fileout_mep,"NULL",4))HIST.FindMEP(fileout_mep,str_pnt,end_pnt,cmt,kT,nummin,distmin);//Compute the MEP
//		PROJECTION END
//---------------------------------------------------------

	file_in.close();

	time(&time_end);
	cout<<endl<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;
  
	return 0;
}
