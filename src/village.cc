//////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2006 Simone Furini <simone.furini@unisi.it>
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

	char filelist[201];
	char fileout[201];
	char fileout2[201];
	char fileene[201];
	char filerestart[201];
	char filerestartout[201];
	char prefix[201];
	char stmp[11];
	ifstream file_in;
	ofstream file_out,file_out2;

	histogram HIST;
	vector<int> col;
	int first,last,stride,freqout,freqrmsd,numit,it,idim,itmp;
	int ind_reiter,num_reiter, num_split, num_bootstrap;
	double tol,rmsd,temperature,kT,timestep_traj;
	bool readfromlist,colFROMfile,extFROMfile,cor,converged,debug,trsFROMfile,massFROMfile;
	bool discout,histout,histprj,biasout;

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

//Initialization
	time(&time_start);
	debug = false;	//Debug output ?
	readfromlist = true;	//Restart file or read from list ?
	timestep_traj = 2e-3; // [ps] Timestep of the traj files
	tol = 1e-3;	//Tollerance of the WHAM procedure
	numit = 0;	//Maximum number of WHAM iterations
	freqout = 1000;	//Frequency for energy and restart output
	freqrmsd = 1000;//Frequency for rmsd checking
	temperature = 300;	//Temperature [k]
	first = 0;	//First time step to read from the data files
	last = 0;	//Last time step to read from the data files, 0 = Until the end
	stride = 1;	//Distance between uncorrelated samples in data files
				// If cor == false, read with this step
				//    cor == true, use this value as the g of the WHAM eqs
				// stride < 0 --> Estimate stride at run-time
				// 	cor needs to be true in this case
	cor = false;	//Use the WHAM equation for uncorrelate/correlated samples
	colFROMfile = false;	//Read the column organization from file.list
	extFROMfile = false;	//Read first/last/stride from file.list
	massFROMfile = false;	//Read masses from file.list
	trsFROMfile = false;	//Read coordinate translations from file.list
	num_split = 1; // Number of iteration for splitted trajectories
	num_bootstrap = 1; // Number of iteration to run ( 1 = Normal Wham)
	histout = false; // Write histogram output file ?
	histprj = false; //Write projections in Histrogram outout file ?
	biasout = false; // Write bias output file ?
	discout = false; // Write discretized coordinate ?
	strcpy(filerestart,"NULL"); // Input restart file
	strcpy(prefix,prm_name);
	strcpy(filelist,prefix);
	strcat(filelist,".list"); // List of the files with the trajectories

//Parameters reading
	if((argv[1][0] == '-')&&(strcmp(argv[1]+1,"h") == 0)){
		printhelp(prm_name);
		return 1;
	}
	if(argc < 3){
		cerr<<"ERROR IN "<<prm_name<<" CALL: parameters missing"<<endl;
		for(ind_arg=0;ind_arg<(argc);ind_arg++)cerr<<argv[ind_arg]<<" ";
		cerr<<endl;
		printhelp(prm_name);
		return 1;
	}
	if((argv[1][0] == '-')&&(strcmp(argv[1]+1,"restart") == 0)){//Read everything from restart file
	      	strcpy(filerestart,argv[2]);
		file_in.open(filerestart);
		if(!file_in){
			cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filerestart<<endl;
			return 1;
		}
		cerr<<"# Reading data from restart file "<<filerestart<<endl;
		file_out<<"# Reading data from restart file "<<filerestart<<endl;
		HIST.readrestart(file_in);//Construct the HIST object reading the restart file
		file_in.close();
		file_in.clear();
		ind_arg = 3; //Where to start reading the other parameters on the command line
		readfromlist = false;	//Do not read list file
	} else {
		if(argc < 6){
			cerr<<"ERROR IN "<<prm_name<<" CALL: parameters missing"<<endl;
			for(ind_arg=0;ind_arg<(argc);ind_arg++)cerr<<argv[ind_arg]<<" ";
			cerr<<endl;
			printhelp(prm_name);
			return 1;
		}
		if((argv[1][0] != '-')||(strcmp(argv[1]+1,"numdim") != 0)){//Dimensionality of the PMF
			cerr<<"ERROR IN "<<prm_name<<" CALL: first parameter must be -numdim"<<endl;
			for(ind_arg=0;ind_arg<(argc);ind_arg++)cerr<<argv[ind_arg]<<" ";
			cerr<<endl;
			printhelp(prm_name);
			return 1;
		}
		HIST.numdim=atoi(argv[2]);//Define number of dimension in the HIST object
		if (argc < (3+HIST.numdim*3)){
			cerr<<"ERROR IN "<<prm_name<<" CALL: parameters missing"<<endl;
			for(ind_arg=0;ind_arg<(argc);ind_arg++)cerr<<argv[ind_arg]<<" ";
			cerr<<endl;
			printhelp(prm_name);
			return 1;
		}
		HIST.AllocateDimensions();//Create the data structures in HIST
		itmp = 3;
		for(idim = 0;idim < HIST.numdim; idim++){
			HIST.DefineBoundaries(idim,atof(argv[itmp+(idim*3)])
					,atof(argv[itmp+(idim*3)+1]),atoi(argv[itmp+(idim*3)+2]));
			if (strcmp(argv[itmp+(idim*3)+3],"P")==0){ // Force periodicity on boundary points
				HIST.periodic[idim] = 1;
				//Number of points used to force P
				HIST.whr_prd[idim] = (double)atoi(argv[itmp+(idim*3)+4]);
				itmp = itmp + 2;
			} else if (strcmp(argv[itmp+(idim*3)+3],"T")==0){ // Force the period T in the PMF
				HIST.periodic[idim] = 2;
				HIST.whr_prd[idim] = atof(argv[itmp+(idim*3)+4]); //This is the period to force
				itmp = itmp + 2;
			} else HIST.periodic[idim] = 0;

		}
		HIST.InitializeDimensions();//Define the data structures in HIST 
		ind_arg = 3 + 3*HIST.numdim; // Where to start reading the other parameters on the command line
	}

	for(;ind_arg<(argc);ind_arg++){
	        if(argv[ind_arg][0] == '-')
	        {
	      	  if(strcmp(argv[ind_arg]+1,"h")==0){//Help
	      		  printhelp(prm_name);
	      		  return 0;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"list")==0){//File containg the list datafiles
	      		  strcpy(filelist,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"restart")==0){//Read a restart file ?
			  cerr<<"ERROR: restart file MUST be the first parameter"<<endl;
			  return 1;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"prefix")==0){//Prefix for the output files
	      		  strcpy(prefix,argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"histout")==0){//Histogram output File
			  histout = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"histprj")==0){//Write projections in Histogram output File
			  histprj = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"biasout")==0){//Biasing potentials output File
			  biasout = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"discout")==0){//Discretized Data output File
			  discout = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"cor")==0){//Use correlated WHAM eqs
			  cor = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"massFROMfile")==0){//Read masses from file.list
			  massFROMfile = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"trsFROMfile")==0){//Read coordinate translations from file.list
			  trsFROMfile = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"extFROMfile")==0){//Read first/last/stride from trj file
			  extFROMfile = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"tol")==0){//Tollerance
	      		  tol=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"numit")==0){//Maximum number of iterations
	      		  numit=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"freqout")==0){//Frequency for energy and restart output
	      		  freqout=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"freqrmsd")==0){//Frequency for rmsd checking
	      		  freqrmsd=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"first")==0){//First time step to read from the data files
	      		  first=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"last")==0){//last time step to read from the data files, 0 = end of file
	      		  last=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"stride")==0){//distance between steps read from the data files
	      		  stride=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"debug")==0){//Debug output ?
			  debug = true;
			  HIST.debug = true;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"temperature")==0){//Temperature
	      		  temperature=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"dt")==0){//Timestep of the traj files
	      		  timestep_traj=atof(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"bootstrap")==0){//Run bootstrap analyses ? How many ?
	      		  num_bootstrap=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"split")==0){//Split the trajectories ? Num of bunches and period
	      		  num_split=atoi(argv[ind_arg+1]);
			  ind_arg++;
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"mass")==0){//Mass of the restrained particles
			  if((ind_arg + HIST.numdim) < argc){
			  	for(idim = 0; idim < HIST.numdim; idim++){
			  		if((atof(argv[ind_arg+1+idim])) <= 0){
			        		cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
			        		printhelp(prm_name);
			        		return 1;
			        	}
	      		  		HIST.mass.push_back(atof(argv[ind_arg+1+idim]));
			  	}
			  }else{
			        cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
			        printhelp(prm_name);
			        return 1;
			  }
			  ind_arg += (HIST.numdim);
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"reference")==0){//Region where the average value of A is set to 0
			  if((ind_arg + 2*HIST.numdim) < argc){
			  	for(idim = 0; idim < HIST.numdim; idim++){
					HIST.DefineReferences(idim,atof(argv[ind_arg+(idim*2)+1])
							,atof(argv[ind_arg+(idim*2)+2]));
			  	}
			  }else{
			        cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
			        printhelp(prm_name);
			        return 1;
			  }
			  ind_arg += 2*(HIST.numdim);
	      	  }
	      	  else if(strcmp(argv[ind_arg]+1,"col")==0){//Columns to read in filedata
			  if(strcmp(argv[ind_arg+1],"listfile")==0){
				//Column organization is different in each file
				//and given as the last values in list file
				colFROMfile = true;
				ind_arg++;
			  } else {
			  	if((ind_arg + HIST.numdim) < argc){//Same column organization in all the files
			  		for(idim = 0; idim < HIST.numdim; idim++){
			  			if((atoi(argv[ind_arg+1+idim])) < 1){
			  	      		cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
			  	      		printhelp(prm_name);
			  	      		return 1;
			  	      	}
	      		  			col.push_back(atoi(argv[ind_arg+1+idim]));
			  		}
			  	}else{
			  	      cerr<<"ERROR IN "<<prm_name<<" CALL"<<endl;
			  	      printhelp(prm_name);
			  	      return 1;
			  	}
			  	ind_arg += (HIST.numdim);
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

	//Derive&Check parameters
	kT = 1.9845e-3 * temperature;	//[kcal/mol]
	if(col.size() == 0){//Col not specified in command line
		for(idim = 0; idim < HIST.numdim; idim++){
			col.push_back(idim+1);
		}
	}
	strcpy(fileout,prefix);
	strcat(fileout,".out"); // Log output of the program
	strcpy(fileene,prefix);
	strcat(fileene,".ene"); // Energy - Probability - Number of samples
	strcpy(filerestartout,prefix); // Output restart file
	strcat(filerestartout,".restart.dat");
	if(HIST.mass.size() == 0){//Masses not specified in command line
		for(idim = 0; idim < HIST.numdim; idim++){
			HIST.mass.push_back(1.0);
		}
	}
	if((num_bootstrap != 1) && (!readfromlist)){
		cerr<<"ERROR: bootstrap analysis does not work with restart file"<<endl;
		return 1;
	}
	if((num_bootstrap != 1) && (num_split != 1)){
		cerr<<"ERROR: choose between bootstrap or splitting !"<<endl;
		return 1;
	}
	if((stride < 0) && (!cor)){
		cerr<<"ERROR: Use correlated wham equations with run-time estimate of autocorrelation time"<<endl;
		return 1;
	}
	num_reiter = (num_bootstrap > num_split)? num_bootstrap : num_split;
	//Output control
  	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
		return 1;
	}
//		START-UP	END
//---------------------------------------------------------

//---------------------------------------------------------
//		WHAM	START
	cout<<"# "<<prm_name<<" starts at: "<<ctime(&time_start)<<endl;

	for(ind_reiter = 0; ind_reiter < num_reiter; ind_reiter++){
		if(num_reiter>1){
			cout<<"# Running PMF calculation "<<ind_reiter+1<<" of "<<num_reiter<<endl;
			file_out<<"# Running PMF calculation "<<ind_reiter+1<<" of "<<num_reiter<<endl;
		}
		//Read trajectories
		if(readfromlist){
  			file_in.open(filelist);
			if(!file_in){
				cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filelist<<endl;
				return 1;
			}
			HIST.readlist(file_in,file_out,col,first,last,cor,stride
					,trsFROMfile,massFROMfile,colFROMfile,extFROMfile
					,ind_reiter,num_split,num_bootstrap
					,kT,discout,prefix,timestep_traj);//Read from list file
			file_in.close();
			file_in.clear();
		}

		//Prepare data structures
		HIST.InitializeWham(file_out);

		//Iterative procedure
		it = 0;
		converged = false;
		rmsd = tol*1e3;	//To be out of tollerance at the beggining
		while((it < numit)&&(!converged)){
			HIST.NewProbabilities(file_out,it);//Compute new probabilities (Punorm)
			HIST.ForcePeriodicity(file_out,it);//If periodic coordinates force the periodicity (Punorm)
			HIST.NewConstants(file_out,it);//Compute new biasing constants, it uses Punorm
			HIST.NormalizeProbability(); // Punorm --> Pnorm, it changes F
			if((it!=0) && ((it%freqrmsd)==0)){//Check convergence
				rmsd = HIST.CheckConvergence(file_out,it,kT);
				cerr<<"# Iteration "<<it<<" RMSD = "<<rmsd<<" on "<<HIST.numgood<<" points"<<endl;
				file_out<<"# Iteration "<<it<<" RMSD = "<<rmsd<<" on "<<HIST.numgood<<" points"<<endl;
			}
			if((it!=0) && ((it%freqout)==0)){//Output
				file_out2.open(filerestartout);
				if(!file_out2){
					cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filerestartout<<endl;
					return 1;
				}
				HIST.WriteRestart(file_out2);
				file_out2.close();
				file_out2.clear();
			}
			if (rmsd < tol)converged = true;//Is it converged ?
			it++;
		} // END while((it < numit)&&(!converged)){

		// Memorize Energy
		HIST.MemorizeEne(kT);

		// Output Energy
		HIST.WriteEne(fileene,ind_reiter,kT);

		// Output Histograms
		if(histout){
			strcpy(fileout2,prefix);
			if(num_reiter != 1){
				sprintf(stmp,".%i",ind_reiter);
				strcat(fileout2,stmp);
			}
			strcat(fileout2,".hist");
			file_out2.open(fileout2);
			if(!file_out2){
				cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout2<<endl;
				return 1;
			}
			HIST.WriteHist(file_out2,histprj,kT);
			file_out2.close();
			file_out2.clear();
		}
	
		// Output Biasing potentials
		if(biasout){
			strcpy(fileout2,prefix);
			if(num_reiter != 1){
				sprintf(stmp,".%i",ind_reiter);
				strcat(fileout2,stmp);
			}
			strcat(fileout2,".bias");
			file_out2.open(fileout2);
			if(!file_out2){
				cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout2<<endl;
				return 1;
			}
			HIST.WriteBias(file_out2,kT);
			file_out2.close();
			file_out2.clear();
		}
	} // END for(ind_reiter = 0; ind_reiter < num_reiter; ind_reiter++){
//		WHAM	END
//---------------------------------------------------------

	// Write Energy, computed by reiterations
	file_out2.open(fileene,fstream::app);
	if(!file_out2){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileene<<endl;
		return 1;
	}
	HIST.WriteEne(file_out2,num_reiter);
	file_out2.close();

	file_out.close();
	time(&time_end);
	cout<<"# "<<prm_name<<" ends at: "<<ctime(&time_start)<<endl;
  
	return 0;
}
