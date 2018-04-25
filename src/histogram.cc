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

////////////////////////////////////////////////////////////
//			CONSTRUCTORS
histogram::histogram(){
	debug = false;
	numdim = numwin = numhist = numgood = numbad = 0;
	histmin = g = histmax = P = A = NULL;
	delta = harmrest = center = numwham = denwham = Punnorm = Pold = Aold = F = sumP = whr_prd = NULL;
	numbin = step = numbinwin = hist = periodic = NULL;
	flaginf = NULL;
	mass.clear();
	backup_energy_reiter.clear();
	return;
}
//			END CONSTRUCTORS
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
//			INPUTS
int histogram::readrestart(ifstream &file_in){
	file_in.read((char *)&numdim,sizeof(int));
	file_in.read((char *)&numwin,sizeof(int));
	AllocateDimensions();
	file_in.read((char *)histmin,numdim*sizeof(double));
	if((unsigned long int)file_in.gcount() != numdim*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)histmax,numdim*sizeof(double));
	if((unsigned long int)file_in.gcount() != numdim*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)numbin,numdim*sizeof(int));
	if((unsigned long int)file_in.gcount() != numdim*sizeof(int)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)periodic,numdim*sizeof(int));
	if((unsigned long int)file_in.gcount() != numdim*sizeof(int)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)whr_prd,numdim*sizeof(double));
	if((unsigned long int)file_in.gcount() != numdim*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	InitializeDimensions();
	AllocateWindows();
	file_in.read((char *)center,numdim*numwin*sizeof(double));
	if((unsigned long int)file_in.gcount() != numdim*numwin*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)harmrest,numdim*numwin*sizeof(double));
	if((unsigned long int)file_in.gcount() != numdim*numwin*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)numbinwin,numwin*sizeof(int));
	if((unsigned long int)file_in.gcount() != numwin*sizeof(int)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)hist,numhist*numwin*sizeof(int));
	if((unsigned long int)file_in.gcount() != numhist*numwin*sizeof(int)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)P,numhist*sizeof(double));
	if((unsigned long int)file_in.gcount() != numhist*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	file_in.read((char *)F,numwin*sizeof(double));
	if((unsigned long int)file_in.gcount() != numwin*sizeof(double)){
		cerr<<"ERROR: reading from a restart file of wrong dimension"<<endl;
		exit(1);
	}
	return 0;
}

int histogram::AllocateDimensions(){
	if (numdim <=0 ){
		cerr<<"ERROR number of dimensions cannot be "<<numdim<<endl;
		exit(1);
	}
	histmin = new double [numdim];
	histmax = new double [numdim];
	refmin = new double [numdim];
	refmax = new double [numdim];
	numbin = new int [numdim];
	periodic = new int [numdim];
	delta = new double [numdim];
	step = new int [numdim];
	whr_prd = new double [numdim];
	if ((!histmin)||(!histmax)||(!numbin)||(!periodic)||(!delta)||(!step)||(!whr_prd)){
		cerr<<"ERROR: memory allocation failed"<<endl;
		exit(1);
	}
	return 0;
}

int histogram::InitializeDimensions(){
	int idim;
	for(idim = 0;idim < numdim; idim++){
		delta[idim] = (histmax[idim] - histmin[idim])/(numbin[idim]-1);
		if(idim == 0) step[idim] = 1;
		else step[idim] = step[idim-1]*numbin[idim-1];
		if(histmax[idim] <= histmin[idim]){
			cerr<<"ERROR: histogram end lower than histogram begin"<<endl;
			exit(1);
		}
		if(numbin[idim] <= 0){
			cerr<<"ERROR: negative number of bins"<<endl;
			exit(1);
		}
	}
	return 0;
}

int histogram::AllocateWindows(){
	int idim, iwin, ihist;
	// Delete previous allocations, if any
	if(numbinwin) delete[] numbinwin;
	if(g) delete[] g;
	if(F) delete[] F;
	if(sumP) delete[] sumP;
	if(center) delete[] center;
	if(harmrest) delete[] harmrest;
	if(A) delete[] A;
	if(Aold) delete[] Aold;
	if(P) delete[] P;
	if(Pold) delete[] Pold;
	if(Punnorm) delete[] Punnorm;
	if(numwham) delete[] numwham;
	if(denwham) delete[] denwham;
	if(flaginf) delete[] flaginf;
	if(hist) delete[] hist;
	// Count the number of grid elements
	numhist = numbin[0];
	for(idim = 1; idim < numdim; idim++) numhist = numhist*numbin[idim];
	// Allocate !
	numbinwin = new int [numwin];
	g = new double [numwin];
	F = new double [numwin];
	sumP = new double [numwin];
	center = new double [numdim*numwin];
	harmrest = new double [numdim*numwin];
	A = new double [numhist];
	Aold = new double [numhist];
	P = new double [numhist];
	Pold = new double [numhist];
	Punnorm = new double [numhist];
	numwham = new double [numhist];
	denwham = new double [numhist];
	flaginf = new bool [numhist];
	hist = new int [numhist*numwin];
	// Check allocation
	if ( (!numbinwin)||(!g)||(!F)||(!sumP)
		||(!center)||(!harmrest)
		||(!A)||(!Aold)||(!P)||(!Punnorm)||(!Pold)||(!numwham)||(!denwham)||(!flaginf)
		||(!hist) ){
		cerr<<"ERROR: memory allocation failed"<<endl;
		exit(1);
	}
	// Reset
	numgood = numbad = 0;
	for(iwin = 0; iwin < numwin; iwin++){
		numbinwin[iwin] = 0;
		g[iwin] = 1.0;
		F[iwin] = sumP[iwin] = 0.0;
		for(idim = 0; idim < numdim; idim++){
			center[idim+numdim*iwin] = harmrest[idim+numdim*iwin] = 0.0;
		}
		for(ihist = 0; ihist < numhist; ihist++){
			hist[ihist+iwin*numhist] = 0;
		}
	}
	for(ihist = 0; ihist < numhist; ihist++){
		A[ihist] = Aold[ihist] = P[ihist] = Punnorm[ihist] = Pold[ihist] = numwham[ihist] = denwham[ihist] = 0.0;
		flaginf[ihist] = false;
	}
	return 0;
}

int histogram::DefineBoundaries(int idim, double min, double max, int num){
	histmin[idim] = min;
	histmax[idim] = max;
	numbin[idim] = num;
	return 0;
}

int histogram::DefineReferences(int idim, double min, double max){
	refmin[idim] = min;
	refmax[idim] = max;
	return 0;
}

int histogram::readlist(ifstream &file_in,ofstream &file_out
		, vector<int> col, int first, int last, bool cor, int stride
		, bool trsFROMfile, bool massFROMfile, bool colFROMfile, bool extFROMfile
		, int ind_reiter, int num_split, int num_bootstrap
		, double kT, bool discout, char *prefix, double timestep_traj){
	char line[401],stmp[201],stmp2[201];
	int ipos,ibin,idim,jdim,iwin,icol,itmp,coin,numcol,numread,split_period,first_2read,last_2read, iread_old;
	double pos,trsl[numdim];
	double ftmp,fread_old;
	ofstream file_out2;
	ifstream file_data;
	string srtmp,filedata;
	vector <long int> grid_index;
	vector<double> vtmp;

	// Count the number of trajectory files = windows
	file_in.getline(line,400);
	numwin = 0;
	while(!file_in.eof()){
		if(line[0]!='#')numwin++;
		file_in.getline(line,400);
	}
	file_out<<"# Reading data from "<<numwin<<" files"<<endl;
	file_in.clear();
	file_in.seekg(ios_base::beg);
	// Allocate the memory for the number of files to read
	AllocateWindows();
	// Expand masses definition to all the windows
	if(massFROMfile){
		mass.clear();
	} else {
		vtmp = mass;
		mass.clear();
		for(iwin=0;iwin<numwin;iwin++){
			for(idim=0;idim<numdim;idim++){
				mass.push_back(vtmp[idim]);
			}
		}
		vtmp.clear();
	}
	// Save the stride in g[iwin] if this is not read from the file.list
	if(!extFROMfile){
		for(iwin=0;iwin<numwin;iwin++)g[iwin]=stride;
	}
	//Output discretized coordinates - initialization
	if(discout){
		strcpy(stmp,prefix);
		if((num_split > 1) || (num_bootstrap != 1)){
			sprintf(stmp2,"%i",ind_reiter);
			strcat(stmp,stmp2);
		}
		file_out2.open(stmp);
		if(!file_out2){
			cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<stmp<<endl;
			return 1;
		}
		WriteHeading(file_out2);
	}
	//Reading data
	iwin = 0;
	file_in.getline(line,400);
	srtmp = line;
	while(!file_in.eof()){
		while(srtmp[0]==' ')srtmp.erase(0,1);
		if((srtmp[0]!='#')&&(srtmp.length()>0)){
			//Open file
			filedata = srtmp.substr(0,srtmp.find_first_of(" \t"));
			file_data.open(filedata.c_str());
			if(!file_data){
				cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<filedata<<endl;
				exit(1);
			}
			//Read translations
			if(trsFROMfile){
				for(idim = 0;idim < numdim; idim++){
					srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
					srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
					if(srtmp.length()==0){
						cerr<<"ERROR reading file "<<endl;
						cerr<<line<<endl;
						return 1;
					}
					trsl[idim] = atof(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
				}
			} else {
				for(idim = 0;idim < numdim; idim++)trsl[idim] = 0.0;
			}
			//Read centers
			for(idim = 0;idim < numdim; idim++){
				srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
				srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
				if(srtmp.length()==0){
					cerr<<"ERROR reading file "<<endl;
					cerr<<line<<endl;
					return 1;
				}
				center[idim+iwin*numdim] = atof(
						srtmp.substr(0,srtmp.find_first_of(" \t")).c_str()) - trsl[idim];
				if((center[idim+iwin*numdim]<histmin[idim])||
					(center[idim+iwin*numdim]>histmax[idim])){
					if(debug){
						cerr<<"WARNING: window center out of the analyzed region"<<endl;
						cerr<<" min = "<<histmin[idim]<<" max = "<<histmax[idim]
							<<" center = "<<center[idim+iwin*numdim]<<endl;
					}
				}
			}
			//Read restraints
			for(idim = 0;idim < numdim; idim++){
				srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
				srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
				if(srtmp.length()==0){
					cerr<<"ERROR reading file "<<endl;
					cerr<<line<<endl;
					exit(1);
				}
				harmrest[idim+iwin*numdim] = atof(
						srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
				//Convert restraints units to [kT]
				harmrest[idim+iwin*numdim] = harmrest[idim+iwin*numdim]/kT;
			}
			//Read column organization
			if(colFROMfile){
				col.erase(col.begin(),col.end());
				for(idim = 0; idim < numdim; idim++){
					srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
					srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
					if(srtmp.length()==0){
						cerr<<"ERROR reading file "<<endl;
						cerr<<line<<endl;
						return 1;
					}
	      			  	col.push_back(atoi(
						srtmp.substr(0,srtmp.find_first_of(" \t")).c_str()));
				}
			}
			//Read masses
			if(massFROMfile){
				for(idim = 0;idim < numdim; idim++){
					srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
					srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
					if(srtmp.length()==0){
						cerr<<"ERROR reading file "<<endl;
						cerr<<line<<endl;
						exit(1);
					}
	      			  	mass.push_back(atof(
						srtmp.substr(0,srtmp.find_first_of(" \t")).c_str()));
				}
			}
			//Read first/last/seed
			if(extFROMfile){
				srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
				srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
				if(srtmp.length()==0){
					cerr<<"ERROR reading file "<<endl;
					cerr<<line<<endl;
					return 1;
				}
				first = atoi(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
				srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
				srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
				if(srtmp.length()==0){
					cerr<<"ERROR reading file "<<endl;
					cerr<<line<<endl;
					return 1;
				}
				last = atoi(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
				srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
				srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
				if(srtmp.length()==0){
					cerr<<"ERROR reading file "<<endl;
					cerr<<line<<endl;
					return 1;
				}
				g[iwin] = atof(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
			}
			// Count the lines to read
			file_data.getline(line,400);
			srtmp = line;
			split_period = 0;
			//iread_old = -1;
			fread_old = -timestep_traj;
			while(!file_data.eof()){
				while(srtmp[0]==' ')srtmp.erase(0,1);
				if((srtmp[0]!='#')&&(srtmp.length()>0)){
					// Read Time Step
					//itmp = atoi(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
					ftmp = atof(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
					//Use just the samples between first and last (if this is defined)
					//if( (itmp != iread_old) && (itmp >= first) && ((last == 0)||(itmp < last)) ){
					if( (abs(ftmp-fread_old) > 0.5*timestep_traj) && (ftmp >= first) && ((last == 0)||(ftmp < last)) ){
						split_period++;
						//iread_old = itmp;
						fread_old = ftmp;
					}
				}
				file_data.getline(line,400);
				srtmp = line;
			}
			file_data.clear();
			file_data.seekg(ios_base::beg);
			// Select bunch, if split analysis requested
			if (num_split > 1) {
				cout<<"# Reading "<<split_period/num_split<<" samples out of "<<split_period<<" from "<<filedata;
				//split_period = split_period / num_split;
				//first_2read = first + split_period * ind_reiter;
				//last_2read = first_2read + split_period;
				split_period = fread_old / num_split;
				first_2read = first + split_period * ind_reiter;
				last_2read = first_2read + split_period;
			} else { // otherwise read from first to last
				cout<<"# Reading "<<split_period<<" samples from "<<filedata;
				first_2read = first;
				last_2read = last;
			}
			if(cor) cout<<", to read ["<<first_2read<<":"<<last_2read<<"] stride = "<<timestep_traj<<" steps";
			else cout<<", to read ["<<first_2read<<":"<<last_2read<<"] stride = "<<g[iwin]<<" steps";
			//Read file data
			file_data.getline(line,400);
			srtmp = line;
			grid_index.clear(); // This is used to memorize the positions of the points
						// read from the file in the grid
						// = -1 --> Outside the grid
						// >= 0 --> Index on the grid
			if(cor){ // if cor read all the frames and keep g[iwin] for the analyses
				stride = 1;
			} else { // if cor==false read only every g[iwin] frames and reset g[iwin] to 1
				stride = g[iwin];
				g[iwin] = 1.0;
			}
			iread_old = 0;
			fread_old = -timestep_traj;
			while(!file_data.eof()){
				while(srtmp[0]==' ')srtmp.erase(0,1);
				if((srtmp[0]!='#')&&(srtmp.length()>0)){
					// Read Time Step
					//itmp = atoi(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
					ftmp = atof(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
					//Use just the samples between first and last (if this is defined)
					//read only timesteps multiple of stride if cor == false
					//read all timesteps if cor == true
					//if( (itmp != iread_old) &&  (itmp >= first_2read) && ((last_2read == 0)||(itmp < last_2read)) && ((itmp % stride) == 0) ){
					if( (abs(ftmp-fread_old) > 0.5*timestep_traj) &&  (ftmp >= first_2read) && ((last_2read == 0)||(ftmp < last_2read)) && ((iread_old % stride) == 0) ){
						//iread_old = itmp;
						fread_old = ftmp;
						iread_old += 1;
						//Delete time step from srtmp
						srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
						//Delete time step from srtmp
						srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
						//Read Positions
						grid_index.push_back(0); // index of the sample in the grid
						numcol = 1;//Number of read column
						numread = 0;//Number of read data
						while((numread < numdim)&&(srtmp.length()!=0)
								&&(grid_index[grid_index.size()-1]>=0)){
						  icol = (find(col.begin(),col.end(),numcol) - col.begin());
						  // Element is in a column to read ? 
						  if(  find(col.begin(),col.end(),numcol) != col.end()){
						    pos = atof(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str())
							    - trsl[icol];
						    ibin = (int)((pos-histmin[icol]+delta[icol]/2.0)/delta[icol]);
						    if((ibin<0)||(ibin>=numbin[icol])){
						      if(debug)file_out<<"# WARNING data out of boundaries"<<endl;
						      // Element outside the grid
						      grid_index[grid_index.size()-1] = -1;
						    }else{
						      // Element inside the grid in this dimension
						      grid_index[grid_index.size()-1] += ibin*step[icol];
						      numread++;
						    }
						  }
						  numcol++;
						  //Delete column
						  srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
						  srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
						}
						// Element inside the grid in all dimensions
						if((numread != numdim)
							&& ( (srtmp.length()==0)&&(numread!=(numdim-1)) ) ){
								cerr<<"ERROR reading "<<filedata<<endl;
								return 1;
						}
					}
				}
				file_data.getline(line,400);
				srtmp = line;
			}//while(!file_data.eof())
			if(g[iwin] < 0){ // Calculate the autocorrelation time
				double mean, norm, jres, kres, tau;
				double ACR[(int)grid_index.size()-1];
				int n_acr;
				bool stop;
				for(jdim = 0; jdim < numdim ; jdim++){
					// Compute mean
					mean = 0.0;
					for(int i=0; i < (int)grid_index.size(); i++){
						ipos = ((int)(grid_index[i]/step[jdim])) % numbin[jdim];
						mean += (histmin[jdim] + ipos*delta[jdim]);
					}
					mean = mean / (int)grid_index.size();
					// Compute autocorrelation
					n_acr = 0;
					stop = false;
					while ((n_acr < (int)grid_index.size()-1) && (!stop) ) {
						ACR[n_acr] = 0.0;
						for(int j=0; j < (int)grid_index.size()-n_acr; j++){
							jres = histmin[jdim] + (((int)(grid_index[j]/step[jdim])) % numbin[jdim])*delta[jdim] - mean;
							kres = histmin[jdim] + (((int)(grid_index[j+n_acr]/step[jdim])) % numbin[jdim])*delta[jdim] - mean;
							ACR[n_acr] += jres*kres;
						}
						ACR[n_acr] = ACR[n_acr] / ((int)grid_index.size() - n_acr);
						if (ACR[n_acr] < 0.0) stop = true; // Stop if negative autocorrelation
						n_acr++;
					}
					// Normalize + Integrate
					tau = 0.0;
					norm = ACR[0];
					for(int i=0; i < n_acr; i++){
						ACR[i] = ACR[i] / norm;
						tau +=  ACR[i]*(1.0 - i/((int)grid_index.size()));
						cout<<"ACR\t"<<i<<"\t"<<ACR[i]<<"\t"<<tau<<"\t"<<((int)grid_index.size())<<endl;
					}
				}
				g[iwin] = (int)tau;
			}
			if(RAND_MAX <= grid_index.size()){
				cerr<<"ERROR: maximum random number too low for the current problem"<<endl;
				exit(1);
			}
			cout<<", total samples = ["<<numgood;
			for(itmp = 0; itmp < (int)grid_index.size(); itmp++){
				if(num_bootstrap > 1)coin = rand() % grid_index.size();
				else coin = itmp;
				// Element inside the grid in all dimensions ?
				if(grid_index[coin] >= 0){
					numgood++;
					hist[grid_index[coin]+iwin*numhist] += 1;
					numbinwin[iwin]++;
					if(discout){
						file_out2<<line<<'\t';
						for(jdim = 0; jdim < numdim ; jdim++){
							ipos = ((int)(grid_index[coin]/step[jdim])) 
								% numbin[jdim];
							pos = histmin[jdim] + ipos*delta[jdim];
							file_out2<<pos<<'\t';
						}
						file_out<<endl;
					}
				}else{
					numbad++;
					if(discout){
						file_out2<<'#'<<line<<endl;
					}
				}
			}
			cout<<":"<<numgood<<"]"<<endl;
			if(discout)file_out2<<endl<<endl;
			//Output file parameters
			file_out<<"# Reading file "<<filedata<<" ["<<first_2read<<":"<<stride<<":"<<last_2read<<"] tau = "<<g[iwin]<<endl;
			if(trsFROMfile){
				file_out<<"#      Translations ";
				for(idim = 0;idim < numdim; idim++){
					file_out<<setw(10)<<trsl[idim];
				}
				file_out<<endl;
			}
			file_out<<"#      Centers ";
			for(idim = 0;idim < numdim; idim++){
				file_out<<setw(10)<<center[idim+iwin*numdim];
			}
			file_out<<endl<<"#   Restraints ";
			for(idim = 0;idim < numdim; idim++){
				file_out<<setw(10)<<harmrest[idim+iwin*numdim];
			}
			file_out<<endl<<"#   Masses ";
			for(idim = 0;idim < numdim; idim++){
				file_out<<setw(10)<<mass[idim+iwin*numdim];
			}
			file_out<<endl;
			//Close file data
			file_data.close();
			file_data.clear();
			iwin++;
		}
		file_in.getline(line,400);
		srtmp = line;
	}//end while(file_in)
	cerr<<"# "<<numbad<<" points outside the analyzed region ("<<numgood<<" inside)"<<endl;
	file_out<<"# "<<numbad<<" points outside the analyzed region ("<<numgood<<" inside)"<<endl;
	//Discretized output
	if(discout)file_out2.close();
	return numbad;
}
//			END INPUTS
////////////////////////////////////////////////////////////
  
////////////////////////////////////////////////////////////
//			WHAM
double histogram::bias(int ihist, int iwin){
        // To use restraints on the center of two coordinates put the 
	// first Force Constant to the actual value, and the second to zero
	// If the restaints is on the first and the last coordinate put the first 
	// Force Constant to zero
	double pos = 0.0;
	double U = 0.0;
	int jdim = 1, idim, ipos;
	if(numdim == 1)jdim = 0; //This is necessary to avoid NULL pointing
	for(idim = 0; idim < numdim ; idim++){
		ipos = ((int)(ihist/step[idim])) % numbin[idim];
		pos = histmin[idim] + ipos*delta[idim];
		if((harmrest[idim+iwin*numdim]==0.0)||(harmrest[jdim+iwin*numdim]==0.0)){
			ipos = ((int)(ihist/step[jdim])) % numbin[jdim];
			pos = (pos*mass[idim+iwin*numdim] + (histmin[jdim] + ipos*delta[jdim])*mass[jdim+iwin*numdim])/(mass[idim+iwin*numdim] + mass[jdim+iwin*numdim]);
		}
		U += 0.5 * harmrest[idim+iwin*numdim] *
			(pos - center[idim+iwin*numdim])*(pos - center[idim+iwin*numdim]);
		jdim++;
		if(jdim==numdim)jdim=0;
	}
	return U;
}

int histogram::InitializeWham(ofstream &file_out){
	int ihist, iwin;
	ihist = 0;
	numbad = numgood = 0;
	WriteHeading(file_out);
	if(debug)file_out<<"WHAM NUMERATORS"<<endl;
	for(ihist = 0; ihist < numhist; ihist++){ 
		for(iwin = 0; iwin<numwin; iwin++){
			numwham[ihist] += (hist[ihist+iwin*numhist]/g[iwin]);
		}
		if(numwham[ihist] < 1){
			flaginf[ihist] = true;
			numbad++;
			if(debug)file_out<<"# WARNING: probabilities too low in bin "<<ihist<<endl;
		}else{
			flaginf[ihist] = false;
			numgood++;
		}
		if(debug)file_out<<ihist<<'\t'<<numwham[ihist]<<endl;
	}
	cerr<<"# WARNING "<<numbad<<" bins with extremely low probability ("
		<<numgood<<" with sufficient probability)"<<endl;
	file_out<<"# WARNING "<<numbad<<" bins with extremely low probability ("
		<<numgood<<" with sufficient probability)"<<endl;
	return numbad;
}

int histogram::NewProbabilities(ofstream &file_out, int it){
	int ihist,iwin,jdim,ipos;
	double U;
	if(debug)file_out<<"WHAM PROBABILITIES IT = "<<it<<endl;
	for(ihist = 0; ihist < numhist; ihist++){ 
		if(!flaginf[ihist]){
		     denwham[ihist] = 0.0;
		     for(iwin = 0; iwin<numwin; iwin++){
		     	U = bias(ihist,iwin);
		     	denwham[ihist] += (numbinwin[iwin]/g[iwin]) * exp(F[iwin] - U);
		     }
		     if(debug){
			     file_out<<ihist<<'\t'<<numwham[ihist]<<'\t'<<denwham[ihist]<<'\t'<<Punnorm[ihist]<<'\t';
			     for(jdim = 0;jdim < numdim; jdim++){
			     	ipos = ((int)(ihist/step[jdim])) % numbin[jdim];
			     	file_out<<histmin[jdim] + ipos*delta[jdim]<<'\t';
			     }
			     file_out<<endl;
		     }
		     Punnorm[ihist] = numwham[ihist]/denwham[ihist];
		}
	}
	return 0;
}

int histogram::NewConstants(ofstream &file_out,int it){
	int iwin,ihist;
	double U;
	//Compute new F
	if(debug)file_out<<"WHAM CONSTANTS IT = "<<it<<endl;
	for(iwin = 0; iwin<numwin; iwin++){
		sumP[iwin] = 0.0;
		for(ihist = 0; ihist < numhist; ihist++){ 
			if(!flaginf[ihist]){
			     U = bias(ihist,iwin);
			     sumP[iwin] += (Punnorm[ihist] * exp(-U));
			}
		}
		F[iwin] = -log(sumP[iwin]);
		if(debug)file_out<<iwin<<'\t'<<sumP[iwin]<<'\t'<<F[iwin]<<endl;
	}
	//for(iwin = 1; iwin<numwin; iwin++){
	//	F[iwin] = F[iwin] - F[0];
	//}
	return 0;
}

int histogram::ForcePeriodicity(ofstream &file_out,int it){
	int idim,jdim,ipos[numdim],ihist,ihist_bnd,itmp1,itmp2;
	double value = 0.0, pos;
	vector<int> ihist_prd;

	for(idim = 0; idim<numdim; idim++){
		if(periodic[idim] == 1){
			for(ihist = 0; ihist < numhist; ihist++){ 
				for(jdim = 0; jdim<numdim; jdim++){
					ipos[jdim] = ((int)(ihist/step[jdim])) % numbin[jdim];
				}
				if(ipos[idim] < whr_prd[idim]){
					ihist_bnd = 0;
					for(jdim = 0; jdim<numdim; jdim++){
						if(jdim != idim)ihist_bnd += ipos[jdim] * step[jdim];
						else ihist_bnd += (numbin[jdim]-1-ipos[jdim]) * step[jdim];
					}
					value = (Punnorm[ihist] + Punnorm[ihist_bnd])/2.0;
					if(debug){
						itmp1 = ((int)(ihist/step[idim])) % numbin[idim];
						itmp2 = ((int)(ihist_bnd/step[idim])) % numbin[idim];
						file_out<<itmp1<<'\t'<<itmp2<<'\t'<<ihist<<'\t'<<ihist_bnd<<'\t'
							<<Punnorm[ihist]<<'\t'<<Punnorm[ihist_bnd]<<'\t'<<value<<endl;
					}
					Punnorm[ihist] = Punnorm[ihist_bnd] = value;
				}
			}
		} else if(periodic[idim] == 2){
			for(ihist = 0; ihist < numhist; ihist++){ 
				for(jdim = 0; jdim<numdim; jdim++){
					ipos[jdim] = ((int)(ihist/step[jdim])) % numbin[jdim];
				}
				pos = histmin[idim] + ipos[idim]*delta[idim];
				if((pos - histmin[idim]) < whr_prd[idim]){
					//apply only to the first period (or you'd apply n times)
					ihist_prd.clear();
					ihist_prd.push_back(ihist);
					value = Punnorm[ihist];
		//cout<<pos<<'\t';
					pos += whr_prd[idim];
					while(pos < histmax[idim]){
						ihist_bnd = 0;
						for(jdim = 0; jdim<numdim; jdim++){
							if(jdim != idim)ihist_bnd += ipos[jdim] * step[jdim];
							else ihist_bnd += (ipos[jdim] 
								+ (int)((whr_prd[jdim]*ihist_prd.size())/delta[jdim]) ) 
								* step[jdim];
						}
						ihist_prd.push_back(ihist_bnd);
						value += Punnorm[ihist_bnd];
		//cout<<pos<<'\t';
		//cout<<" T = "<<histmin[idim] + delta[idim]*(((int)(ihist_bnd/step[idim])) % numbin[idim])<<'\t';
						pos += whr_prd[idim];
					}
					value = value / ihist_prd.size();
		//cout<<"V = "<<value<<" # = "<<ihist_prd.size()<<endl;
					for(itmp1=0;itmp1<(int)ihist_prd.size();itmp1++){
						Punnorm[ihist_prd[itmp1]] = value;
					}
				}
			}
		}
	}
	
	return 0;
}

double histogram::CheckConvergence(ofstream &file_out,int it ,double kT){
	int ihist;
	double rmsd = 0.0;
	if(debug)file_out<<"WHAM CONVERGENCE IT = "<<it<<endl;
	ComputeEnergy(kT);
	for(ihist = 0; ihist < numhist; ihist++){ 
		if(!flaginf[ihist]){
			//rmsd += (A[ihist] - Aold[ihist]) * (A[ihist] - Aold[ihist]);
			rmsd += (P[ihist] - Pold[ihist]) * (P[ihist] - Pold[ihist]);
			Aold[ihist] = A[ihist];
			Pold[ihist] = P[ihist];
		}
		if(debug)file_out<<ihist<<" "<<Punnorm[ihist]<<" "<<P[ihist]<<" "<<numwham[ihist]
			     <<" "<<denwham[ihist]<<" "<<A[ihist]<<" "<<endl;
	}
	//rmsd = sqrt(kT*kT*rmsd/numgood);
	rmsd = sqrt(rmsd/numgood);
	return rmsd;
}

double histogram::NormalizeProbability(){
	double sum;
	int ihist,iwin;
	//Calculate normalization factor
	sum = 0.0;
	for(ihist = 0; ihist < numhist; ihist++)if(!flaginf[ihist])sum += Punnorm[ihist];
	//Normalize probabilities
	for(ihist = 0; ihist < numhist; ihist++){
		if(flaginf[ihist])P[ihist] = 0.0;
		else P[ihist] = Punnorm[ihist]/sum;
	}
	//Normalize coefficient
	for(iwin = 0;iwin < numwin; iwin++){
		F[iwin] += log(sum);
	}
	return sum;
}

int histogram::ComputeEnergy(double kT){
	int ihist;
	numgood = numbad = 0;
	for(ihist = 0; ihist < numhist; ihist++)A[ihist] = -kT*log(P[ihist]);
	for(ihist = 0; ihist < numhist; ihist++){ 
		if(!flaginf[ihist]){
		     A[ihist] = -kT*log(P[ihist]);
		     flaginf[ihist] = ( isinf(A[ihist]) || isnan(A[ihist]) );
		     if(flaginf[ihist])numbad++;
		     else numgood++;
		} else {
			numbad++;
		}
	}
	return numgood;
}

int histogram::TranslateEnergy(double kT){
	int ihist,ipos,idim,num;
	double pos,mean;
	bool translate = false,inside_ref;
	num = 0;
	mean = 0.0;
	for(idim = 0;idim < numdim; idim++){
		if(refmin[idim] != refmax[idim]) translate = true;
	}
	if (translate) {
		for(ihist = 0; ihist < numhist; ihist++){
			inside_ref = true;
			for(idim = 0;idim < numdim; idim++){
				ipos = ((int)(ihist/step[idim])) % numbin[idim];
				pos = histmin[idim] + ipos*delta[idim];
				if ( (pos < refmin[idim]) || (pos > refmax[idim]) ){
					inside_ref = false;
				}
			}
			if (inside_ref == true){
				mean += A[ihist];
				num += 1;
			}
		}
		mean = mean / num;
		for(ihist = 0; ihist < numhist; ihist++)A[ihist] -= mean;
	}
	return 0; 
}
//			END WHAM
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
//			ANALYSES
double newbarrier(vector<double> &peak, vector<double> &barrier,vector<double> &value
		,vector<bool> &raising,int ind_new,int ind_old, bool &raise){
	if(raising[ind_old]){
		if(value[ind_new] > value[ind_old]){
			raise = true;
			return  barrier[ind_old] + (value[ind_new] - value[ind_old]);
		}else{
			raise = false;
			return barrier[ind_old];
		}
	} else {
		if(value[ind_new] > value[ind_old]){
			raise = true;
			if(value[ind_new] > peak[ind_old]) return barrier[ind_old] + (value[ind_new] - value[ind_old]);
			else return (value[ind_new] - value[ind_old]) > barrier[ind_old] ? (value[ind_new] - value[ind_old]) : barrier[ind_old];
		}else{
			raise = false;
			return barrier[ind_old];
		}
	}
	cerr<<"ERROR"<<endl;
	exit(1);
}

int histogram::defineMEPmatrix4(vector<double> &A_LCL,ofstream &file_out,char cmt
		,unsigned long int ihist_str, unsigned long int ihist_end
		,bool flag_forward){
	vector<double> barrier; // MEP to that point
	vector<unsigned long int> fromwhere; // Where the MEP cames from
	vector<bool> raising; // If the energy along the MEP increasing ?
	vector<double> peak;
	bool raise;
	unsigned long int ihist,ihist_from[4];
	int i0_tmp,i1_tmp,i2_tmp,i3_tmp,i0,i1,i2,i3,ind_from;
	int i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	int i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	int i2_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	int i3_str= (int)( ihist_str / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
	int i0_end= (ihist_end % (numbin[0]+2) ) - 1;
	int i1_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	int i2_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	int i3_end= (int)( ihist_end / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
	double barrier_tmp;

	//----------------------------------------------------------------------------
	// INITIALIZATION
	for(i3=-1;i3<(numbin[3]+1);i3++){
		for(i2=-1;i2<(numbin[2]+1);i2++){
			for(i1=-1;i1<(numbin[1]+1);i1++){
				for(i0=-1;i0<(numbin[0]+1);i0++){
					peak.push_back(0.0);
					fromwhere.push_back(0);
					raising.push_back(true);
					if((i3 >= i3_str) && (i3 <= i3_end) &&
					   (i2 >= i2_str) && (i2 <= i2_end) &&
					   (i1 >= i1_str) && (i1 <= i1_end) &&
					   (i0 >= i0_str) && (i0 <= i0_end)){
						barrier.push_back(0.0);
					} else {
						barrier.push_back(INF);
					}
				}
			}
		}
	}
	peak[ihist_str] = A_LCL[ihist_str];
	barrier[ihist_str-1] = 0.0;
	// END INITIALIZATION
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// COMPUTE MEP
	cout<<"A["<<ihist_str<<"] = "<<A_LCL[ihist_str]<<endl; cout.flush();
	cout<<"A["<<ihist_end<<"] = "<<A_LCL[ihist_end]<<endl; cout.flush();
	for(i3=0;i3<numbin[3];i3++){
		for(i2=0;i2<numbin[2];i2++){
			for(i1=0;i1<numbin[1];i1++){
				for(i0=0;i0<numbin[0];i0++){
				   ihist = (i0+1) + (i1+1)*(numbin[0]+2) + (i2+1)*(numbin[0]+2)*(numbin[1]+2) 
				   	+ (i3+1)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
				   if(ihist >= ihist_str){
				      ihist_from[0] = (i0) + (i1+1)*(numbin[0]+2) +   (i2+1)*(numbin[0]+2)*(numbin[1]+2) 
				      	+   (i3+1)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
				      ihist_from[1] = (i0+1) + (i1)*(numbin[0]+2) +   (i2+1)*(numbin[0]+2)*(numbin[1]+2) 
				      	+   (i3+1)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
				      ihist_from[2] = (i0+1) +   (i1+1)*(numbin[0]+2) + (i2)*(numbin[0]+2)*(numbin[1]+2) 
				      	+   (i3+1)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
				      ihist_from[3] = (i0+1) +   (i1+1)*(numbin[0]+2) +   (i2+1)*(numbin[0]+2)*(numbin[1]+2) 
				      	+ (i3)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
				      barrier[ihist] = newbarrier(peak,barrier,A_LCL,raising,ihist,ihist_from[0],raise);
				      raising[ihist] = raise;
				      fromwhere[ihist] = ihist_from[0];
				      //cout<<ihist<<'\t'<<barrier[ihist];
				      for(ind_from = 1;ind_from < 4; ind_from++){
barrier_tmp = newbarrier(peak,barrier,A_LCL,raising,ihist,ihist_from[ind_from],raise);
i0_tmp= (ihist_from[ind_from] % (numbin[0]+2) ) - 1;
i1_tmp= (int)( (ihist_from[ind_from] % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
i2_tmp= (int)( (ihist_from[ind_from] % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
i3_tmp= (int)( ihist_from[ind_from] / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
if((i3_tmp < i3_str) || (i3_tmp > i3_end) ||
   (i2_tmp < i2_str) || (i2_tmp > i2_end) ||
   (i1_tmp < i1_str) || (i1_tmp > i1_end) ||
   (i0_tmp < i0_str) || (i0_tmp > i0_end)) barrier_tmp = INF;
	   	if(barrier_tmp < barrier[ihist]){
	   		barrier[ihist] = barrier_tmp;
	   		raising[ihist] = raise;
	   		fromwhere[ihist] = ihist_from[ind_from];
	   	}
	   	//cout<<'\t'<<barrier[ihist];
	   }
	   if(A_LCL[ihist]>peak[fromwhere[ihist]]) peak[ihist] = A_LCL[ihist];
	   else peak[ihist] = peak[fromwhere[ihist]];
	   //cout<<" FINAL"<<'\t'<<A[ihist]<<'\t'<<barrier[ihist]<<'\t'<<peak[ihist]
	   //       <<'\t'<<raising[ihist]<<'\t'<<fromwhere[ihist];
i0_tmp= (fromwhere[ihist] % (numbin[0]+2) ) - 1;
i1_tmp= (int)( (fromwhere[ihist] % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
i2_tmp= (int)( (fromwhere[ihist] % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
i3_tmp= (int)( fromwhere[ihist] / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
//if((i3_tmp>= i3_str) && (i3_tmp<= i3_end) &&
//   (i2_tmp>= i2_str) && (i2_tmp<= i2_end) &&
//   (i1_tmp>= i1_str) && (i1_tmp<= i1_end) &&
//   (i0_tmp>= i0_str) && (i0_tmp<= i0_end)){
//	cout<<" INSIDE"<<endl;
//} else {
//	cout<<" OUTSIDE"<<endl;
//}
				   }// END if(ihist >= ihist_str){
				}// END for(i0=0;i0<numbin[0];i0++){
			}// END for(i1=0;i1<numbin[1];i1++){
		}// END for(i2=0;i2<numbin[2];i2++){
	}// END for(i3=0;i3<numbin[3];i3++){
	// END COMPUTE MEP
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// OUTPUT
	file_out<<setprecision(4);
	file_out.setf(ios::fixed);
	ihist = ihist_end;
	if(flag_forward){
		file_out<<cmt<<" Forward direction"<<endl;
		if(cmt == '%') file_out<<"FRW = ["<<endl;
		cout<<"Energetic Barrier - Forward direction = "<<barrier[ihist]<<endl;
	}
	else {
		file_out<<cmt<<" Backward direction"<<endl;
		if(cmt == '%') file_out<<"BKW = ["<<endl;
		cout<<"Energetic Barrier - Backward direction = "<<barrier[ihist]<<endl;
	}
	file_out<<cmt<<" Highest energy barrier along the Minimum Energy Path = "<<barrier[ihist]<<endl;
	file_out<<cmt<<setw(9)<<" ihist"<<setw(10)
		<<"ion_0"<<setw(10)<<"ion_1"<<setw(10)<<"ion_2"<<setw(10)<<"ion_3"<<setw(10)<<"distance"
		<<setw(10)<<"A"<<setw(10)<<"barrier"<<setw(10)<<"fromwhere"<<endl;
	i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	i2_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	i3_str= (int)( ihist_str / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
	if(!flag_forward){
		i0_str = numbin[0] - i0_str - 1;
		i1_str = numbin[1] - i1_str - 1;
		i2_str = numbin[2] - i2_str - 1;
		i3_str = numbin[3] - i3_str - 1;
	}
	bool flag_endmep = false;
	bool failed = false; // This is used to avoid infinite cycle in un-connected maps
	while((!flag_endmep)&&(!failed)){
		if(ihist == ihist_str)flag_endmep = true;
		i0 = (ihist % (numbin[0]+2) ) - 1;
		i1 = (int)( (ihist % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
		i2 = (int)( (ihist % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
		i3 = (int)( ihist / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
		if(!flag_forward){
			i0 = numbin[0] - i0 - 1;
			i1 = numbin[1] - i1 - 1;
			i2 = numbin[2] - i2 - 1;
			i3 = numbin[3] - i3 - 1;
		}
		file_out<<setw(10)<<ihist<<setw(10)<<histmin[0]+i0*delta[0]<<setw(10)<<histmin[1]+i1*delta[1]
			<<setw(10)<<histmin[2]+i2*delta[2]<<setw(10)<<histmin[3]+i3*delta[3]
			<<setw(10)<<sqrt(
((histmin[0]+i0*delta[0])-(histmin[0]+i0_str*delta[0]))*((histmin[0]+i0*delta[0])-(histmin[0]+i0_str*delta[0])) + 
((histmin[1]+i1*delta[1])-(histmin[1]+i1_str*delta[1]))*((histmin[1]+i1*delta[1])-(histmin[1]+i1_str*delta[1])) +
((histmin[2]+i2*delta[2])-(histmin[2]+i2_str*delta[2]))*((histmin[2]+i2*delta[2])-(histmin[2]+i2_str*delta[2])) +
((histmin[3]+i3*delta[3])-(histmin[3]+i3_str*delta[3]))*((histmin[3]+i3*delta[3])-(histmin[3]+i3_str*delta[3])) )
			<<setw(10)<<A_LCL[ihist]<<setw(10)<<barrier[ihist]<<setw(10)<<fromwhere[ihist]<<endl;
		ihist = fromwhere[ihist];
		if(ihist <= 0){
			file_out<<cmt<<" MEP search failed !"<<endl;
			failed = true;
		}
	}
	if(cmt == '%') file_out<<"];"<<endl;
	file_out<<endl<<endl;
	//----------------------------------------------------------------------------
	  
	return 0;
}

int histogram::defineMEPmatrix3(vector<double> &A_LCL,ofstream &file_out, char cmt
		,unsigned long int ihist_str, unsigned long int ihist_end
		,bool flag_forward){
	vector<double> barrier; // MEP to that point
	vector<unsigned long int> fromwhere; // Where the MEP cames from
	vector<bool> raising; // If the energy along the MEP increasing ?
	vector<double> peak;
	bool raise;
	unsigned long int ihist,ihist_from[3];
	int i0_tmp,i1_tmp,i2_tmp,i0,i1,i2,ind_from;
	int i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	int i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	int i2_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	int i0_end= (ihist_end % (numbin[0]+2) ) - 1;
	int i1_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	int i2_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	double barrier_tmp;
	double dist_2_center_min,dist_2_center;
	int iwin_closer,iwin_closer_old,iwin;

	//----------------------------------------------------------------------------
	// INITIALIZATION
	for(i2=-1;i2<(numbin[2]+1);i2++){
		for(i1=-1;i1<(numbin[1]+1);i1++){
			for(i0=-1;i0<(numbin[0]+1);i0++){
				peak.push_back(0.0);
				fromwhere.push_back(0);
				raising.push_back(true);
				if((i2 >= i2_str) && (i2 <= i2_end) &&
				   (i1 >= i1_str) && (i1 <= i1_end) &&
				   (i0 >= i0_str) && (i0 <= i0_end)){
					barrier.push_back(0.0);
				} else {
					barrier.push_back(INF);
				}
			}
		}
	}
	peak[ihist_str] = A_LCL[ihist_str];
	barrier[ihist_str-1] = 0.0;
	// END INITIALIZATION
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// COMPUTE MEP
	cout<<"A["<<ihist_str<<"] = "<<A_LCL[ihist_str]<<endl; cout.flush();
	cout<<"A["<<ihist_end<<"] = "<<A_LCL[ihist_end]<<endl; cout.flush();
	for(i2=0;i2<numbin[2];i2++){
		for(i1=0;i1<numbin[1];i1++){
			for(i0=0;i0<numbin[0];i0++){
			   ihist = (i0+1) + (i1+1)*(numbin[0]+2) + (i2+1)*(numbin[0]+2)*(numbin[1]+2);
			   //cout<<"DEBUG1: "<<ihist<<"\t"<<ihist_str<<"\t"<<ihist_end<<"\t"<<numbin[0]<<"\t"<<numbin[1]<<"\t"<<i0<<"\t"<<i1<<"\t"<<i2<<endl;
			   if(ihist >= ihist_str){
			      ihist_from[0] = (i0)   + (i1+1)*(numbin[0]+2) + (i2+1)*(numbin[0]+2)*(numbin[1]+2);
			      ihist_from[1] = (i0+1) + (i1)  *(numbin[0]+2) + (i2+1)*(numbin[0]+2)*(numbin[1]+2);
			      ihist_from[2] = (i0+1) + (i1+1)*(numbin[0]+2) + (i2)  *(numbin[0]+2)*(numbin[1]+2);
			      barrier[ihist] = newbarrier(peak,barrier,A_LCL,raising,ihist,ihist_from[0],raise);
			      raising[ihist] = raise;
			      fromwhere[ihist] = ihist_from[0];
			      //cout<<"DEBUG2: "<<ihist<<"\t"<<ihist_from[0]<<"\t"<<ihist_from[1]<<"\t"<<ihist_from[2]<<endl;
			      //cout<<"DEBUG3: "<<ihist<<"\t"<<barrier[ihist]<<endl;
			      for(ind_from = 1;ind_from < 3; ind_from++){
					barrier_tmp = newbarrier(peak,barrier,A_LCL,raising,ihist,ihist_from[ind_from],raise);
					i0_tmp= (ihist_from[ind_from] % (numbin[0]+2) ) - 1;
					i1_tmp= (int)( (ihist_from[ind_from] % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
					i2_tmp= (int)( (ihist_from[ind_from] % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
					if((i2_tmp < i2_str) || (i2_tmp > i2_end) ||
					   (i1_tmp < i1_str) || (i1_tmp > i1_end) ||
					   (i0_tmp < i0_str) || (i0_tmp > i0_end)) barrier_tmp = INF;
					if(barrier_tmp < barrier[ihist]){
						barrier[ihist] = barrier_tmp;
						raising[ihist] = raise;
						fromwhere[ihist] = ihist_from[ind_from];
					}
			      		//cout<<"DEBUG3: "<<ihist<<"\t"<<barrier[ihist]<<"\t"<<i0_tmp<<"\t"<<i1_tmp<<"\t"<<i2_tmp<<endl;
	   			}
	   			if(A_LCL[ihist]>peak[fromwhere[ihist]]) peak[ihist] = A_LCL[ihist];
	   			else peak[ihist] = peak[fromwhere[ihist]];
				i0_tmp= (fromwhere[ihist] % (numbin[0]+2) ) - 1;
				i1_tmp= (int)( (fromwhere[ihist] % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
				i2_tmp= (int)( (fromwhere[ihist] % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
			   }// END if(ihist >= ihist_str){
			   //cout<<"DEBUG4: "<<ihist<<"\t"<<barrier[ihist]<<"\t"<<fromwhere[ihist]<<endl;
			}// END for(i0=0;i0<numbin[0];i0++){
		}// END for(i1=0;i1<numbin[1];i1++){
	}// END for(i2=0;i2<numbin[2];i2++){
	// END COMPUTE MEP
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// OUTPUT
	file_out<<setprecision(4);
	file_out.setf(ios::fixed);
	ihist = ihist_end;
	iwin_closer_old = -1;
	if(flag_forward){
		file_out<<cmt<<" Forward direction"<<endl;
		if(cmt == '%')file_out<<"FRW = ["<<endl;
		cout<<"Energetic Barrier - Forward direction = "<<barrier[ihist]<<endl;
	}
	else {
		file_out<<cmt<<" Backward direction"<<endl;
		if(cmt == '%')file_out<<"BKW = ["<<endl;
		cout<<"Energetic Barrier - Backward direction = "<<barrier[ihist]<<endl;
	}
	file_out<<cmt<<" Highest energy barrier along the Minimum Energy Path = "<<barrier[ihist]<<endl;
	file_out<<cmt<<setw(9)<<" ihist"<<setw(10)<<"ion_0"<<setw(10)<<"ion_1"<<setw(10)<<"ion_2"
		<<setw(10)<<"distance"<<setw(10)<<"A"<<setw(10)<<"barrier"<<setw(10)<<"fromwhere"<<endl;
	i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	i2_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	if(!flag_forward){
		i0_str = numbin[0] - i0_str - 1;
		i1_str = numbin[1] - i1_str - 1;
		i2_str = numbin[2] - i2_str - 1;
	}
	bool flag_endmep = false;
	bool failed = false; // This is used to avoid infinite cycle in un-connected maps
	while((!flag_endmep)&&(!failed)){
		if(ihist == ihist_str)flag_endmep = true;
		i0 = (ihist % (numbin[0]+2) ) - 1;
		i1 = (int)( (ihist % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
		i2 = (int)( (ihist % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
		if(!flag_forward){
			i0 = numbin[0] - i0 -1;
			i1 = numbin[1] - i1 -1;
			i2 = numbin[2] - i2 -1;
		}
		// Search for the closest window
		
		dist_2_center_min = INF;
		for(iwin = 0;iwin < numwin; iwin++){
			dist_2_center =
				((histmin[0] + i0*delta[0]) - center[0+iwin*numdim]) *
					((histmin[0] + i0*delta[0]) - center[0+iwin*numdim]) +
				((histmin[1] + i1*delta[1]) - center[1+iwin*numdim]) *
					((histmin[1] + i1*delta[1]) - center[1+iwin*numdim]) +
				((histmin[2] + i2*delta[2]) - center[2+iwin*numdim]) *
					((histmin[2] + i2*delta[2]) - center[2+iwin*numdim]);
			if(dist_2_center < dist_2_center_min){
				dist_2_center_min = dist_2_center;
				iwin_closer = iwin;
			}
		}
		if (iwin_closer != iwin_closer_old){
			if (flag_forward) file_out<<"#TRJ_FRW ";
			else  file_out<<"#TRJ_BKW ";
			file_out<<setw(10)<<center[0+iwin_closer*numdim]<<setw(10)<<center[1+iwin_closer*numdim]<<setw(10)<<center[2+iwin_closer*numdim]<<endl;
			iwin_closer_old = iwin_closer;
		}

		file_out<<setw(10)<<ihist<<setw(10)<<histmin[0]+i0*delta[0]<<setw(10)<<histmin[1]+i1*delta[1]
			<<setw(10)<<histmin[2]+i2*delta[2]
			<<setw(10)<<sqrt(
((histmin[0]+i0*delta[0])-(histmin[0]+i0_str*delta[0]))*((histmin[0]+i0*delta[0])-(histmin[0]+i0_str*delta[0])) + 
((histmin[1]+i1*delta[1])-(histmin[1]+i1_str*delta[1]))*((histmin[1]+i1*delta[1])-(histmin[1]+i1_str*delta[1])) +
((histmin[2]+i2*delta[2])-(histmin[2]+i2_str*delta[2]))*((histmin[2]+i2*delta[2])-(histmin[2]+i2_str*delta[2])) )
			<<setw(10)<<A_LCL[ihist]<<setw(10)<<barrier[ihist]<<setw(10)<<fromwhere[ihist]<<endl;
		ihist = fromwhere[ihist];
		if(ihist <= 0){
			file_out<<cmt<<" MEP search failed !"<<endl;
			failed = true;
		}
	}
	if(cmt == '%')file_out<<"];"<<endl;
	file_out<<endl<<endl;
	//----------------------------------------------------------------------------
	  
	return 0;
}

int histogram::defineMEPmatrix2(vector<double> &A_LCL,ofstream &file_out, char cmt
		,unsigned long int ihist_str, unsigned long int ihist_end
		,bool flag_forward){
	vector<double> barrier; // MEP to that point
	vector<unsigned long int> fromwhere; // Where the MEP cames from
	vector<bool> raising; // If the energy along the MEP increasing ?
	vector<double> peak;
	bool raise;
	unsigned long int ihist,ihist_from[3];
	int i0_tmp,i1_tmp,i0,i1,ind_from;
	int i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	int i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	int i0_end= (ihist_end % (numbin[0]+2) ) - 1;
	int i1_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	double barrier_tmp;
	double dist_2_center_min,dist_2_center;
	int iwin_closer,iwin_closer_old,iwin;

	//----------------------------------------------------------------------------
	// INITIALIZATION
	for(i1=-1;i1<(numbin[1]+1);i1++){
		for(i0=-1;i0<(numbin[0]+1);i0++){
			peak.push_back(0.0);
			fromwhere.push_back(0);
			raising.push_back(true);
			if((i1 >= i1_str) && (i1 <= i1_end) && (i0 >= i0_str) && (i0 <= i0_end)){
				barrier.push_back(0.0);
			} else {
				barrier.push_back(INF);
			}
		}
	}
	peak[ihist_str] = A_LCL[ihist_str];
	barrier[ihist_str-1] = 0.0;
	// END INITIALIZATION
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// COMPUTE MEP
	cout<<"A["<<ihist_str<<"] = "<<A_LCL[ihist_str]<<endl; cout.flush();
	cout<<"A["<<ihist_end<<"] = "<<A_LCL[ihist_end]<<endl; cout.flush();
	for(i1=0;i1<numbin[1];i1++){
		for(i0=0;i0<numbin[0];i0++){
		   ihist = (i0+1) + (i1+1)*(numbin[0]+2);
		   if(ihist >= ihist_str){
		      ihist_from[0] = (i0)   + (i1+1)*(numbin[0]+2) ;
		      ihist_from[1] = (i0+1) + (i1)  *(numbin[0]+2) ;
		      barrier[ihist] = newbarrier(peak,barrier,A_LCL,raising,ihist,ihist_from[0],raise);
		      raising[ihist] = raise;
		      fromwhere[ihist] = ihist_from[0];
		      for(ind_from = 1;ind_from < 2; ind_from++){
				barrier_tmp = newbarrier(peak,barrier,A_LCL,raising,ihist,ihist_from[ind_from],raise);
				i0_tmp= (ihist_from[ind_from] % (numbin[0]+2) ) - 1;
				i1_tmp= (int)( (ihist_from[ind_from] % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
				if((i1_tmp < i1_str) || (i1_tmp > i1_end) ||
				   (i0_tmp < i0_str) || (i0_tmp > i0_end)) barrier_tmp = INF;
				if(barrier_tmp < barrier[ihist]){
					barrier[ihist] = barrier_tmp;
					raising[ihist] = raise;
					fromwhere[ihist] = ihist_from[ind_from];
				}
			}
			if(A_LCL[ihist]>peak[fromwhere[ihist]]) peak[ihist] = A_LCL[ihist];
			else peak[ihist] = peak[fromwhere[ihist]];
			i0_tmp= (fromwhere[ihist] % (numbin[0]+2) ) - 1;
			i1_tmp= (int)( (fromwhere[ihist] % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	   		}// END if(ihist >= ihist_str){
		}// END for(i0=0;i0<numbin[0];i0++){
	}// END for(i1=0;i1<numbin[1];i1++){
	// END COMPUTE MEP
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// OUTPUT
	file_out<<setprecision(4);
	file_out.setf(ios::fixed);
	ihist = ihist_end;
	iwin_closer_old = -1;
	if(flag_forward){
		file_out<<cmt<<" Forward direction"<<endl;
		if(cmt == '%')file_out<<"FRW = ["<<endl;
		cout<<"Energetic Barrier - Forward direction = "<<barrier[ihist]<<endl;
	}
	else {
		file_out<<cmt<<" Backward direction"<<endl;
		if(cmt == '%')file_out<<"BKW = ["<<endl;
		cout<<"Energetic Barrier - Backward direction = "<<barrier[ihist]<<endl;
	}
	file_out<<cmt<<" Highest energy barrier along the Minimum Energy Path = "<<barrier[ihist]<<endl;
	file_out<<cmt<<setw(9)<<" ihist"<<setw(10)<<"ion_0"<<setw(10)<<"ion_1"<<setw(10)<<"ion_2"
		<<setw(10)<<"distance"<<setw(10)<<"A"<<setw(10)<<"barrier"<<setw(10)<<"fromwhere"<<endl;
	i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	if(!flag_forward){
		i0_str = numbin[0] - i0_str - 1;
		i1_str = numbin[1] - i1_str - 1;
	}
	bool flag_endmep = false;
	bool failed = false; // This is used to avoid infinite cycle in un-connected maps
	while((!flag_endmep)&&(!failed)){
		if(ihist == ihist_str)flag_endmep = true;
		i0 = (ihist % (numbin[0]+2) ) - 1;
		i1 = (int)( (ihist % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
		if(!flag_forward){
			i0 = numbin[0] - i0 -1;
			i1 = numbin[1] - i1 -1;
		}
		// Search for the closest window
		
		dist_2_center_min = INF;
		for(iwin = 0;iwin < numwin; iwin++){
			dist_2_center =
				((histmin[0] + i0*delta[0]) - center[0+iwin*numdim]) *
					((histmin[0] + i0*delta[0]) - center[0+iwin*numdim]) +
				((histmin[1] + i1*delta[1]) - center[1+iwin*numdim]) *
					((histmin[1] + i1*delta[1]) - center[1+iwin*numdim]);
			if(dist_2_center < dist_2_center_min){
				dist_2_center_min = dist_2_center;
				iwin_closer = iwin;
			}
		}
		if (iwin_closer != iwin_closer_old){
			if (flag_forward) file_out<<"#TRJ_FRW ";
			else  file_out<<"#TRJ_BKW ";
			file_out<<setw(10)<<center[0+iwin_closer*numdim]<<setw(10)<<center[1+iwin_closer*numdim]<<endl;
			iwin_closer_old = iwin_closer;
		}

		file_out<<setw(10)<<ihist<<setw(10)<<histmin[0]+i0*delta[0]<<setw(10)<<histmin[1]+i1*delta[1]
			<<setw(10)<<sqrt(
((histmin[0]+i0*delta[0])-(histmin[0]+i0_str*delta[0]))*((histmin[0]+i0*delta[0])-(histmin[0]+i0_str*delta[0])) + 
((histmin[1]+i1*delta[1])-(histmin[1]+i1_str*delta[1]))*((histmin[1]+i1*delta[1])-(histmin[1]+i1_str*delta[1])))
			<<setw(10)<<A_LCL[ihist]<<setw(10)<<barrier[ihist]<<setw(10)<<fromwhere[ihist]<<endl;
		ihist = fromwhere[ihist];
		if(ihist <= 0){
			file_out<<cmt<<" MEP search failed !"<<endl;
			failed = true;
		}
	}
	if(cmt == '%')file_out<<"];"<<endl;
	file_out<<endl<<endl;
	//----------------------------------------------------------------------------
	  
	return 0;
}

int histogram::FindMEP(char *fileout,vector<double> pnt_str,vector<double> pnt_end,char cmt
		, double kT, int nummin, double dist_min){
	ofstream file_out;

	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
	        return 1;
	}
	if(((int)pnt_str.size()!=numdim) || ((int)pnt_str.size()!=numdim)){
		cerr<<"ERROR: WRONG DIMENSION OF THE EXTREME POINTS"<<endl;
	        return 1;
	}
	ComputeEnergy(kT);
	if(numdim == 4){
		(*this).FindMEP4(file_out,pnt_str,pnt_end,cmt,kT,nummin,dist_min);
	} else if (numdim == 3) {
		(*this).FindMEP3(file_out,pnt_str,pnt_end,cmt,kT,nummin,dist_min);
	} else if (numdim == 2) {
		(*this).FindMEP2(file_out,pnt_str,pnt_end,cmt,kT,nummin,dist_min);
	} else {
		cerr<<"Routine FindMEP only works in 2-3-4 dimensions (because I'm lazy...)"<<endl;
		exit(1);
	}
	file_out.close();
	return 0;
}

int histogram::FindMEP4(ofstream &file_out,vector<double> pnt_str,vector<double> pnt_end,char cmt
		, double kT, int nummin, double dist_min){
	vector<double> A_fw; // Energy, with elements in the standard order
	vector<double> A_bk; // Energy, with elements in the reverse order
	int i0,i1,i2,i3,ind_min,ind_tmp,ind_min_end = 0,ind_min_str = 0;
	int i0_str,i1_str,i2_str,i3_str,i0_end,i1_end,i2_end,i3_end;
	double dist_tmp,dist_str,dist_end,MIN[nummin],i0_min[nummin],i1_min[nummin],i2_min[nummin],i3_min[nummin];
	unsigned long int ihist,ihist_min[nummin],ihist_str = 0,ihist_end = 0;
	bool flagalreadyfound;

	//---------------------------------------------------------------------
	// Initialization A_fw A_bk
	for(i3=-1;i3<(numbin[3]+1);i3++){
		for(i2=-1;i2<(numbin[2]+1);i2++){
			for(i1=-1;i1<(numbin[1]+1);i1++){
				for(i0=-1;i0<(numbin[0]+1);i0++){
					if((i3 > -1) && (i3 < numbin[3]) &&
						(i2 > -1) && (i2 < numbin[2]) &&
						(i1 > -1) && (i1 < numbin[1]) &&
						(i0 > -1) && (i0 < numbin[0])){
ihist = i0 + i1*numbin[0] + i2*numbin[0]*numbin[1] + i3*numbin[0]*numbin[1]*numbin[2];
A_fw.push_back(A[ihist]);
ihist = (numbin[0]-1-i0) + (numbin[1]-1-i1)*numbin[0] + (numbin[2]-1-i2)*numbin[0]*numbin[1] 
	+ (numbin[3]-1-i3)*numbin[0]*numbin[1]*numbin[2];
A_bk.push_back(A[ihist]);
					} else {
						A_fw.push_back(INF);
						A_bk.push_back(INF);
					}
				}
			}
		}
	}
	// END Initialization A_fw A_bk
	//---------------------------------------------------------------------
	  
	//---------------------------------------------------------------------
	// Minima search
	for(ind_min = 0; ind_min < nummin; ind_min++){
		MIN[ind_min] = INF;
		for(i3=-1;i3<(numbin[3]+1);i3++){
			for(i2=-1;i2<(numbin[2]+1);i2++){
				for(i1=-1;i1<(numbin[1]+1);i1++){
					for(i0=-1;i0<(numbin[0]+1);i0++){
						if((i3 > -1) && (i3 < numbin[3]) &&
							(i2 > -1) && (i2 < numbin[2]) &&
							(i1 > -1) && (i1 < numbin[1]) &&
							(i0 > -1) && (i0 < numbin[0])){
ihist = i0 + i1*numbin[0] + i2*numbin[0]*numbin[1] + i3*numbin[0]*numbin[1]*numbin[2];
if(A[ihist] < MIN[ind_min]){
	flagalreadyfound = false;
	for(ind_tmp = 0; ind_tmp < ind_min; ind_tmp++){
		dist_tmp = sqrt(
			 (histmin[0] + i0*delta[0] - i0_min[ind_tmp])*(histmin[0] + i0*delta[0] - i0_min[ind_tmp])
			+(histmin[1] + i1*delta[1] - i1_min[ind_tmp])*(histmin[1] + i1*delta[1] - i1_min[ind_tmp])
			+(histmin[2] + i2*delta[2] - i2_min[ind_tmp])*(histmin[2] + i2*delta[2] - i2_min[ind_tmp])
			+(histmin[3] + i3*delta[3] - i3_min[ind_tmp])*(histmin[3] + i3*delta[3] - i3_min[ind_tmp])
				);
		if(dist_tmp < dist_min){
			flagalreadyfound = true;
		}
	}
	if(!flagalreadyfound){
		MIN[ind_min] = A[ihist];
		i0_min[ind_min] = histmin[0] + i0*delta[0];
		i1_min[ind_min] = histmin[1] + i1*delta[1];
		i2_min[ind_min] = histmin[2] + i2*delta[2];
		i3_min[ind_min] = histmin[3] + i3*delta[3];
		ihist_min[ind_min] = (i0+1) + (i1+1)*(numbin[0]+2) + (i2+1)*(numbin[0]+2)*(numbin[1]+2) 
			+ (i3+1)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
	}
}
						}
					}
				}
			}
		}
	}
	file_out<<cmt<<" Energy minima further than "<<dist_min<<" A to each other"<<endl;
	file_out<<setprecision(4);
	file_out.setf(ios::fixed);
	if(cmt == '%')file_out<<"MIN = ["<<endl;
	dist_str = dist_end = INF;
	for(ind_min = 0; ind_min < nummin; ind_min++){
		file_out<<setw(10)<<i0_min[ind_min]
			<<setw(10)<<i1_min[ind_min]
			<<setw(10)<<i2_min[ind_min]
			<<setw(10)<<i3_min[ind_min]
			<<setw(10)<<MIN[ind_min]
			<<endl;
		dist_tmp = sqrt(
			 (pnt_str[0] - i0_min[ind_min])*(pnt_str[0] - i0_min[ind_min])
			+(pnt_str[1] - i1_min[ind_min])*(pnt_str[1] - i1_min[ind_min])
			+(pnt_str[2] - i2_min[ind_min])*(pnt_str[2] - i2_min[ind_min])
			+(pnt_str[3] - i3_min[ind_min])*(pnt_str[3] - i3_min[ind_min])
				);
		if(dist_tmp < dist_str){
			ihist_str = ihist_min[ind_min];
			ind_min_str = ind_min;
			dist_str = dist_tmp;
		}
		dist_tmp = sqrt(
			 (pnt_end[0] - i0_min[ind_min])*(pnt_end[0] - i0_min[ind_min])
			+(pnt_end[1] - i1_min[ind_min])*(pnt_end[1] - i1_min[ind_min])
			+(pnt_end[2] - i2_min[ind_min])*(pnt_end[2] - i2_min[ind_min])
			+(pnt_end[3] - i3_min[ind_min])*(pnt_end[3] - i3_min[ind_min])
				);
		if(dist_tmp < dist_end){
			ihist_end = ihist_min[ind_min];
			ind_min_end = ind_min;
			dist_end = dist_tmp;
		}
	}
	if(cmt == '%')file_out<<"];"<<endl;
	file_out<<cmt<<" Minima "<<ihist_str<<" ( A = "<<A_fw[ihist_str]
		<<" kcal/mol) is the closest to the requested start point"<<endl;
	file_out<<cmt<<"\t"<<i0_min[ind_min_str]<<'\t'<<i1_min[ind_min_str]
		<<'\t'<<i2_min[ind_min_str]<<'\t'<<i3_min[ind_min_str]<<endl;
	file_out<<cmt<<" Requested start point = "<<pnt_str[0]<<'\t'<<pnt_str[1]<<'\t'<<pnt_str[2]<<'\t'<<pnt_str[3]<<endl;
	file_out<<cmt<<" Minima "<<ihist_end<<" ( A = "<<A_fw[ihist_end]
		<<" kcal/mol) is the closest to the requested end point"<<endl;
	file_out<<cmt<<"\t"<<i0_min[ind_min_end]<<'\t'<<i1_min[ind_min_end]
		<<'\t'<<i2_min[ind_min_end]<<'\t'<<i3_min[ind_min_end]<<endl;
	file_out<<cmt<<" Requested end point = "<<pnt_end[0]<<'\t'<<pnt_end[1]<<'\t'<<pnt_end[2]<<'\t'<<pnt_end[3]<<endl;
	file_out<<endl<<endl;

	(*this).defineMEPmatrix4(A_fw,file_out,cmt,ihist_str,ihist_end,true); // MEP for forward transitions

	i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	i2_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	i3_str= (int)( ihist_str / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
	i0_end= (ihist_end % (numbin[0]+2) ) - 1;
	i1_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	i2_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	i3_end= (int)( ihist_end / ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2)) ) - 1;
	ihist_str = (i0_str+1) + (i1_str+1)*(numbin[0]+2) + (i2_str+1)*(numbin[0]+2)*(numbin[1]+2) 
		   	+ (i3_str+1)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
	ihist_end = (i0_end+1) + (i1_end+1)*(numbin[0]+2) + (i2_end+1)*(numbin[0]+2)*(numbin[1]+2) 
		   	+ (i3_end+1)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
	cout<<"A["<<ihist_str<<"] = "<<A_fw[ihist_str]<<" A["<<ihist_end<<"] = "<<A_fw[ihist_end]<<endl;
	ihist_str = (numbin[0]-i0_end) + (numbin[1]-i1_end)*(numbin[0]+2) 
		+ (numbin[2]-i2_end)*(numbin[0]+2)*(numbin[1]+2) 
		+ (numbin[3]-i3_end)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
	ihist_end = (numbin[0]-i0_str) + (numbin[1]-i1_str)*(numbin[0]+2) 
		+ (numbin[2]-i2_str)*(numbin[0]+2)*(numbin[1]+2) 
		+ (numbin[3]-i3_str)*(numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2);
	cout<<"A["<<ihist_end<<"] = "<<A_bk[ihist_end]<<" A["<<ihist_str<<"] = "<<A_bk[ihist_str]<<endl;
	(*this).defineMEPmatrix4(A_bk,file_out,cmt,ihist_str,ihist_end,false); // MEP for backward transitions
	return 0;
}

int histogram::FindMEP3(ofstream &file_out,vector<double> pnt_str,vector<double> pnt_end,char cmt
		, double kT, int nummin, double dist_min){
	vector<double> A_fw; // Energy, with elements in the standard order
	vector<double> A_bk; // Energy, with elements in the reverse order
	int i0,i1,i2,ind_min,ind_tmp,ind_min_end = 0,ind_min_str = 0;
	int i0_str,i1_str,i2_str,i0_end,i1_end,i2_end;
	double dist_tmp,dist_str,dist_end,MIN[nummin],i0_min[nummin],i1_min[nummin],i2_min[nummin];
	unsigned long int ihist,ihist_min[nummin],ihist_str = 0,ihist_end = 0;
	bool flagalreadyfound;
	//---------------------------------------------------------------------
	// Initialization A_fw A_bk
	for(i2=-1;i2<(numbin[2]+1);i2++){
		for(i1=-1;i1<(numbin[1]+1);i1++){
			for(i0=-1;i0<(numbin[0]+1);i0++){
				if((i2 > -1) && (i2 < numbin[2]) &&
					(i1 > -1) && (i1 < numbin[1]) &&
					(i0 > -1) && (i0 < numbin[0])){
ihist = i0 + i1*numbin[0] + i2*numbin[0]*numbin[1];
A_fw.push_back(A[ihist]);
ihist = (numbin[0]-1-i0) + (numbin[1]-1-i1)*numbin[0] + (numbin[2]-1-i2)*numbin[0]*numbin[1];
A_bk.push_back(A[ihist]);
				} else {
					A_fw.push_back(INF);
					A_bk.push_back(INF);
				}
			}
		}
	}
	// END Initialization A_fw A_bk
	//---------------------------------------------------------------------
	  
	//---------------------------------------------------------------------
	// Minima search
	for(ind_min = 0; ind_min < nummin; ind_min++){
		MIN[ind_min] = INF;
			for(i2=-1;i2<(numbin[2]+1);i2++){
				for(i1=-1;i1<(numbin[1]+1);i1++){
					for(i0=-1;i0<(numbin[0]+1);i0++){
						if((i2 > -1) && (i2 < numbin[2]) &&
							(i1 > -1) && (i1 < numbin[1]) &&
							(i0 > -1) && (i0 < numbin[0])){
ihist = i0 + i1*numbin[0] + i2*numbin[0]*numbin[1];
if(A[ihist] < MIN[ind_min]){
	flagalreadyfound = false;
	for(ind_tmp = 0; ind_tmp < ind_min; ind_tmp++){
		dist_tmp = sqrt(
			 (histmin[0] + i0*delta[0] - i0_min[ind_tmp])*(histmin[0] + i0*delta[0] - i0_min[ind_tmp])
			+(histmin[1] + i1*delta[1] - i1_min[ind_tmp])*(histmin[1] + i1*delta[1] - i1_min[ind_tmp])
			+(histmin[2] + i2*delta[2] - i2_min[ind_tmp])*(histmin[2] + i2*delta[2] - i2_min[ind_tmp])
				);
		if(dist_tmp < dist_min){
			flagalreadyfound = true;
		}
	}
	if(!flagalreadyfound){
		MIN[ind_min] = A[ihist];
		i0_min[ind_min] = histmin[0] + i0*delta[0];
		i1_min[ind_min] = histmin[1] + i1*delta[1];
		i2_min[ind_min] = histmin[2] + i2*delta[2];
		ihist_min[ind_min] = (i0+1) + (i1+1)*(numbin[0]+2) + (i2+1)*(numbin[0]+2)*(numbin[1]+2);
	}
}
						}
					}
				}
			}
	}
	file_out<<cmt<<" Energy minima further than "<<dist_min<<" A to each other"<<endl;
	dist_str = dist_end = INF;
	if(cmt == '%')file_out<<"MIN = ["<<endl;
	for(ind_min = 0; ind_min < nummin; ind_min++){
		file_out<<'\t'<<i0_min[ind_min]
			<<'\t'<<i1_min[ind_min]
			<<'\t'<<i2_min[ind_min]
			<<'\t'<<MIN[ind_min]
			<<endl;
		dist_tmp = sqrt(
			 (pnt_str[0] - i0_min[ind_min])*(pnt_str[0] - i0_min[ind_min])
			+(pnt_str[1] - i1_min[ind_min])*(pnt_str[1] - i1_min[ind_min])
			+(pnt_str[2] - i2_min[ind_min])*(pnt_str[2] - i2_min[ind_min])
				);
		if(dist_tmp < dist_str){
			ihist_str = ihist_min[ind_min];
			ind_min_str = ind_min;
			dist_str = dist_tmp;
		}
		dist_tmp = sqrt(
			 (pnt_end[0] - i0_min[ind_min])*(pnt_end[0] - i0_min[ind_min])
			+(pnt_end[1] - i1_min[ind_min])*(pnt_end[1] - i1_min[ind_min])
			+(pnt_end[2] - i2_min[ind_min])*(pnt_end[2] - i2_min[ind_min])
				);
		if(dist_tmp < dist_end){
			ihist_end = ihist_min[ind_min];
			ind_min_end = ind_min;
			dist_end = dist_tmp;
		}
	}
	if(cmt == '%')file_out<<"];"<<endl;
	file_out<<cmt<<" Minima "<<ihist_str<<" ( A = "<<A_fw[ihist_str]
		<<" kcal/mol) is the closest to the requested start point"<<endl;
	file_out<<cmt<<"\t"<<i0_min[ind_min_str]<<'\t'<<i1_min[ind_min_str]
		<<'\t'<<i2_min[ind_min_str]<<endl;
	file_out<<cmt<<" Requested start point = "<<pnt_str[0]<<'\t'<<pnt_str[1]<<'\t'<<pnt_str[2]<<endl;
	file_out<<cmt<<" Minima "<<ihist_end<<" ( A = "<<A_fw[ihist_end]
		<<" kcal/mol) is the closest to the requested end point"<<endl;
	file_out<<cmt<<"\t"<<i0_min[ind_min_end]<<'\t'<<i1_min[ind_min_end]
		<<'\t'<<i2_min[ind_min_end]<<endl;
	file_out<<cmt<<" Requested end point = "<<pnt_end[0]<<'\t'<<pnt_end[1]<<'\t'<<pnt_end[2]<<endl;
	file_out<<endl<<endl;

	(*this).defineMEPmatrix3(A_fw,file_out,cmt,ihist_str,ihist_end,true); // MEP for forward transitions

	i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	i2_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	i0_end= (ihist_end % (numbin[0]+2) ) - 1;
	i1_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	i2_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2)*(numbin[2]+2))) / ((numbin[1]+2)*(numbin[0]+2)) ) - 1;
	ihist_str = (i0_str+1) + (i1_str+1)*(numbin[0]+2) + (i2_str+1)*(numbin[0]+2)*(numbin[1]+2);
	ihist_end = (i0_end+1) + (i1_end+1)*(numbin[0]+2) + (i2_end+1)*(numbin[0]+2)*(numbin[1]+2); 
	cout<<"A["<<ihist_str<<"] = "<<A_fw[ihist_str]<<" A["<<ihist_end<<"] = "<<A_fw[ihist_end]<<endl;
	ihist_str = (numbin[0]-i0_end) + (numbin[1]-i1_end)*(numbin[0]+2) 
		+ (numbin[2]-i2_end)*(numbin[0]+2)*(numbin[1]+2); 
	ihist_end = (numbin[0]-i0_str) + (numbin[1]-i1_str)*(numbin[0]+2) 
		+ (numbin[2]-i2_str)*(numbin[0]+2)*(numbin[1]+2);
	cout<<"A["<<ihist_end<<"] = "<<A_bk[ihist_end]<<" A["<<ihist_str<<"] = "<<A_bk[ihist_str]<<endl;
	(*this).defineMEPmatrix3(A_bk,file_out,cmt,ihist_str,ihist_end,false); // MEP for backward transitions
	return 0;
}

int histogram::FindMEP2(ofstream &file_out,vector<double> pnt_str,vector<double> pnt_end,char cmt
		, double kT, int nummin, double dist_min){
	vector<double> A_fw; // Energy, with elements in the standard order
	vector<double> A_bk; // Energy, with elements in the reverse order
	int i0,i1,i2,ind_min,ind_tmp,ind_min_end = 0,ind_min_str = 0;
	int i0_str,i1_str,i2_str,i0_end,i1_end,i2_end;
	double dist_tmp,dist_str,dist_end,MIN[nummin],i0_min[nummin],i1_min[nummin],i2_min[nummin];
	unsigned long int ihist,ihist_min[nummin],ihist_str = 0,ihist_end = 0;
	bool flagalreadyfound;
	//---------------------------------------------------------------------
	// Initialization A_fw A_bk
	for(i1=-1;i1<(numbin[1]+1);i1++){
		for(i0=-1;i0<(numbin[0]+1);i0++){
			if((i1 > -1) && (i1 < numbin[1]) && (i0 > -1) && (i0 < numbin[0])){
				ihist = i0 + i1*numbin[0];
				A_fw.push_back(A[ihist]);
				ihist = (numbin[0]-1-i0) + (numbin[1]-1-i1)*numbin[0];
				A_bk.push_back(A[ihist]);
			} else {
				A_fw.push_back(INF);
				A_bk.push_back(INF);
			}
		}
	}
	// END Initialization A_fw A_bk
	//---------------------------------------------------------------------
	  
	//---------------------------------------------------------------------
	// Minima search
	for(ind_min = 0; ind_min < nummin; ind_min++){
		MIN[ind_min] = INF;
		for(i1=-1;i1<(numbin[1]+1);i1++){
			for(i0=-1;i0<(numbin[0]+1);i0++){
				if ((i1 > -1) && (i1 < numbin[1]) && (i0 > -1) && (i0 < numbin[0])){
ihist = i0 + i1*numbin[0] ;
if(A[ihist] < MIN[ind_min]){
	flagalreadyfound = false;
	for(ind_tmp = 0; ind_tmp < ind_min; ind_tmp++){
		dist_tmp = sqrt(
			 (histmin[0] + i0*delta[0] - i0_min[ind_tmp])*(histmin[0] + i0*delta[0] - i0_min[ind_tmp])
			+(histmin[1] + i1*delta[1] - i1_min[ind_tmp])*(histmin[1] + i1*delta[1] - i1_min[ind_tmp])
				);
		if(dist_tmp < dist_min){
			flagalreadyfound = true;
		}
	}
	if(!flagalreadyfound){
		MIN[ind_min] = A[ihist];
		i0_min[ind_min] = histmin[0] + i0*delta[0];
		i1_min[ind_min] = histmin[1] + i1*delta[1];
		ihist_min[ind_min] = (i0+1) + (i1+1)*(numbin[0]+2);
	}
}
				}
			}
		}
	}
	file_out<<cmt<<" Energy minima further than "<<dist_min<<" A to each other"<<endl;
	dist_str = dist_end = INF;
	if(cmt == '%')file_out<<"MIN = ["<<endl;
	for(ind_min = 0; ind_min < nummin; ind_min++){
		file_out<<'\t'<<i0_min[ind_min]
			<<'\t'<<i1_min[ind_min]
			<<'\t'<<MIN[ind_min]
			<<endl;
		dist_tmp = sqrt(
			 (pnt_str[0] - i0_min[ind_min])*(pnt_str[0] - i0_min[ind_min])
			+(pnt_str[1] - i1_min[ind_min])*(pnt_str[1] - i1_min[ind_min])
				);
		if(dist_tmp < dist_str){
			ihist_str = ihist_min[ind_min];
			ind_min_str = ind_min;
			dist_str = dist_tmp;
		}
		dist_tmp = sqrt(
			 (pnt_end[0] - i0_min[ind_min])*(pnt_end[0] - i0_min[ind_min])
			+(pnt_end[1] - i1_min[ind_min])*(pnt_end[1] - i1_min[ind_min])
				);
		if(dist_tmp < dist_end){
			ihist_end = ihist_min[ind_min];
			ind_min_end = ind_min;
			dist_end = dist_tmp;
		}
	}
	if(cmt == '%')file_out<<"];"<<endl;
	file_out<<cmt<<" Minima "<<ihist_str<<" ( A = "<<A_fw[ihist_str]
		<<" kcal/mol) is the closest to the requested start point"<<endl;
	file_out<<cmt<<"\t"<<i0_min[ind_min_str]<<'\t'<<i1_min[ind_min_str]<<endl;
	file_out<<cmt<<" Requested start point = "<<pnt_str[0]<<'\t'<<pnt_str[1]<<endl;
	file_out<<cmt<<" Minima "<<ihist_end<<" ( A = "<<A_fw[ihist_end]
		<<" kcal/mol) is the closest to the requested end point"<<endl;
	file_out<<cmt<<"\t"<<i0_min[ind_min_end]<<'\t'<<i1_min[ind_min_end]<<endl;
	file_out<<cmt<<" Requested end point = "<<pnt_end[0]<<'\t'<<pnt_end[1]<<endl;
	file_out<<endl<<endl;

	(*this).defineMEPmatrix2(A_fw,file_out,cmt,ihist_str,ihist_end,true); // MEP for forward transitions

	i0_str= (ihist_str % (numbin[0]+2) ) - 1;
	i1_str= (int)( (ihist_str % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	i0_end= (ihist_end % (numbin[0]+2) ) - 1;
	i1_end= (int)( (ihist_end % ((numbin[0]+2)*(numbin[1]+2))) / (numbin[0]+2) ) - 1;
	ihist_str = (i0_str+1) + (i1_str+1)*(numbin[0]+2);
	ihist_end = (i0_end+1) + (i1_end+1)*(numbin[0]+2); 
	cout<<"A["<<ihist_str<<"] = "<<A_fw[ihist_str]<<" A["<<ihist_end<<"] = "<<A_fw[ihist_end]<<endl;
	ihist_str = (numbin[0]-i0_end) + (numbin[1]-i1_end)*(numbin[0]+2); 
	ihist_end = (numbin[0]-i0_str) + (numbin[1]-i1_str)*(numbin[0]+2);
	cout<<"A["<<ihist_end<<"] = "<<A_bk[ihist_end]<<" A["<<ihist_str<<"] = "<<A_bk[ihist_str]<<endl;
	(*this).defineMEPmatrix2(A_bk,file_out,cmt,ihist_str,ihist_end,false); // MEP for backward transitions
	return 0;
}
//			END ANALYSES
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
//			OUTPUT

////////////////////////////////////////////////////////////
//			OUTPUT
int histogram::WriteHeading(ofstream &file_out){
	int idim, ibin;
	file_out<<"# numdim = "<<numdim<<endl;
	for(idim = 0;idim < numdim; idim++){
		file_out<<"# dim = "<<idim
			<<" min = "<<histmin[idim]
			<<" max = "<<histmax[idim]
			<<" numbin = "<<numbin[idim]<<endl;
		if(debug){
			file_out<<"# delta = "<<delta[idim]<<'\t';
			for(ibin = 0; ibin <numbin[idim]; ibin++){
				file_out<<histmin[idim] - delta[idim]/2.0 + ibin*delta[idim]<<'\t';
			}
			file_out<<endl;
		}

	}
	file_out<<"# numwin = "<<numwin<<endl;
	return 0;
}

int histogram::WriteRestart(ofstream &file_out){
	file_out.write((char *)&numdim,sizeof(int));
	file_out.write((char *)&numwin,sizeof(int));
	file_out.write((char *)histmin,numdim*sizeof(double));
	file_out.write((char *)histmax,numdim*sizeof(double));
	file_out.write((char *)numbin,numdim*sizeof(int));
	file_out.write((char *)periodic,numdim*sizeof(int));
	file_out.write((char *)whr_prd,numdim*sizeof(double));
	file_out.write((char *)center,numdim*numwin*sizeof(double));
	file_out.write((char *)harmrest,numdim*numwin*sizeof(double));
	file_out.write((char *)numbinwin,numwin*sizeof(int));
	file_out.write((char *)hist,numhist*numwin*sizeof(int));
	file_out.write((char *)P,numhist*sizeof(double));
	file_out.write((char *)F,numwin*sizeof(double));
	return 0;
}

int histogram::MemorizeEne(double kT){
	int ihist;
	vector <double> vtmp;
	//NormalizeProbability();
	ComputeEnergy(kT);
	//Translate energy
	TranslateEnergy(kT);
	for(ihist = 0; ihist < numhist; ihist++){
		vtmp.push_back(A[ihist]);
	}
	backup_energy_reiter.push_back(vtmp);
	return 0;
}

// Output of the energy resulting from a Bootstrap/Split analysis
int histogram::WriteEne(ofstream &file_out,int num_reiter){
	int idim,ihist,ipos,ind_reiter;
	vector <double> A_ave,A_std;
	//Calculate averages
	for(ihist = 0; ihist < numhist; ihist++){
		A_ave.push_back(0.0);
		A_std.push_back(0.0);
		for(ind_reiter = 0; ind_reiter < num_reiter; ind_reiter++){
			A_ave[ihist] += backup_energy_reiter[ind_reiter][ihist];
		}
		A_ave[ihist] = A_ave[ihist] / num_reiter;
		for(ind_reiter = 0; ind_reiter < num_reiter; ind_reiter++){
			A_std[ihist] += ( (A_ave[ihist] - backup_energy_reiter[ind_reiter][ihist])
					  * (A_ave[ihist] - backup_energy_reiter[ind_reiter][ihist]) );
		}
		A_std[ihist] = sqrt(A_std[ihist]/(num_reiter-1));
	}
	//Output Data
	file_out<<endl<<endl<<"# Number of iterations = "<<num_reiter<<endl;
	file_out<<"# First "<<numdim<<" columns are the coordinates"<<endl;
	file_out<<"# Then: Average_Energy Standard_Deviation"<<endl;
	for(ihist = 0; ihist < numhist; ihist++){
		for(idim = 0;idim < numdim; idim++){
			ipos = ((int)(ihist/step[idim])) % numbin[idim];
			file_out<<histmin[idim] + ipos*delta[idim]<<'\t';
		}
		file_out<<A_ave[ihist]<<'\t'<<A_std[ihist]<<endl;
		if((ihist!=0) && (((ihist+1)%numbin[0]) == 0))file_out<<endl;
	}
	return 0;
}

// Output of the energy of each iteration (or of a single run)
int histogram::WriteEne(char *fileene, int ind_reiter, double kT){
	int MAXCHAR = 10000;
	char stmp[MAXCHAR+1];
	int idim,ihist,ipos,iwin,iwin_closer;
	double stdP_FB,stdP_FS,dist_2_center,dist_2_center_min,sum;
	ifstream file_in;
	ofstream file_out;

	if(ind_reiter > 0){// The file already exists so copy it
		file_in.open(fileene);
           	file_out.open(".copy_village");
                if(!file_out) {
                	cerr <<"ERROR in opening temporary file for writing"<<endl;
                        exit(1);
                }
		while(!file_in.eof()){
			file_in.getline(stmp,MAXCHAR);
			if( (stmp[0]!='#') && (strlen(stmp)>1) ){ // copy all the data lines, not the comment lines
				file_out<<stmp<<endl;
			}
		}
		file_in.close();
		file_in.clear();
		if(remove(fileene)){
                	cerr <<"ERROR deleting file "<<fileene<<endl;
                        exit(1);
		}
		file_in.open(".copy_village"); // open to read it
		if(!file_in){
			cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE TEMPORATY FILE .copy_village "<<endl;
			exit(1);
		}
		file_out.close();
		file_out.clear();
	}
	// Write the heading
	file_out.open(fileene);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileene<<" FOR OUTPUT"<<endl;
		exit(1);
	}
	WriteHeading(file_out);
	file_out<<"# First "<<numdim<<" columns are the coordinates"<<endl;
	file_out<<"# kT = "<<kT<<" kcal/mol"<<endl;
	file_out<<"# Then: Data in iteration i with the format:"<<endl;
	file_out<<"#       A[ihist]"<<endl;
	file_out<<"#       Number of samples"<<endl;
	file_out<<"#       -kT*log(P[ihist]+stdP_FB)"<<endl;
	file_out<<"#       -kT*log(P[ihist]-stdP_FB)"<<endl;
	file_out<<"#       P[ihist]"<<endl;
	file_out<<"#       stdP_FB"<<endl;
	file_out<<"#       stdP_FS"<<endl;
	// Prepare the energy
	ComputeEnergy(kT);
	//Output Data
	for(ihist = 0; ihist < numhist; ihist++){
		if(ind_reiter > 0){ // if exist, copy the previous part of the line
			file_in.getline(stmp,MAXCHAR);
			file_out<<stmp<<'\t';
		} else { // otherwise, write the coordinates (ONE FOR ALL)
			for(idim = 0;idim < numdim; idim++){
				ipos = ((int)(ihist/step[idim])) % numbin[idim];
				file_out<<histmin[idim] + ipos*delta[idim]<<'\t';
			}
		}
		// Search for the closer window
		dist_2_center_min = INF;
		iwin_closer = 0;
		sum = 0;
		for(iwin = 0;iwin < numwin; iwin++){
			dist_2_center = 0;
			for(idim = 0;idim < numdim; idim++){
				ipos = ((int)(ihist/step[idim])) % numbin[idim];
				dist_2_center +=
					((histmin[idim] + ipos*delta[idim]) - center[idim+iwin*numdim]) *
					((histmin[idim] + ipos*delta[idim]) - center[idim+iwin*numdim]);
			}
			if(dist_2_center < dist_2_center_min){
				dist_2_center_min = dist_2_center;
				iwin_closer = iwin;
			}
			sum += hist[ihist+iwin*numhist];
		}
		stdP_FB = CalculateDeviationProb(ihist,iwin_closer); // Formula Furini-Barbini
		stdP_FS = CalculateDeviationProb(ihist); // Formula Ferrenberg-Swendsen
		file_out<<A[ihist]<<'\t'<<sum<<'\t'<<-kT*log(P[ihist]+stdP_FB)<<'\t'<<-kT*log(P[ihist]-stdP_FB)
			<<'\t'<<P[ihist]<<'\t'<<stdP_FB<<'\t'<<stdP_FS<<endl;
		if((ihist!=0) && (((ihist+1)%numbin[0]) == 0))file_out<<endl;
	}
	file_in.close();
	file_out.close();
	remove(".copy_village");
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Output of the histogram in DX format
int histogram::WriteDxHist(char *fileout){
	if (numdim != 3){
		cerr<<"WARNING: DX output only for 3D files "<<fileout<<endl;
		return 1;
	}
	ofstream file_out;
	int i,j,k,ihist,sum,iwin;
	ComputeEnergy();
	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
	        return 1;
	}
	file_out<<setiosflags(ios::fixed);
	file_out<<setprecision(3);
	file_out<<"object 1 class gridpositions counts "<<numbin[2]<<" "<<numbin[1]<<" "<<numbin[0]<<endl;
	file_out<<"origin "<<histmin[0]<<" "<<histmin[1]<<" "<<histmin[2]<<endl;
	file_out<<"delta 0.000000e+00 0.000000e+00 "<<delta[2]<<endl;
	file_out<<"delta 0.000000e+00 "<<delta[1]<<" 0.000000e+00"<<endl;
	file_out<<"delta "<<delta[0]<<" 0.000000e+00 0.000000e+00"<<endl;
	file_out<<"object 2 class gridconnections counts "<<numbin[2]<<" "<<numbin[1]<<" "<<numbin[0]<<endl;
	file_out<<"object 3 class array type double rank 0 items "<<numbin[0]*numbin[1]*numbin[2]<<" data follows"<<endl;
	ihist = 0;
	for(k = 0; k < numbin[2]; k++){
		for(j = 0; j < numbin[1]; j++){
			for(i = 0; i < numbin[0]; i++){
				sum = 0;
				for(iwin = 0;iwin < numwin; iwin++){
					sum += hist[ihist+iwin*numhist];
				}
				file_out<<setw(8)<<sum<<"\t";
				ihist++;
				if((ihist%3) == 0)file_out<<endl;
			}
		}
	}
	file_out<<endl;
	file_out<<"attribute \"dep\" string \"positions\""<<endl;
	file_out<<"object \"regular positions regular connections\" class field"<<endl;
	file_out<<"component \"positions\" value 1"<<endl;
	file_out<<"component \"connections\" value 2"<<endl;
	file_out<<"component \"data\" value 3"<<endl;
	file_out<<endl;
	file_out.close();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Output of the histogram in DX format
int histogram::WriteDxEne(char *fileout,double kT){
	if (numdim != 3){
		cerr<<"WARNING: DX output only for 3D files "<<fileout<<endl;
		return 1;
	}
	ofstream file_out;
	int i,j,k,ihist;
	ComputeEnergy(kT);
	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
	        return 1;
	}
	file_out<<setiosflags(ios::fixed);
	file_out<<setprecision(3);
	file_out<<"object 1 class gridpositions counts "<<numbin[2]<<" "<<numbin[1]<<" "<<numbin[0]<<endl;
	file_out<<"origin "<<histmin[0]<<" "<<histmin[1]<<" "<<histmin[2]<<endl;
	file_out<<"delta 0.000000e+00 0.000000e+00 "<<delta[2]<<endl;
	file_out<<"delta 0.000000e+00 "<<delta[1]<<" 0.000000e+00"<<endl;
	file_out<<"delta "<<delta[0]<<" 0.000000e+00 0.000000e+00"<<endl;
	file_out<<"object 2 class gridconnections counts "<<numbin[2]<<" "<<numbin[1]<<" "<<numbin[0]<<endl;
	file_out<<"object 3 class array type double rank 0 items "<<numbin[0]*numbin[1]*numbin[2]<<" data follows"<<endl;
	ihist = 0;
	for(k = 0; k < numbin[2]; k++){
		for(j = 0; j < numbin[1]; j++){
			for(i = 0; i < numbin[0]; i++){
		     		if ( isinf(A[ihist]) || isnan(A[ihist]) ) {
					file_out<<setw(8)<<100.0<<"\t";
				} else {
					file_out<<setw(8)<<A[ihist]<<"\t";
				}
				ihist++;
				if((ihist%3) == 0)file_out<<endl;
			}
		}
	}
	file_out<<endl;
	file_out<<"attribute \"dep\" string \"positions\""<<endl;
	file_out<<"object \"regular positions regular connections\" class field"<<endl;
	file_out<<"component \"positions\" value 1"<<endl;
	file_out<<"component \"connections\" value 2"<<endl;
	file_out<<"component \"data\" value 3"<<endl;
	file_out<<endl;
	file_out.close();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Output of the energy in Matlab format
// Xi  N-dimension array, with N = number of point in dimension i
// PMF(N1,N2,N3,..) N1xN2xN3x... Multimendional matrix of the energy
int histogram::WriteEneMat(char *fileout, double kT){
	int idim,ihist,ipos,ibin, iwin;
	ofstream file_out;
	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
	        return 1;
	}
	ComputeEnergy(kT);
	for(ihist = 0; ihist < numhist; ihist++){ 
		for(iwin = 0; iwin<numwin; iwin++){
			numwham[ihist] += hist[ihist+iwin*numhist];
		}
	}
	if(numdim < 3){ // Output data for less that 3 dimensions
		for(idim = 0;idim < numdim; idim++){
			file_out<<"X"<<idim<<" = ["<<endl;
			for(ibin = 0;ibin < numbin[idim]; ibin++){
				file_out<<'\t'<<histmin[idim] + ibin*delta[idim];
			}
			file_out<<"];"<<endl;
		}
		file_out<<"POS = ["<<endl;
		for(ihist = 0; ihist < numhist; ihist++){
			for(idim = 0;idim < numdim; idim++){
				ipos = ((int)(ihist/step[idim])) % numbin[idim];
				file_out<<histmin[idim] + ipos*delta[idim]<<'\t';
			}
			file_out<<endl;
		}
		file_out<<"];"<<endl<<"PMF = ["<<endl;
		for(ihist = 0; ihist < numhist; ihist++){
			file_out<<A[ihist]<<'\t';
			if((ihist!=0) && (((ihist+1)%numbin[0]) == 0))file_out<<endl;
		}
		file_out<<"];"<<endl;
		return 0;
	} else {
		for(idim = 0;idim < numdim; idim++){
			file_out<<"X"<<idim<<" = ["<<endl;
			for(ibin = 0;ibin < numbin[idim]; ibin++){
				file_out<<'\t'<<histmin[idim] + ibin*delta[idim];
			}
			file_out<<"];"<<endl;
		}
		// Writing PMF data
		for(ihist = 0; ihist < numhist; ihist++){
			if(((ihist)%(numbin[1]*numbin[0])) == 0){// begin of a block of the first two dimensions
				file_out<<"PMF(:,:";
				for(idim = 2;idim < numdim; idim++){
					ipos = ((int)(ihist/step[idim])) % numbin[idim];
					file_out<<","<<ipos+1;
				}
				file_out<<") = ["<<endl;
			}
			file_out<<'\t'<<A[ihist];
			if(((ihist+1)%(numbin[0])) == 0)file_out<<endl;
			if(((ihist+1)%(numbin[1]*numbin[0])) == 0){ // end of a block of the first two dimensions
				file_out<<"];"<<endl;
			}
		}
		// Writing number of samples data
		for(ihist = 0; ihist < numhist; ihist++){
			if(((ihist)%(numbin[1]*numbin[0])) == 0){// begin of a block of the first two dimensions
				file_out<<"SMP(:,:";
				for(idim = 2;idim < numdim; idim++){
					ipos = ((int)(ihist/step[idim])) % numbin[idim];
					file_out<<","<<ipos+1;
				}
				file_out<<") = ["<<endl;
			}
			file_out<<'\t'<<numwham[ihist];
			if(((ihist+1)%(numbin[0])) == 0)file_out<<endl;
			if(((ihist+1)%(numbin[1]*numbin[0])) == 0){ // end of a block of the first two dimensions
				file_out<<"];"<<endl;
			}
		}
	}

	file_out.close();
	return 0;
}

// Standard deviation on unbiased probability calculated with 
// the formula by Furini-Barbini
double histogram::CalculateDeviationProb(int ihist,int iwin){
	int jwin;
	double U,den,num;
	den = 0.0;
	for(jwin = 0; jwin<numwin; jwin++){
		U = bias(ihist,jwin);
		den += ( (numbinwin[jwin]/g[jwin]) * exp(F[jwin] - U) );
	}
	U = bias(ihist,iwin);
	num = exp(U - F[iwin]) * ( 1.0*hist[ihist+iwin*numhist] / (1.0*numbinwin[iwin]) );
	return sqrt ( num/den );
}

// Standard deviation on unbiased probability calculated with 
// the formula by Ferrenberg-Swendesen
double histogram::CalculateDeviationProb(int ihist){
	int jwin;
	double den;
	den = 0.0;
	for(jwin = 0; jwin<numwin; jwin++){
		den += (hist[ihist+jwin*numhist]/g[jwin]);
	}
	return  P[ihist] / sqrt(den);
}

int histogram::WriteHist(ofstream &file_out,bool histprj,double kT){
	int idim,jdim,ihist,ipos,iwin,ibin,jbin,sum;
	double pos,P_total;

	WriteHeading(file_out);
	if(histprj){
		for(idim = 0;idim < numdim; idim++){
			for(ibin = 0;ibin < numbin[idim]; ibin++){
				pos = histmin[idim]+ibin*delta[idim];
				file_out<<pos<<'\t';
				sum = 0;
				for(ihist = 0; ihist< numhist; ihist++){
					jbin = ((int)(ihist/step[idim])) % numbin[idim];
					if(jbin == ibin){
						sum += hist[ihist];
					}
				}
				file_out<<sum<<endl;
			}
			file_out<<endl<<endl;
		}
	} else {
		file_out<<"# numdim = "<<numdim<<endl;
		for(jdim = 0;jdim < numdim; jdim++){
			file_out<<"# dim = "<<jdim
				<<" min = "<<histmin[jdim]
				<<" max = "<<histmax[jdim]
				<<" numbin = "<<numbin[jdim]<<endl;

		}
		file_out<<"# numwin = "<<numwin<<endl;
		P_total = 0.0;
		for(ihist = 0; ihist < numhist; ihist++){
			P_total += P[ihist];
		}
		file_out<<"# P_total = "<<P_total<<endl;
		for(iwin = 0;iwin < numwin; iwin++){
			P_total = 0.0;
			for(ihist = 0; ihist < numhist; ihist++){
				P_total += (1.0*hist[ihist+iwin*numhist])/(1.0*numbinwin[iwin]);
			}
			file_out<<"# P_total["<<iwin<<"] = "<<P_total<<" F["<<iwin<<"] = "<<F[iwin]<<endl;
		}
		for(ihist = 0; ihist < numhist; ihist++){
			if(((ihist%numbin[0])==0)&&(ihist!=0))file_out<<endl;
			sum = 0;
			for(iwin = 0;iwin < numwin; iwin++){
				sum += hist[ihist+iwin*numhist];
			}
			// Output coordinates
			for(jdim = 0;jdim < numdim; jdim++){
				ipos = ((int)(ihist/step[jdim])) % numbin[jdim];
				file_out<<histmin[jdim]+ipos*delta[jdim]<<'\t';
			}
			// Output cumulative values
			file_out<<sum<<'\t'<<P[ihist]<<'\t';
			// Output window's values
			for(iwin = 0;iwin < numwin; iwin++){
				file_out<<hist[ihist+iwin*numhist]<<'\t';
				file_out<<CalculateDeviationProb(ihist,iwin)<<'\t';
			}
			file_out<<endl;
		}
	}
	return 0;
}

int histogram::WriteBias(ofstream &file_out, double kT){
	int iwin,jdim,ihist,ipos;
	double U;
	for(iwin = 0; iwin < numwin; iwin++){
		file_out<<"# Centers = ";
		for(jdim = 0;jdim < numdim; jdim++){
			file_out<<center[jdim+iwin*numdim]<<'\t';
		}
		file_out<<endl<<"# Restraints = ";
		for(jdim = 0;jdim < numdim; jdim++){
			file_out<<harmrest[jdim+iwin*numdim]<<'\t';
		}
		file_out<<endl;
		for(ihist = 0; ihist < numhist; ihist++){
			if(((ihist%numbin[0])==0)&&(ihist!=0)){
				file_out<<endl;
			}
			for(jdim = 0;jdim < numdim; jdim++){
				ipos = ((int)(ihist/step[jdim])) % numbin[jdim];
				file_out<<histmin[jdim]+ipos*delta[jdim]<<'\t';
			}
			U = bias(ihist,iwin);
			file_out<<U<<endl;
		}
		file_out<<endl<<endl;
	}
	return 0;
}

int histogram::WriteProjection(char *fileout, int prjdim, bool collect_samples,  double kT){
	ofstream file_out;
	file_out.open(fileout);
	if(!file_out){
		cerr<<"ERROR: IT IS NOT POSSIBLE TO OPEN THE FILE: "<<fileout<<endl;
	        return 1;
	}
	if     ((numdim == 4)&&(prjdim == 2)) prj4D2D(file_out,collect_samples,kT);
	else if((numdim == 4)&&(prjdim == 1)) prj4D1D(file_out,kT);
	else if((numdim == 3)&&(prjdim == 2)) prj3D2D(file_out,kT);
	else if((numdim == 3)&&(prjdim == 1)) prj3D1D(file_out,kT);
	else if((numdim == 2)&&(prjdim == 1)) prj2D1D(file_out,kT);
	else {
		cerr<<"ERROR: impossible projection"<<endl;
		exit(1);
	}
	file_out.close();
	return 0;
}
//			END OUTPUT
////////////////////////////////////////////////////////////
  
////////////////////////////////////////////////////////////
//			ANALYSES START
int histogram::prj4D2D(ofstream &file_out, bool collect_samples, double kT){
	int dprj,p1,p2,d1,d2,i,k,j,l,ind,indprj,iwin,*SMP=NULL;
	double *A, *Mprj, Mtot, x, y;

	p1 = p2 = d1 = d2 = ind = -1;
	for(dprj=0;dprj<3;dprj++){
		switch(dprj) {
			case 0:
				p1 = 0; p2 = 1; d1 = 2; d2 = 3;
				break;
			case 1:
				p1 = 1; p2 = 2; d1 = 0; d2 = 3;
				break;
			case 2:
				p1 = 2; p2 = 3; d1 = 0; d2 = 1;
				break;
		}
		//Allocate memory
		Mprj = new double [ numbin[p1]*numbin[p2] ];
		A = new double [ numbin[p1]*numbin[p2] ]; // Free energy
		if(collect_samples)SMP = new int [ numbin[p1]*numbin[p2] ]; // Number of samples

		//Projecting data
		for(k = 0; k < numbin[p2]; k++){
			for(j = 0; j < numbin[p1]; j++){
				indprj = k*numbin[p1] + j;
				Mprj[indprj] = A[indprj] = 0.0; //initialize
				if(collect_samples)SMP[indprj] = 0; //initialize
				for(l = 0; l < numbin[d2]; l++){
					for(i = 0; i < numbin[d1]; i++){
						switch(dprj){
							case 0:
								ind = l*numbin[2]*numbin[1]*numbin[0] +
									i*numbin[1]*numbin[0] +
									k*numbin[0] + j;
								break;
							case 1:
								ind = l*numbin[2]*numbin[1]*numbin[0] +
									k*numbin[1]*numbin[0] +
									j*numbin[0] + i;
								break;
							case 2:
								ind = k*numbin[2]*numbin[1]*numbin[0] +
									j*numbin[1]*numbin[0] +
									l*numbin[0] + i;
								break;
						}
						if(collect_samples){
							for(iwin = 0; iwin < numwin; iwin++){
								SMP[indprj] += hist[ind+iwin*numhist];
							}
						}
						if(!isinf(P[ind])){
							Mprj[indprj] += P[ind];
						}
					}
				}
			}
		}

		//Normalize
		Mtot = 0.0;
		indprj = 0;
		for(k = 0; k < numbin[p2]; k++){
			for(j = 0; j < numbin[p1]; j++){
				Mtot += Mprj[indprj];
				indprj++;
			}
		}
		cout<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		file_out<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		indprj = 0;
		for(k = 0; k < numbin[p2]; k++){
			for(j = 0; j < numbin[p1]; j++){
				Mprj[indprj] = Mprj[indprj]/Mtot;
				A[indprj] = -kT*log(Mprj[indprj]);
				indprj++;
			}
		}

		//Output
		indprj = 0;
		for(k = 0; k < numbin[p2]; k++){
			y = histmin[p2] + delta[p2]*k;
			for(j = 0; j < numbin[p1]; j++){
				x = histmin[p1] + delta[p1]*j;
				file_out<<x<<'\t'<<y<<'\t'<<Mprj[indprj]<<'\t'<<A[indprj];
				if(collect_samples)file_out<<'\t'<<SMP[indprj]<<endl;
				else file_out<<endl;
				indprj++;
			}
			file_out<<endl;
		}
		file_out<<endl<<endl<<endl;

		//deallocate
		delete[] A;
		if(collect_samples)delete[] SMP;
		delete[] Mprj;
	}
	return 0;
}

int histogram::prj4D1D(ofstream &file_out, double kT){
	int dprj,iprj,d1,d2,d3,i,k,j,ind,ntmp,itmp;
	double *A, *Mprj, Mtot, x, min, max, deltatmp, left, right;
	double *Pprjs,*Ptmp;

	ntmp = 50;
	Pprjs = new double [numbin[0]+numbin[1]+numbin[2]+numbin[3]];
	Ptmp = new double [ntmp];

	min = INF;
	max = -INF;
	d1 = d2 = d3 = ind = -1;
	for(dprj=0;dprj<4;dprj++){
		//allocate
		Mprj = new double [ numbin[dprj] ];
		A = new double [ numbin[dprj] ];
		//Projecting data
		for(iprj=0;iprj<numbin[dprj];iprj++){
			Mprj[iprj] = A[iprj] = 0.0; //initialize
			switch (dprj) {
				case 0:
					d1 = 1; d2 = 2; d3 = 3;
					break;
				case 1:
					d1 = 0; d2 = 2; d3 = 3;
					break;
				case 2:
					d1 = 0; d2 = 1; d3 = 3;
					break;
				case 3:
					d1 = 0; d2 = 1; d3 = 2;
					break;
			}
			for(k = 0; k < numbin[d3]; k++){
				for(j = 0; j < numbin[d2]; j++){
					for(i = 0; i < numbin[d1]; i++){
						switch (dprj) {
							case 0:
								ind = k*numbin[2]*numbin[1]*numbin[0] +
									j*numbin[1]*numbin[0] +
									i*numbin[0] + iprj;
								break;
							case 1:
								ind = k*numbin[2]*numbin[1]*numbin[0] +
									j*numbin[1]*numbin[0] +
									iprj*numbin[0] + i;
								break;
							case 2:
								ind = k*numbin[2]*numbin[1]*numbin[0] +
									iprj*numbin[1]*numbin[0] +
									j*numbin[0] + i;
								break;
							case 3:
								ind = iprj*numbin[2]*numbin[1]*numbin[0] +
									k*numbin[1]*numbin[0] +
									j*numbin[0] + i;
								break;
						}
						if(!isinf(P[ind])){
							Mprj[iprj] += P[ind];
						}
					}
				}
			}
		}
		//Normalize
		Mtot = 0.0;
		for(iprj=0;iprj<numbin[dprj];iprj++)Mtot += Mprj[iprj];
		cout<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		file_out<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		for(iprj=0;iprj<numbin[dprj];iprj++){
			Mprj[iprj] = Mprj[iprj]/Mtot;
			A[iprj] = -kT*log(Mprj[iprj]);
		}
		//Output
		file_out<<"# Projection along dimension "<<dprj<<endl;
		for(iprj=0;iprj<numbin[dprj];iprj++){
			x = histmin[dprj] + delta[dprj]*iprj;
			if(x < min) min = x;
			if(x > max) max = x;
			file_out<<x<<'\t'<<Mprj[iprj]<<'\t'<<A[iprj]<<endl;
		}
		for(iprj = 0; iprj < numbin[dprj]; iprj++){
			switch(dprj){
				case 0:
					ind = iprj;
					break;
				case 1:
					ind = numbin[0]+iprj;
					break;
				case 2:
					ind = numbin[0]+numbin[1]+iprj;
					break;
				case 3:
					ind = numbin[0]+numbin[1]+numbin[2]+iprj;
					break;
			}
			Pprjs[ind] = Mprj[iprj];
		}
		file_out<<endl<<endl;
		//deallocate
		delete[] A;
		delete[] Mprj;
	}

	deltatmp = (max - min)/(ntmp-1);
	for(itmp = 0;itmp < ntmp; itmp++){
		Ptmp[itmp] = 0.0;
		left = min - (deltatmp/2.0) + itmp * deltatmp;
		right = min + (deltatmp/2.0) + itmp * deltatmp;
		for(dprj = 0; dprj < 4; dprj++){
			for(iprj = 0; iprj < numbin[dprj]; iprj++){
				x = histmin[dprj] + delta[dprj]*iprj;
				if((x >= left) && (x<=right) ){ 
					switch(dprj){
						case 0:
							ind = iprj;
							break;
						case 1:
							ind = numbin[0]+iprj;
							break;
						case 2:
							ind = numbin[0]+numbin[1]+iprj;
							break;
						case 3:
							ind = numbin[0]+numbin[1]+numbin[2]+iprj;
							break;
					}
					Ptmp[itmp] += Pprjs[ind];
				}
			}
		}
	}
	Mtot =  0.0;
	for(itmp = 0;itmp < ntmp; itmp++){
		Mtot += Ptmp[itmp];
	}
	cout<<"# Overall probability for 4 dimensions  = "<<Mtot<<endl;
	file_out<<"# Overall probability for 4 dimensions  = "<<Mtot<<endl;
	file_out<<"# Overall Projeciton from "<<min<<" to "<<max<<endl;
	for(itmp = 0;itmp < ntmp; itmp++){
		Ptmp[itmp] = Ptmp[itmp]/Mtot;
		x = min + itmp * deltatmp;
		file_out<<x<<'\t'<<Ptmp[itmp]<<'\t'<<-kT*log(Ptmp[itmp])<<endl;
	}

	delete[] Ptmp;
	delete[] Pprjs;
	return 0;
}

int histogram::prj3D2D(ofstream &file_out, double kT){
	int dprj,p1,p2,d1,i,k,j,ind,indprj;
	double *A, *Mprj, Mtot, x, y;

	p1 = p2 = d1 = ind = -1;
	for(dprj=0;dprj<2;dprj++){
		//Allocate memory
		switch(dprj) {
			case 0:
				p1 = 0; p2 = 1; d1 = 2;
				break;
			case 1:
				p1 = 1; p2 = 2; d1 = 0;
				break;
		}
		Mprj = new double [ numbin[p1]*numbin[p2] ];
		A = new double [ numbin[p1]*numbin[p2] ];

		//Projecting data
		for(k = 0; k < numbin[p2]; k++){
			for(j = 0; j < numbin[p1]; j++){
				indprj = k*numbin[p1] + j;
				Mprj[indprj] = A[indprj] = 0.0; //initialize
				for(i = 0; i < numbin[d1]; i++){
					switch(dprj){
						case 0:
							ind = i*numbin[1]*numbin[0] +
								k*numbin[0] + j;
							break;
						case 1:
							ind = k*numbin[1]*numbin[0] +
								j*numbin[0] + i;
							break;
					}
					if(!isinf(P[ind])){
						Mprj[indprj] += P[ind];
					}
				}
			}
		}

		//Normalize
		Mtot = 0.0;
		indprj = 0;
		for(k = 0; k < numbin[p2]; k++){
			for(j = 0; j < numbin[p1]; j++){
				Mtot += Mprj[indprj];
				indprj++;
			}
		}
		cout<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		file_out<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		indprj = 0;
		for(k = 0; k < numbin[p2]; k++){
			for(j = 0; j < numbin[p1]; j++){
				Mprj[indprj] = Mprj[indprj]/Mtot;
				A[indprj] = -kT*log(Mprj[indprj]);
				indprj++;
			}
		}

		//Output
		indprj = 0;
		for(k = 0; k < numbin[p2]; k++){
			y = histmin[p2] + delta[p2]*k;
			for(j = 0; j < numbin[p1]; j++){
				x = histmin[p1] + delta[p1]*j;
				file_out<<x<<'\t'<<y<<'\t'<<Mprj[indprj]<<'\t'<<A[indprj]<<endl;
				indprj++;
			}
			file_out<<endl;
		}
		file_out<<endl<<endl<<endl;

		//deallocate
		delete[] A;
		delete[] Mprj;
	}
	return 0;
}

int histogram::prj3D1D(ofstream &file_out, double kT){
	int dprj,p1,d1,d2,i,k,j,ind,indprj;
	double *A, *Mprj, Mtot, x;

	p1 = d1 = d2 = ind = -1;
	for(dprj=0;dprj<3;dprj++){
		//Allocate memory
		switch(dprj) {
			case 0:
				p1 = 0; d1 = 1; d2 = 2;
				break;
			case 1:
				p1 = 1; d1 = 0; d2 = 2;
				break;
			case 2:
				p1 = 2; d1 = 0; d2 = 1;
				break;
		}
		Mprj = new double [ numbin[p1] ];
		A = new double [ numbin[p1] ];

		//Projecting data
		for(k = 0; k < numbin[p1]; k++){
			indprj = k;
			Mprj[indprj] = A[indprj] = 0.0; //initialize
			for(i = 0; i < numbin[d1]; i++){
				for(j = 0; j < numbin[d2]; j++){
					switch(dprj){
						case 0:
							ind = j*numbin[1]*numbin[0] +
								i*numbin[0] + k;
							break;
						case 1:
							ind = j*numbin[1]*numbin[0] +
								k*numbin[0] + i;
							break;
						case 2:
							ind = k*numbin[1]*numbin[0] +
								j*numbin[0] + i;
							break;
					}
					if(!isinf(P[ind])){
						Mprj[indprj] += P[ind];
					}
				}
			}
		}

		//Normalize
		Mtot = 0.0;
		indprj = 0;
		for(j = 0; j < numbin[p1]; j++){
			Mtot += Mprj[indprj];
			indprj++;
		}
		cout<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		file_out<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		indprj = 0;
		for(j = 0; j < numbin[p1]; j++){
			Mprj[indprj] = Mprj[indprj]/Mtot;
			A[indprj] = -kT*log(Mprj[indprj]);
			indprj++;
		}

		//Output
		indprj = 0;
		for(j = 0; j < numbin[p1]; j++){
			x = histmin[p1] + delta[p1]*j;
			file_out<<x<<'\t'<<Mprj[indprj]<<'\t'<<A[indprj]<<endl;
			indprj++;
		}
		file_out<<endl<<endl<<endl;

		//deallocate
		delete[] A;
		delete[] Mprj;
	}
	return 0;
}

int histogram::prj2D1D(ofstream &file_out, double kT){
	/* For the projection along the merged variable I assumed that the 
	 * two variables are identical (i.e. two ions along the same coordinate)
	 * In other situations this projection doesn't make sense*/
	int dprj,p1,d1,i,j,ind;
	double *A, *Mprj, Mtot, x;
	int numbin_merged,j_merged;
	double histmin_merged,histmax_merged,delta_merged,*P_merged,sum_merged;

	histmin_merged = histmin[0] < histmin[1] ? histmin[0] : histmin[1];
	histmax_merged = histmax[0] > histmax[1] ? histmax[0] : histmax[1];
	delta_merged = delta[0] < delta[1] ? delta[0] : delta[1];
	numbin_merged = int((histmax_merged - histmin_merged) / delta_merged) + 2;
	cout<<"# Histrogram for merged variable = "<<histmin_merged<<":"<<delta_merged<<":"<<histmax_merged<<" ["<<numbin_merged<<"] "<<endl;
	P_merged = new double [ numbin_merged ];
	for(j = 0; j < numbin_merged; j++) P_merged[j] = 0.0;
	for(dprj=0;dprj<2;dprj++){
		//Allocate memory
		switch(dprj) {
			case 0:
				p1 = 0; d1 = 1;
				break;
			case 1:
				p1 = 1; d1 = 0;
				break;
		}
		Mprj = new double [ numbin[p1] ];
		A = new double [ numbin[p1] ];

		//Projecting data
		for(j = 0; j < numbin[p1]; j++){
			Mprj[j] = A[j] = 0.0; //initialize
			for(i = 0; i < numbin[d1]; i++){
				switch(dprj){
					case 0:
						ind = i*numbin[0] + j;
						break;
					case 1:
						ind = j*numbin[0] + i;
						break;
				}
				if(!isinf(P[ind])){
					Mprj[j] += P[ind];
				}
			}
		}

		//Normalize
		Mtot = 0.0;
		for(j = 0; j < numbin[p1]; j++) Mtot += Mprj[j];
		cout<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		file_out<<"# Overall probability for projection along dimension "<<dprj<<" = "<<Mtot<<endl;
		for(j = 0; j < numbin[p1]; j++){
			Mprj[j] = Mprj[j]/Mtot;
			A[j] = -kT*log(Mprj[j]);
		}

		//Output
		for(j = 0; j < numbin[p1]; j++){
			x = histmin[p1] + delta[p1]*j;
			file_out<<x<<'\t'<<Mprj[j]<<'\t'<<A[j]<<endl;
			j_merged = (int)((x-histmin_merged+delta_merged/2.0)/delta_merged);
			//cout<<j_merged<<"\t"<<numbin_merged<<endl;
			P_merged[j_merged] += (Mprj[j]/delta[p1]);
		}
		file_out<<endl<<endl;

		//deallocate
		delete[] A;
		delete[] Mprj;
	}

	//Normalize merged variables
	sum_merged = 0.0;
	for(j = 0; j < numbin_merged; j++) sum_merged += P_merged[j];
	cout<<"# Overall probability along merged directions = "<<sum_merged<<endl;
	file_out<<"# Overall probability along merged directions = "<<sum_merged<<endl;
	//Output merged variables
	for(j = 0; j < numbin_merged; j++){
		P_merged[j] /= sum_merged;
		x = histmin_merged + delta_merged*j;
		file_out<<x<<'\t'<<P_merged[j]<<'\t'<<-kT*log(P_merged[j])<<endl;
	}
	delete[] P_merged;
	return 0;
}
//			ANALYSES END
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
//			DE-CONSTRUCTORS
histogram::~histogram(){
	if (histmin)   delete[] histmin;
	if (histmax)   delete[] histmax;
	if (refmin)    delete[] refmin;
	if (refmax)    delete[] refmax;
	if (hist)      delete[] hist;
	if (P)         delete[] P;
	if (Pold)      delete[] Pold;
	if (A)         delete[] A;
	if (delta)     delete[] delta;
	if (harmrest)  delete[] harmrest;
	if (center)    delete[] center;
	if (numwham)   delete[] numwham;
	if (denwham)   delete[] denwham;
	if (Punnorm)   delete[] Punnorm;
	if (Aold)      delete[] Aold;
	if (F)         delete[] F;
	if (g)         delete[] g;
	if (sumP)      delete[] sumP;
	if (numbin)    delete[] numbin;
	if (step)      delete[] step;
	if (whr_prd)   delete[] whr_prd;
	if (numbinwin) delete[] numbinwin;
	if (flaginf)   delete[] flaginf;
	if (periodic)  delete[] periodic;
	return;
}
//			END DE-CONTRUCTOR
////////////////////////////////////////////////////////////
