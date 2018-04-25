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

int WriteEne(ofstream &file_out, double kT, int numdim, int numwin, double *histmin, double *histmax, int *numbin,
		double *P, bool *flaginf, double *numwham, double *denwham){
	int numhist = numbin[0];
	int idim,ihist,ipos;
	double dtmp,A;
	double Pnorm[numhist*numwin];
	double delta[numdim];
	int step[numdim];

	//Initialize
	for(idim = 0;idim< numdim; idim++){
		delta[idim] = (histmax[idim] - histmin[idim])/(numbin[idim]-1);
		if(idim == 0) step[idim] = 1;
		else step[idim] = step[idim-1]*numbin[idim-1];
	}
	for(idim = 1; idim < numdim; idim++) numhist = numhist*numbin[idim];

	//Output Heading
	file_out<<"# numdim = "<<numdim<<endl;
	for(idim = 0;idim < numdim; idim++){
		file_out<<"# dim = "<<idim
			<<" min = "<<histmin[idim]
			<<" max = "<<histmax[idim]
			<<" numbin = "<<numbin[idim]<<endl;

	}
	file_out<<"# numwin = "<<numwin<<endl;
	file_out<<"# Coordinates... P[ihist] ihist numwham[ihist] denwham[ihist] kT*A[ihist] "<<endl;
	//Normalize probabilities
	dtmp = 0.0;
	for(ihist = 0; ihist < numhist; ihist++)if(!flaginf[ihist])dtmp += P[ihist];
	for(ihist = 0; ihist < numhist; ihist++){
		if(flaginf[ihist])Pnorm[ihist]=0.0;
		else Pnorm[ihist] = P[ihist]/dtmp;
	}
	//Output Data
	for(ihist = 0; ihist < numhist; ihist++){
		A = -kT*log(Pnorm[ihist]);
		for(idim = 0;idim < numdim; idim++){
			ipos = ((int)(ihist/step[idim])) % numbin[idim];
			file_out<<histmin[idim] + ipos*delta[idim]<<'\t';
		}
		file_out<<Pnorm[ihist]<<'\t'<<ihist<<'\t'<<numwham[ihist]<<'\t'<<denwham[ihist]<<'\t'<<A<<endl;
		if((ihist!=0) && (((ihist+1)%numbin[0]) == 0))file_out<<endl;
	}

	return 0;
}

void start_3Dout(ofstream &file_out, double* min, int *n, double *d){
	file_out<<setiosflags(ios::fixed);
	file_out<<setprecision(3);
	file_out<<"object 1 class gridpositions counts "<<n[2]<<" "<<n[1]<<" "<<n[0]<<endl;
	file_out<<"origin "<<min[0]<<" "<<min[1]<<" "<<min[2]<<endl;
	file_out<<"delta 0.000000e+00 0.000000e+00 "<<d[2]<<endl;
	file_out<<"delta 0.000000e+00 "<<d[1]<<" 0.000000e+00"<<endl;
	file_out<<"delta "<<d[0]<<" 0.000000e+00 0.000000e+00"<<endl;
	file_out<<"object 2 class gridconnections counts "<<n[2]<<" "<<n[1]<<" "<<n[0]<<endl;
	file_out<<"object 3 class array type double rank 0 items "<<n[0]*n[1]*n[2]<<" data follows"<<endl;
	return;
}

void end_3Dout(ofstream &file_out){
	file_out<<endl;
	file_out<<"attribute \"dep\" string \"positions\""<<endl;
	file_out<<"object \"regular positions regular connections\" class field"<<endl;
	file_out<<"component \"positions\" value 1"<<endl;
	file_out<<"component \"connections\" value 2"<<endl;
	file_out<<"component \"data\" value 3"<<endl;
	file_out<<endl;
	return;
}

int readDATAMATRIX(ifstream &file_in, int &numdim, double *min, double *max, int *n, vector<double> &M){
	char line[10001];
	string srtmp;
	double dtmp;
	int nread,idim,npoints;

	//Read file heading
	npoints = 1;
	file_in.getline(line,400);
	nread = sscanf(line,"# numdim = %d\n",&numdim);
	if(nread != 1)return 1;
	for(idim=0;idim<numdim;idim++){
		file_in.getline(line,400);
		nread = sscanf(line,"# dim = %*d min = %lf max = %lf numbin = %d\n",min+idim,max+idim,n+idim);
		if(nread != 3){
			cerr<<"ERROR"<<endl<<line<<endl;
			cerr<<"nread = "<<nread<<endl;
			return 1;
		}
		npoints = npoints * n[idim];
	}

	file_in.getline(line,10000);
	srtmp = line;
	while(!file_in.eof()){
		if((srtmp.length()!=0)&&(srtmp[0]!='#')){
			for(idim=0;idim<numdim;idim++){
				srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));
				srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
				if(srtmp.length()==0){
					cerr<<line<<endl;
					return 1;
				}
			}
			while(srtmp.length()!=0){
				if ( srtmp.substr(0,srtmp.find_first_of(" \t")).find("inf") == string::npos ){
					dtmp = atof(srtmp.substr(0,srtmp.find_first_of(" \t")).c_str());
				}
				else{
					if ( srtmp.substr(0,srtmp.find_first_of(" \t")).find("-") == string::npos ){
						dtmp = 1e5;
					}else{
						dtmp = -1e5;
					}
					//cerr<<"WARNING: infinite value"<<endl;
				}
				srtmp.erase(0,srtmp.find_first_of(" \t\n\0"));//delete read
				srtmp.erase(0,srtmp.find_first_not_of(" \t\n\0"));
			}
			M.push_back(dtmp);
		}
		file_in.getline(line,10000);
		srtmp = line;
	}

	if((int)M.size() != npoints)return 1;
	else return 0;
}

int writeDATAMATRIX(ofstream &file_out, int &numdim, double *histmin, double *histmax, int *numbin, double *M){
	double delta[numdim],pos;
	int step[numdim],numhist,idim,ihist,ipos;

	//Initialize delta & step
	for(idim = 0; idim < numdim; idim++){
		delta[idim] = (histmax[idim] - histmin[idim])/(numbin[idim]-1);
		if(idim == 0) step[idim] = 1;
		else step[idim] = step[idim-1]*numbin[idim-1];
	}
	numhist = numbin[0];
	for(idim = 1; idim < numdim; idim++) numhist = numhist*numbin[idim];

	//Write file heading
	file_out<<"# numdim = "<<numdim<<endl;
	for(idim = 0; idim < numdim; idim++){
		file_out<<"# dim = "<<idim<<" min = "<<histmin[idim]<<" max = "<<histmax[idim]<<" numbin = "<<numbin[idim]<<endl;
	}

	//Write data
	for(ihist = 0; ihist< numhist; ihist++){
		for(idim = 0;idim < numdim; idim++){
			ipos = ((int)(ihist/step[idim])) % numbin[idim];
			pos = histmin[idim]+ipos*delta[idim];
			file_out<<pos<<'\t';
		}
		file_out<<M[ihist]<<endl;
	}

	return 0;
}

