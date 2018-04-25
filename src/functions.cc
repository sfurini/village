////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
//PRINTHELP
//Output messaggi di help
int printhelp(char *prm_name){
if(strcmp(prm_name,"village")==0){
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                                   VILLAGE                                           "<<endl;
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                               INPUT FILE                                            "<<endl;
cout<<"-list       [village.list.txt]	It includes the list of the trajectory files        "<<endl;
cout<<"            The format of this file is:"<<endl;
cout<<"            FILE_NAME [TRAS] CENTER_HARMONIC_RESTRAINT K_HARMONIC [COL] [MASS] [EXTR]"<<endl;
cout<<"            If an harmonic potential acts on the center of mass of two CV, put them"<<endl;
cout<<"            after the other, and set the second harmonic constant to 0                  "<<endl;
cout<<"            COL = Columns were the CVs are in the trajectory files                      "<<endl;
cout<<"            Use this option if different files have CVs in different columns            "<<endl;
cout<<"            EXTR = There are the first and the last snapshot to read from the trajectory file"<<endl;
cout<<"-restart    Recover a previous calculations from this binary file"<<endl;
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                              OUTPUT FILE                                            "<<endl;
cout<<"-prefix     [village] Default prefix for output files"<<endl;
cout<<"-o          [village.out.txt] Textual output of the program"<<endl;
cout<<"-ene        [village.ene.txt] Energy"<<endl;
cout<<"-restartout [village.restart.dat] Binary file. It can be used for analyses/restart"<<endl;
cout<<"-histout    write histograms"<<endl;
cout<<"-histprj    write histogram projections"<<endl;
cout<<"-biasout    write harmonic potentials"<<endl;
cout<<"-discout    write discretized trajectoriesthe grid"<<endl;
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                              OTHER PARAMETERS                                       "<<endl;
cout<<"-h          Print this help"<<endl;
cout<<"-numdim     This must be the first parameter"<<endl;
cout<<"            #_DIMS MIN_DIM_1 MAX_DIM_1 STEP_DIM_1  [P #_P_DIM_1 / T T_DIM_1]"<<endl;
cout<<"            [MIN_DIM_i MAX_DIM_i STEP_DIM_i] [P #_P_DIM_1 / T #T_DIM_1]"<<endl;
cout<<"            P Force the same energy for the number of points #_P at the boundaries"<<endl;
cout<<"            T Force periodicity in energy along dimension i with period T_DIM_1"<<endl;
cout<<"-mass       [M_1 ..] Mass associated with the CVs"<<endl;
cout<<"            Use this for restrains on the centre of mass of two variables"<<endl;
cout<<"-col        [C_DIM_1 C_DIM_2 ...] Column were the CVs are in the trajectory files"<<endl;
cout<<"              if equal to 'listfile' read the column organization in the trajectory files"<<endl;
cout<<"              from the -list file"<<endl;
cout<<"-massFROMfile	Read masses from list file"<<endl;
cout<<"-trsFROMfile	Read translation from list file"<<endl;
cout<<"-debug      Execution with debug output (A LOT OF DEBUG OUTPUT !!)"<<endl;
cout<<"-tol        [1e-3] WHAM tollerance to stop iterations"<<endl;
cout<<"-numit      [0] WHAM maximum number of iterations"<<endl;
cout<<"-freqout    Number of step between successive output of restart/ene files"<<endl;
cout<<"-freqrmsd   Number of step between successive checks of RMSD"<<endl;
cout<<"-temperature [300] Kelvin"<<endl;
cout<<"-first      [0] First timestep to read from trajectory files (it checks the value on the first column of the traj file)"<<endl;
cout<<"-last       [0] Last timestep to read from trajectory files (it checks the value on the first column of the traj file)"<<endl;
cout<<"                0 = Read until the last one"<<endl;
cout<<"-stride     [1] IF cor == false, read the trajectory files with this step"<<endl;
cout<<"                IF cor == true, use this value as g parameter of the WHAM eqs"<<endl;
cout<<"-dt         [2e-3] Timestep in the traj files"<<endl;
cout<<"-extFROMfile	Read first/last/stride from list file"<<endl;
cout<<"-cor [false] IF true use the WHAM eqs. from correlated samples"<<endl;
cout<<"-extFROMfile [false] Read first/stride/last of the trajectory files"<<endl;
cout<<"             from the -list file"<<endl;
cout<<"-bootstrap  [1] Run this number of bootstrap iterations"<<endl;
cout<<"            1 = No bootstrap, normal WHAM"<<endl;
cout<<"-split      [1] Divide each trajectory into this number of fragments"<<endl;
cout<<"-reference  [] Define the region set to zero in the free energy"<<endl;
cout<<"            This is useful only for average calculations"<<endl;
}
else if(strcmp(prm_name,"project")==0){
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                                   PROJECT                                           "<<endl;
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                               INPUT FILE                                            "<<endl;
cout<<"-i	Read data this binary file"<<endl;
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                              OUTPUT FILE                                            "<<endl;
cout<<"-o	Textual output of the program"<<endl;
cout<<"-------------------------------------------------------------------------------------"<<endl;
cout<<"                              OTHER PARAMETERS                                       "<<endl;
cout<<"-h          Print this help"<<endl;
cout<<"-enemat	Output energy in MATLAB format"<<endl;
cout<<"-DXene	Output energy in DX format"<<endl;
cout<<"-DXhist	Output histograms in DX format"<<endl;
cout<<"-mep	Calculate MEP"<<endl;
cout<<"-mepstart	Starting point of the MEP"<<endl;
cout<<"-mepend	Ending point of the MEP"<<endl;
cout<<"-prjdim	Dimensionality of the projection"<<endl;
cout<<"-prjsamples	Project also the number of samples"<<endl;
cout<<"-nummin	Number of minima searched"<<endl;
cout<<"-distmin	Distance between minima"<<endl;
cout<<"-temperature	Temperature"<<endl;

} else {
	cerr<<"ERROR IN SUBROUTINE: int printhelp(char *prm_name)"<<endl;
	cerr<<"PROGRAM "<<prm_name<<" does not exist"<<endl;
	exit(1);
}
cout<<"-------------------------------------------------------------------------------------"<<endl;
  	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// HELIX FUNCTIONS
// Versor of a helical trajectory
int helix_versor(double &dx, double &dy, double &dz,double rad, double step,int clock, int dir, double teta, bool normalize = true){
	double norm;
	dx = - dir * rad * sin(teta); //There should be also omega, but it's a moltiplicative constant for all the components so it doesn't really matter for the versor
	dy = dir * rad * cos(teta);
	dz = dir * clock * (step/(2.0*PI));
	if(normalize){
		norm = sqrt(dx*dx + dy*dy + dz*dz);
		dx = dx/norm;
		dy = dy/norm;
		dz = dz/norm;
	}
	return 0;
}
// Radius and start angles of helical trajectory
int helix_parameters(double xhel, double yhel, double xcen, double ycen, double &rad, double &teta){
	rad = sqrt((xhel-xcen)*(xhel-xcen) + (yhel-ycen)*(yhel-ycen));
	teta = atan2(yhel-ycen,xhel-xcen);
	return 0;
}
// Length of the helix
double helix_length(double rad, double step, double numturn){
	return numturn*PI*sqrt(rad*rad + (step/(2.0*PI))*(step/(2.0*PI)));
}
double helix_length(double rad, double step, double xstart, double ystart, double zstart, double xcen, double ycen, double x, double y, double z){
	double tetastart,teta;
	tetastart = atan2(ystart-ycen,xstart-xcen);
	teta = atan2(y-ycen,x-xcen);
	if( (z - zstart) >= step ){
		teta += 2.0*PI;
		z -= step;
	}
	return (teta-tetastart)*sqrt(rad*rad + (step/(2.0*PI))*(step/(2.0*PI)));
}
