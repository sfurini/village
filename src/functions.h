// General Routines
int printhelp(char *prm_name);
int helix_versor(double &dx, double &dy, double &dz,double rad
		, double step,int clock, int dir, double teta, bool normalize = true);
int helix_parameters(double xhel, double yhel, double xcen, double ycen, double &rad, double &teta);
double helix_length(double rad, double step, double numturn);
double helix_length(double rad, double step, double xstart, double ystart, double zstart, double xcen, double ycen, double x, double y, double z);

//readwrite.cc
int WriteEne(ofstream &file_out, double kT, int numdim, int numwin, double *histmin, double *histmax, int *numbin,
		double *P, bool *flaginf, double *numwham, double *denwham);
void start_3Dout(ofstream &file_out, double* min, int *n, double *d);
void end_3Dout(ofstream &file_out);
int readDATAMATRIX(ifstream &file_in, int &numdim, double *min, double *max, int *n, vector<double> &M);
int writeDATAMATRIX(ofstream &file_out, int &numdim, double *histmin, double *histmax, int *numbin, double *M);
