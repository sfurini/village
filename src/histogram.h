#ifndef HISTOGRAM
#define HISTOGRAM

class histogram {
	private:
		// numwham = cumulative number of points in each bin
		// denwham = sum(number_of_points_i*exp(-beta*(Wi-Fi)))
		// F = Energy reference for each window
		// numbinwin = cumulative number of point in each window
		// hist = number of points for each window in each bin
		int numwin,numhist;
		double *histmin,*histmax,*g,*refmin,*refmax,*delta,*harmrest,*center;
		int *numbin,*numbinwin,*step,*hist;
		double *P,*Punnorm,*Pold,*A,*Aold,*numwham,*denwham,*F,*sumP;
		bool *flaginf;
		vector < vector <double> > backup_energy_reiter;

		int prj4D2D(ofstream &file_out, bool collect_samples = false,double kT = 1.0);
		int prj4D1D(ofstream &file_out, double kT = 1.0);
		int prj3D2D(ofstream &file_out, double kT = 1.0);
		int prj3D1D(ofstream &file_out, double kT = 1.0);
		int prj2D1D(ofstream &file_out, double kT = 1.0);

		//Analyses
		int defineMEPmatrix4(vector<double> &A,ofstream &file_out, char cmt
				,unsigned long int ihist_str, unsigned long int ihist_end
				,bool flag_forward);
		int FindMEP4(ofstream &file_out, vector<double> pnt_str, vector<double> pnt_end, char cmt = '#'
				,double kT = 1.0,int nummin = 100, double dist_min = 2.0);
		int defineMEPmatrix3(vector<double> &A,ofstream &file_out, char cmt
				,unsigned long int ihist_str, unsigned long int ihist_end
				,bool flag_forward);
		int FindMEP3(ofstream &file_out, vector<double> pnt_str, vector<double> pnt_end, char cmt = '#'
				,double kT = 1.0,int nummin = 100, double dist_min = 2.0);
		int defineMEPmatrix2(vector<double> &A,ofstream &file_out, char cmt
				,unsigned long int ihist_str, unsigned long int ihist_end
				,bool flag_forward);
		int FindMEP2(ofstream &file_out, vector<double> pnt_str, vector<double> pnt_end, char cmt = '#'
				,double kT = 1.0,int nummin = 100, double dist_min = 2.0);
		double CalculateDeviationProb(int ihist,int iwin); // Formula Furini-Barbini
		double CalculateDeviationProb(int ihist); // Formula Ferrenber-Swendesen
	public:
		double *whr_prd;
		vector<double> mass;
		int numdim,numbad,numgood,*periodic;
		bool debug;

		//Contructor
		histogram();

		//Input
		int readrestart(ifstream &file_in);
		int AllocateDimensions();
		int InitializeDimensions(); 
		int AllocateWindows();
		int DefineBoundaries(int idim, double min, double max, int num);
		int DefineReferences(int idim, double min, double max);
		int readlist(ifstream &file_in,ofstream &file_out, vector<int> col, int first = 0, int last = 0, bool cor = false, int stride = 1, bool trsFROMfile = false, bool massFROMfile = false, bool colFROMfile = false, bool extFROMfile = false, int ind_reiter = 1, int num_split = 1, int num_bootstrap = 1, double kT = 1.0, bool discout = false, char* prefix = NULL, double timestep_traj = 0.002);

		//Wham
		double bias(int ihist, int iwin);
		int InitializeWham(ofstream &file_out);
		int NewProbabilities(ofstream &file_out, int it = 0);
		int NewConstants(ofstream &file_out, int it = 0);
		int ForcePeriodicity(ofstream &file_out, int it = 0);
		double CheckConvergence(ofstream &file_out, int it = 0,double kT = 1.0);
		double NormalizeProbability();
		int ComputeEnergy(double kT = 1.0);
		int TranslateEnergy(double kT = 1.0);

		//Analyses
		int FindMEP(char *fileout, vector<double> pnt_str, vector<double> pnt_end,char cmt = '#',double kT = 1.0,int nummin = 100, double dist_min = 2.0);

		//Output
		int WriteHeading(ofstream &file_out);
		int WriteRestart(ofstream &file_out);
		int MemorizeEne(double kT = 1.0);
		int WriteEne(char *fileene, int ind_reiter = 0, double kT = 1.0);
		int WriteEne(ofstream &file_out,int num_reiter);
		int WriteEneMat(char *fileout, double kT = 1.0);
		int WriteHist(ofstream &file_out, bool histprj = false, double kT = 1.0);
		int WriteBias(ofstream &file_out, double kT = 1.0);
		int WriteProjection(char *fileout, int prjdim, bool collect_samples = false, double kT = 1.0);
		int WriteDxHist(char *fileout);
		int WriteDxEne(char *fileout, double kT = 1.0);

		//De-Constructor
		~histogram();
};

#endif
