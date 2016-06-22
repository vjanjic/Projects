#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "ParticleFilter.h"
#include "TrackingModel.h"
#include "system.h"


using namespace std;

//templated conversion from string
template<class T>
bool num(const string s, T &n)
{	istringstream ss(s);
	ss >> n;
	return !ss.fail();
}

//write a given pose to a stream
inline void WritePose(ostream &f, vector<float> &pose)
{	for(int i = 0; i < (int)pose.size(); i++)
		f << pose[i] << " ";
	f << endl;
}

bool ProcessCmdLine(int argc, char **argv, string &path, int &cameras, int &frames, int &particles, int &layers, int &threads, int &threadModel, bool &OutputBMP)
{
	string    usage("Usage : Track (Dataset Path) (# of cameras) (# of frames to process)\n");
	usage += string("              (# of particles) (# of annealing layers) \n");
	usage += string("              [thread model] [# of threads] [write .bmp output (nonzero = yes)]\n\n");
	usage += string("        Thread model : 0 = Auto-select from available models\n");
        usage += string("                       1 = Intel TBB                 ");
        usage += string("(unavailable)\n");
        usage += string("                       2 = Posix / Windows threads   ");
        usage += string("(unavailable)\n");
        usage += string("                       3 = OpenMP                    ");
        usage += string("(unavailable)\n");
        usage += string("                       4 = Serial\n");

	string errmsg("Error : invalid argument - ");
	if(argc < 6 || argc > 9)															//check for valid number of arguments
	{	cout << "Error : Invalid number of arguments" << endl << usage << endl;
		return false;
	}
	path = string(argv[1]);																//get dataset path and add backslash if needed
	if(path[path.size() - 1] != DIR_SEPARATOR[0])
		path.push_back(DIR_SEPARATOR[0]);
	if(!num(string(argv[2]), cameras))													//parse each argument
	{	cout << errmsg << "number of cameras" << endl << usage << endl; 
		return false; 
	}
	if(!num(string(argv[3]), frames))													
	{	cout << errmsg << "number of frames" << endl << usage << endl; 
		return false; 
	}
	if(!num(string(argv[4]), particles))												
	{	cout << errmsg << "number of particles" << endl << usage << endl;
		return false;
	}
	if(!num(string(argv[5]), layers))													
	{	cout << errmsg << "number of annealing layers" << endl << usage << endl;
		return false;
	}
	threads = -1;
	threadModel = 0;
	if(argc < 7) 																		//use default single thread mode if no threading arguments present
		return true;
	if(!num(string(argv[6]), threadModel))
	{	cout << errmsg << "Thread Model" << endl << usage << endl;
		return false;
	}
	if(argc > 7)
		if(!num(string(argv[7]), threads))
		{	cout << errmsg << "number of threads" << endl << usage << endl;
			return false;
		}
	int n;
	OutputBMP = true;																	//do not output bmp results by default
	if(argc > 8)
	{	if(!num(string(argv[8]), n))
		{	cout << errmsg << "Output BMP flag" << endl << usage << endl;
			return false;
		}
		OutputBMP = (n != 0);
	}
	return true;
}

//Body tracking Single Threaded
int mainSingleThread(string path, int cameras, int frames, int particles, int layers, bool OutputBMP)
{
	cout << endl << "Running Single Threaded" << endl << endl;

	TrackingModel model;
	if(!model.Initialize(path, cameras, layers))										//Initialize model parameters
	{	cout << endl << "Error loading initialization data." << endl;
		return 0;
	}
	model.GetObservation(0);															//load data for first frame
	ParticleFilter<TrackingModel> pf;													//particle filter instantiated with body tracking model type
	pf.SetModel(model);																	//set the particle filter model
	pf.InitializeParticles(particles);													//generate initial set of particles and evaluate the log-likelihoods

	cout << "Using dataset : " << path << endl;
	cout << particles << " particles with " << layers << " annealing layers" << endl << endl;
	ofstream outputFileAvg((path + "poses.txt").c_str());

	vector<float> estimate;																//expected pose from particle distribution




	for(int i = 0; i < frames; i++)														//process each set of frames
	{	cout << "Processing frame " << i << endl;
		if(!pf.Update((float)i))														//Run particle filter step
		{	cout << "Error loading observation data" << endl;
			return 0;
		}		
		pf.Estimate(estimate);															//get average pose of the particle distribution
		WritePose(outputFileAvg, estimate);
		if(OutputBMP)
			pf.Model().OutputBMP(estimate, i);											//save output bitmap file
	}




	return 1;
}

int main(int argc, char **argv)
{
	string path;
	bool OutputBMP;
	int cameras, frames, particles, layers, threads, threadModel;								//process command line parameters to get path, cameras, and frames
  
  cout << "PARSEC Benchmark Suite" << endl << flush;

	if(!ProcessCmdLine(argc, argv, path, cameras, frames, particles, layers, threads, threadModel, OutputBMP))	
		return 0;
  if(threadModel == 0) {
    threadModel = 4;
  }
	switch(threadModel)
	{
		case 0 : 
      //This case should never happen, we auto-select the thread model before this switch
      cout << "Internal error. Aborting." << endl;
      exit(1);
			break;
  case 1 :
    cout << "Not compiled with Intel TBB support. " << endl;
    cout << "If the environment supports it, rebuild with USE_TBB #defined." << endl;
    break;  
  case 2 :
    cout << "Not compiled with Posix threads support. " << endl;
    cout << "If the environment supports it, rebuild with USE_THREADS #defined." << endl;
    break;
  case 3 : 
    cout << "Not compiled with OpenMP support. " << endl;
    cout << "If the environment supports OpenMP, rebuild with USE_OPENMP #defined." << endl;
    break;
  case 4 :
    mainSingleThread(path, cameras, frames, particles, layers, OutputBMP);                          //single threaded tracking
    break;
  default : 
    cout << "Invalid thread model argument. " << endl;
    cout << "Thread model : 0 = Auto-select thread model" << endl;
    cout << "               1 = Intel TBB" << endl;
    cout << "               2 = Posix / Windows Threads" << endl;
    cout << "               3 = OpenMP" << endl;
    cout << "               4 = Serial" << endl;
    break;
	}
	return 0;
}

