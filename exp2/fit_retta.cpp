#include <TCanvas.h>
//#include <TLegend.h>
//#include <TRandom3.h>
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
//#include "TF1.h"
//#include "TMath.h"
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class hit
{
	public:
	
	double x;
	double y;
	
	//vector<hit> cluster;
	
	hit(double x, double y): x(x), y(y) {}
	
	double distance_from(hit h)
	{
		double Dx2 = pow((this->x - h.x), 2);
		double Dy2 = pow((this->y - h.y), 2);
		return sqrt(Dx2 + Dy2);
	}
	
	bool is_inside(double xmax, double ymax) {
		return (abs(x)<=xmax && abs(y)<=ymax);
	}
};

class eeevent
{
	int num1;
	
	int num_trig;
	int num_acq;
	
	int secs;
	int nsecs;
	
	int clk1;
	int clk2;
	
	double hv1;
	double hv2;
	double hv3;
	
	public:
	
	vector<hit> hits [3];  // array of vector of hit
	
	eeevent(int num1, int num_trig, int num_acq, int secs, int nsecs, int clk1, int clk2, double hv1, double hv2, double hv3):
		num1(num1), num_trig(num_trig), num_acq(num_acq), secs(secs), nsecs(nsecs),
		clk1(clk1), clk2(clk2), hv1(hv1), hv2(hv2), hv3(hv3) {}
	
	void gethit(double x, double y, int z) {
		hit newhit(x, y);
		switch (z) {
			case 149:
				this->hits[0].push_back(newhit);
				break;
			case 95:
				this->hits[1].push_back(newhit);
				break;
			case 49:
				this->hits[2].push_back(newhit);
				break;
			default:
				cout << "Unrecognized chamber" << endl;
		}
	}
};

typedef vector<eeevent> eeevector;
typedef vector<hit> hitvector;

eeevector parse(string path)
{
	cout << "Path given: " << path << endl;
	
	ifstream file;
	file.open(path);
	if (file==NULL) {
		cout << "File does not exist" << endl;
		exit(1);
	}
	const int linesize=65536;
	char* line = (char*)malloc(linesize*sizeof(char));  // get buffer for every line
	
	vector<eeevent> evector;               // vector for the events
	double hv1, hv2, hv3;                  // variables for high voltage values
	
	while (!file.eof()) {
		file.getline(line, linesize);       // fill line
		
		char* piece = strtok (line," \t");  // get initial number
		if (piece==NULL) continue;          // the file is finished at this point, but better run the eof() check
		int num1 = atoi(piece);             // inital number
		if (num1 == 0) continue;            // line needs to start with a number, not characters (like the name line)
		
		piece = strtok (NULL, " \t");       // get descriptor
		string descriptor = piece;
		
		if (descriptor=="STATUS") {
			char* next=piece+7;              // hack to get after "STATUS"
			sscanf(next, "HV %lf %*f %lf %*f %lf %*f", &hv1, &hv2, &hv3); 
			//cout << "hv: " << hv1 << " - " << hv2 << " - " << hv3 << endl;
		}
		else if (descriptor=="EVENT") {
			
			piece = strtok (NULL, " \t");    // get num_trig
			int num_trig = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get num_acq
			int num_acq = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get seconds
			int secs = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get nseconds
			int nsecs = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get clk1
			int clk1 = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get clk2
			int clk2 = atoi(piece);
			
			eeevent neweeevent(num1, num_trig, num_acq, secs, nsecs, clk1, clk2, hv1, hv2, hv3);
			
			piece = strtok (NULL, " \t");    // get events
			while (piece != NULL) {
				double x = atof(piece);
				piece = strtok (NULL, " \t");    // get more pieces
				
				if (piece!=NULL) {
					double y = atof(piece);
					piece = strtok (NULL, " \t");    // get more pieces
					
					int z = atoi(piece);
					piece = strtok (NULL, " \t");    // get more pieces
					
					neweeevent.gethit(x, y, z);
					//cout << x << " - " << y << " - " << z << endl;
				}
			}
			
			evector.push_back(neweeevent);
		}
	}
	cout << "Parsing completed: " << evector.size() << " events captured" << endl << endl;
	
	return evector;
}

void hits_from_mc(hitvector hits[3], int N)
{
	double x1[N], x2[N], x3[N], y1[N], y2[N], y3[N];
	// FIXME: io li definisco qui, ma tu dovresti già avere array di eventi generati col MC
	// devi anche quantizzarli in x e sporcarli in y prima di darli in pasto a questo programma
	for (int ev=0; ev<N; ++ev) {
		hit hit0(x1[ev], y1[ev]);
		hit hit1(x2[ev], y2[ev]);
		hit hit2(x3[ev], y3[ev]);
		hits[0].push_back(hit0);
		hits[1].push_back(hit1);
		hits[2].push_back(hit2);
	}
}

const double xmin= -46.5;
const double xmax= 43.63;
const double xbin = (xmax-xmin)/24;
const double ymax= 84;

int parser(string path)
{
	
	//eeevector evector;
	//evector = parse(path);
	
	const int N=10;  // FIXME: numero di eventi generati;
	
	vector<hit> hits[3];  // will be filled with hits
	hits_from_mc(hits, N);
	
	for (int ev=0; ev<N; ++ev) {  // for every event
		
	}
	
	return 0;
}

int main(int, char** argv)
{
	string path = argv[1];
	parser(path);
	
	return 0;
}
