#include <TCanvas.h>
//#include <TLegend.h>
//#include <TRandom3.h>
#include "TH1.h"
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
	double x;
	double y;
	
	public:
	
	vector<hit> cluster;
	
	hit(double x, double y): x(x), y(y) {}
	
	double distance_from(hit h)
	{
		double Dx2 = pow((this->x-h.x), 2);
		double Dy2 = pow((this->y-h.y), 2);
		return sqrt(Dx2 + Dy2);
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
	
	vector<hit> hits1;
	vector<hit> hits2;
	vector<hit> hits3;
	
	eeevent(int num1, int num_trig, int num_acq, int secs, int nsecs, int clk1, int clk2, double hv1, double hv2, double hv3):
		num1(num1), num_trig(num_trig), num_acq(num_acq), secs(secs), nsecs(nsecs),
		clk1(clk1), clk2(clk2), hv1(hv1), hv2(hv2), hv3(hv3) {}
	
	void gethit(double x, double y, int z) {
		hit newhit(x, y);
		switch (z) {
			case 149:
				this->hits1.push_back(newhit);
				break;
			case 95:
				this->hits2.push_back(newhit);
				break;
			case 49:
				this->hits3.push_back(newhit);
				break;
			default:
				cout << "Unrecognized chamber" << endl;
		}
	}
};

typedef vector<eeevent> eeevector;
typedef vector<hit> hitvector;

vector<eeevent> parse(string path)
{
	cout << "Path given: " << path << endl;
	
	ifstream file;
	file.open(path);
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
	cout << "Parsing complete: " << evector.size() << " events captured" << endl;
	
	return evector;
}

void cluster_size_study(eeevector evector, double step, double maxradius)
{
	// Create the canvas
	TCanvas *canv = new TCanvas("canvas","canvas title",800,600);
	canv->cd();
	
	TH1D *hist = new TH1D("", "Hist title", 1000, 0, 100);
	
	double radius;
	int counter;
	
	for (radius=0; radius<maxradius; radius+=step) {          // for every radius
		counter=0;
		for (size_t i=0; i<evector.size(); ++i) {              // for every event
			eeevent event1 = evector[i];
			if (event1.hits1.size() > 1) {                      // if there are hits on ch1
				for (size_t j=0; j<event1.hits1.size(); ++j) {   // for every hit in ch1
					hit hit1 = event1.hits1[j];
					for (size_t k=j+1; k<event1.hits1.size(); ++k) {   // for other hit in ch1
						hit hit2 = event1.hits1[k];
						double dist = hit2.distance_from(hit1);
						if (dist < radius) {
							if (hit1.cluster.size() == 0 && hit1.cluster.size() == 0) counter++;
							hit1.cluster.push_back(hit2);
							hit2.cluster.push_back(hit1);
						}
					}
				}
			}
			if (event1.hits2.size() > 1) {                      // if there are hits on ch2
				for (size_t j=0; j<event1.hits2.size(); ++j) {   // for every hit in ch2
					hit hit1 = event1.hits2[j];
					for (size_t k=j+1; k<event1.hits2.size(); ++k) {   // for other hit in ch2
						hit hit2 = event1.hits2[k];
						double dist = hit2.distance_from(hit1);
						if (dist < radius) {
							if (hit1.cluster.size() == 0 && hit1.cluster.size() == 0) counter++;
							hit1.cluster.push_back(hit2);
							hit2.cluster.push_back(hit1);
						}
					}
				}
			}
			if (event1.hits3.size() > 1) {                      // if there are hits on ch3
				for (size_t j=0; j<event1.hits3.size(); ++j) {   // for every hit in ch3
					hit hit1 = event1.hits3[j];
					for (size_t k=j+1; k<event1.hits3.size(); ++k) {   // for other hit in ch3
						hit hit2 = event1.hits3[k];
						double dist = hit2.distance_from(hit1);
						if (dist < radius) {
							if (hit1.cluster.size() == 0 && hit1.cluster.size() == 0) counter++;
							hit1.cluster.push_back(hit2);
							hit2.cluster.push_back(hit1);
						}
					}
				}
			}
		}
		for (int i=0; i<counter; ++i) {
			hist->Fill(radius);
		}
	}
	hist->Draw();
	cout << "Done" << endl;
}

int parser(string path)
{
	const double step = 0.1;
	const double maxradius = 100;
	
	eeevector evector;
	evector = parse(path);
	
	cluster_size_study(evector, step, maxradius);
	
	//compute_eff(1, evector);
	
	return 0;
}

int main(int, char** argv)
{
	string path = argv[1];
	parser(path);
	
	return 0;
}

/* Code graveyard
*
* 
* 
	hitvector hvector;
	switch (nchamber) {
		case 1: hvector=evector.hits1; break;
		case 2: hvector=evector.hits2; break;
		case 3: hvector=evector.hits3; break;
		default: cout << "Unexpected chamber number" << endl;
	}
* 
* 
*/
