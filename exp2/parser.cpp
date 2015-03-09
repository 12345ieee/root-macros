//#include <TCanvas.h>
//#include <TLegend.h>
//#include <TRandom3.h>
//#include "TH1.h"
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
	
	hit(double x, double y): x(x), y(y) {}
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
	
	vector<hit> events1;
	vector<hit> events2;
	vector<hit> events3;
	
	double hv1;
	double hv2;
	double hv3;
	
	public:
	
	eeevent(int num1, int num_trig, int num_acq, int secs, int nsecs, int clk1, int clk2, double hv1, double hv2, double hv3):
		num1(num1), num_trig(num_trig), num_acq(num_acq), secs(secs), nsecs(nsecs),
		clk1(clk1), clk2(clk2), hv1(hv1), hv2(hv2), hv3(hv3) {}
	
	void gethit(double x, double y, int z) {
		hit newhit(x, y);
		switch (z) {
			case 149:
				this->events1.push_back(newhit);
				break;
			case 95:
				this->events2.push_back(newhit);
				break;
			case 49:
				this->events3.push_back(newhit);
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

void cluster_size_study(eeevector evector)
{
	
	
}

int main(int, char** argv)
{
	eeevector evector;
	evector = parse(argv[1]);
	
	cluster_size_study(evector);
	
	//compute_eff(1, evector);
	
	return 0;
}


/* Code graveyard
*
* 
* 
	hitvector hvector;
	switch (nchamber) {
		case 1: hvector=evector.events1; break;
		case 2: hvector=evector.events2; break;
		case 3: hvector=evector.events3; break;
		default: cout << "Unexpected chamber number" << endl;
	}
* 
* 
*/
