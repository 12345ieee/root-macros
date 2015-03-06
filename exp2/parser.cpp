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

class eeevent
{
	int num1;
	
	int num_trig;
	int num_acq;
	
	int secs;
	int nsecs;
	
	int clk1;
	int clk2;
	
	vector<int[2]> event1;
	vector<int[2]> event2;
	vector<int[2]> event3;
	
	double hv1;
	double hv2;
	double hv3;
	
	eeevent(int num1, int num_trig, int num_acq, int secs, int nsecs, int clk1, int clk2, int hv1, int hv2, int hv3):
		num1(num1), num_trig(num_trig), num_acq(num_acq), secs(secs), nsecs(nsecs),
		clk1(clk1), clk2(clk2), hv1(hv1), hv2(hv2), hv3(hv3) {}
};

int parse(string path)
{
	cout << "Path given: " << path << endl;
	
	ifstream file;
	file.open(path);
	const int linesize=65536;
	char* line = (char*)malloc(linesize*sizeof(char));  // get buffer for every line
	
	vector<eeevent> evector;
	
	while (!file.eof()) {
		file.getline(line, linesize);       // fill line
		
		char* piece = strtok (line," \t");  // get initial number
		if (piece==NULL) continue;          // the file is finished at this point, but better run the eof() check
		int num1 = atoi(piece);             // inital number
		if (num1 == 0) continue;            // line needs to start with a number, not characters (like the name line)
		
		piece = strtok (NULL, " \t");       // get descriptor
		string descriptor = piece;
		if (descriptor=="STATUS") {
			
		}
		else if (descriptor=="EVENT") {
			
		}
		else continue;
		
		while (piece != NULL) {
			piece = strtok (NULL, " \t");    // get more pieces
			//cout << piece << endl;
		}
	}

	
	return 0;
}

int main(int, char** argv)
{
	parse(argv[1]);
	return 0;
}
