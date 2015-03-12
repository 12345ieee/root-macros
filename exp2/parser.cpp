#include <TCanvas.h>
//#include <TLegend.h>
//#include <TRandom3.h>
#include "TH1.h"
#include "TH2.h"
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

const double xmax=40;
const double ymax=84;

void plot_events(int ch, eeevector evector)
{
	char name[256];
	sprintf(name, "canv_ch%d", ch);
	
	char title[256];
	sprintf(title, "Events ch%d", ch);
	
	TCanvas *canv = new TCanvas(name,"canvas eff", 800, 600);
	TH2D *hist = new TH2D("", title, 100, -xmax, xmax, 1000, -ymax, ymax);
	
	
	for (size_t i=0; i<evector.size(); ++i) {              // for every event
		eeevent event1 = evector[i];
		for (size_t j=0; j<event1.hits[ch].size(); ++j) {   // for every hit in ch
			hit hit1 = event1.hits[ch][j];
			hist->Fill(hit1.x, hit1.y);
		}
	}
	
	canv->cd();
	hist->Draw("colz");
	
	cout << "Chamber " << ch << ": " << hist->Integral() << " events" << endl << endl;
}

void cluster_size_study(eeevector evector)
{
	const double step = 0.01;
	const double maxradius = 20;
	
	TCanvas *canv = new TCanvas("canvas","canvas title", 800, 600);
	canv->cd();
	TH1D *hist = new TH1D("", "Cluster", maxradius/step-1, 0, maxradius);
	
	for (size_t i=0; i<evector.size(); ++i) {              // for every event
		eeevent event1 = evector[i];
		for (int ch=0; ch<3; ++ch) {                        // for every chamber
			if (event1.hits[ch].size() > 1) {                      // if there are hits in ch
				for (size_t j=0; j<event1.hits[ch].size(); ++j) {   // for every hit in ch
					hit hit1 = event1.hits[ch][j];
					for (size_t k=j+1; k<event1.hits[ch].size(); ++k) {   // for every other hit in ch
						hit hit2 = event1.hits[ch][k];
						double dist = hit2.distance_from(hit1);
						for (double radius=maxradius; radius > dist; radius-=step) {   // for every radius
							hist->Fill(radius);
						}
					}
				}
			}
		}
	}
	hist->SetStats(kFALSE);
	hist->Draw();
	cout << "Done clusters" << endl << endl;
}

void chamber_assign(int chtest, int* ch1, int* ch2)
{
	switch (chtest) {
		case 0: *ch1 = 1; *ch2 = 2; break;
		case 1: *ch1 = 0; *ch2 = 2; break;
		case 2: *ch1 = 0; *ch2 = 1; break;
		default: cout << "Unrecognized chamber" << endl;
	}
}

hit compute_hitpoint(int chtest, hit hit1, hit hit2)
{
	double z, z1, z2;
	switch (chtest) {
		case 0: z = 149.00; z1 =  95.50; z2 =  49.00; break;
		case 1: z =  95.50; z1 = 149.00; z2 =  49.00; break;
		case 2: z =  49.00; z1 = 149.00; z2 =  95.50; break;
		default: cout << "Unrecognized chamber" << endl; exit(1);
	}
	// line equation is: (x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)=A
	
	double x1 = hit1.x;
	double x2 = hit2.x;
	double y1 = hit1.y;
	double y2 = hit2.y;
	
	double A = (z-z1)/(z2-z1);
	double x = (x2-x1)*A +x1;
	double y = (y2-y1)*A +y1;
	return {x, y};
} 

void compute_eff(int chtest, eeevector evector)
{
	const double dist_cutoff=3;
	
	char name_exp[256];
	sprintf(name_exp, "canv_exp_ch%d", chtest);
	char name_det[256];
	sprintf(name_det, "canv_det_ch%d", chtest);
	char name_eff[256];
	sprintf(name_eff, "canv_eff_ch%d", chtest);
	
	char title_exp[256];
	sprintf(title_exp, "Expected ch%d", chtest);
	char title_det[256];
	sprintf(title_det, "Detected ch%d", chtest);
	char title_eff[256];
	sprintf(title_eff, "Efficiency ch%d", chtest);
	
	TCanvas *canv_exp = new TCanvas(name_exp,"canvas exp", 800, 600);
	TH2D *hist_exp = new TH2D("", title_exp, 100, -xmax, xmax, 1000, -ymax, ymax);
	TCanvas *canv_det = new TCanvas(name_det,"canvas det", 800, 600);
	TH2D *hist_det = new TH2D("", title_det, 100, -xmax, xmax, 1000, -ymax, ymax);
	TCanvas *canv_eff = new TCanvas(name_eff,"canvas eff", 800, 600);
	TH2D *hist_eff = new TH2D("", title_eff, 100, -xmax, xmax, 1000, -ymax, ymax);
	
	int expected=0;
	int detected=0;
	
	int ch1, ch2;
	chamber_assign(chtest, &ch1, &ch2);
	
	char name_nhits[256];
	sprintf(name_nhits, "canv_nhits_ch%d", chtest);
	char title_nhits[256];
	sprintf(title_eff, "Hits number ch%d", chtest);
	
	TCanvas *canv_nhits = new TCanvas(name_nhits,"canvas nits", 800, 600);
	TH1D *hist_nhits = new TH1D("", title_nhits, 10, 0, 10);
	
	for (size_t i=0; i<evector.size(); ++i) {              // for every event
		eeevent event1 = evector[i];
		for (size_t j=0; j<event1.hits[ch1].size(); ++j) {   // for every hit in ch1
			hit hit1 = event1.hits[ch1][j];
			for (size_t k=0; k<event1.hits[ch2].size(); ++k) {   // for every hit in ch2
				hit hit2 = event1.hits[ch2][k];
				hit hit_test = compute_hitpoint(chtest, hit1, hit2);
				if (hit_test.is_inside(xmax, ymax)) {
					int counter=0;
					expected++;
					hist_exp->Fill(hit_test.x, hit_test.y);
					for (size_t h=0; h<event1.hits[chtest].size(); ++h) {   // for every hit in chtest
						hit hitf = event1.hits[chtest][h];
						double dist = hit_test.distance_from(hitf);
						if (dist < dist_cutoff) {
							detected++;
							hist_det->Fill(hitf.x, hitf.y);
							counter++;
						}
					}
					hist_nhits->Fill(counter);
				}
			}
		}
	}
	
	canv_nhits->cd();
	hist_nhits->Draw();
	
	canv_exp->cd();
	hist_exp->Draw();
	
	canv_det->cd();
	hist_det->Draw();
	
	hist_eff->Divide(hist_det, hist_exp);
	canv_eff->cd();
	hist_eff->Draw("colz");
	
	cout << "Chamber " << chtest << ":" << endl;
	cout << "Distance cutoff chosen: " << dist_cutoff << endl;
	cout << "Expected events: " << expected << endl;
	cout << "Detected events: " << detected << endl;
	cout << "Efficiency: " << (double)detected/expected << endl << endl;
}

void single_events_fitter(eeevector evector)
{
	eeevector goodevents[3];  // events good for chamber ch
	eeevector perfevents;     // events with a single hit/ch
	for (size_t i=0; i != evector.size(); ++i) {  // for every event
		int counter=0;
		eeevent event1=evector[i];
		if (event1.hits[0].size()==1) counter+=1;
		if (event1.hits[1].size()==1) counter+=2;
		if (event1.hits[2].size()==1) counter+=4;
		if (counter==3 || counter==7) goodevents[2].push_back(event1);
		if (counter==5 || counter==7) goodevents[1].push_back(event1);
		if (counter==6 || counter==7) goodevents[0].push_back(event1);
		if (counter==7) perfevents.push_back(event1);
		if (counter>7) cout << "Unrecognized counter" << endl;
	}
	for (int ch=0; ch<3; ++ch) {
		cout << "Good events for ch" << ch << ": " << goodevents[ch].size() << endl;
	}
	cout << "Perfect events: " << perfevents.size() << endl << endl;
	
	
	//~ for (size_t i=0; i<perfevents.size(); ++i) {              // for every event
		//~ eeevent event1 = perfrvents[i];
		//~ hit hit0 = event1.hits[0][0];
		//~ hit hit1 = event1.hits[1][0];
		//~ hit hit2 = event1.hits[2][0];
		//~ hit hit_test = compute_hitpoint(chtest, hit1, hit2);
		//~ if (hit_test.is_inside(xmax, ymax)) {
			//~ int counter=0;
			//~ expected++;
			//~ hist_exp->Fill(hit_test.x, hit_test.y);
			//~ for (size_t h=0; h<event1.hits[chtest].size(); ++h) {   // for every hit in chtest
				//~ hit hitf = event1.hits[chtest][h];
				//~ double dist = hit_test.distance_from(hitf);
				//~ if (dist < dist_cutoff) {
					//~ detected++;
					//~ hist_det->Fill(hitf.x, hitf.y);
					//~ counter++;
				//~ }
			//~ }
			//~ hist_nhits->Fill(counter);
		//~ }
	//~ }
}

void average_hit_number(eeevector evector)
{
	for (int ch=0; ch<3; ++ch) {
		char name[256];
		sprintf(name, "canv_hits_ch%d", ch);
		
		char title[256];
		sprintf(title, "Hit number ch%d", ch);
		
		TCanvas *canv = new TCanvas(name,"canvas", 800, 600);
		int nbins=20;
		TH1I *hist = new TH1I("", title, nbins, 0, nbins);
		
		for (size_t i=0; i != evector.size(); ++i) {  // for every event
			int nhits=evector[i].hits[ch].size();
			hist->Fill(nhits);
		}
		
		canv->cd();
		hist->Draw();
	}
}

int parser(string path)
{
	
	eeevector evector;
	evector = parse(path);
	
	//for (int i = 0; i < 3; i++) plot_events(i, evector);
	 
	//cluster_size_study(evector);
	//evector = clusterize(evector, cluster_size);
	 
	//for (int i = 0; i < 3; i++) compute_eff(i, evector);
	
	average_hit_number(evector);
	
	single_events_fitter(evector);
	
	
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
	if (dist < radius) {
		if (hit1.cluster.size() == 0 && hit1.cluster.size() == 0) counter++;
		hit1.cluster.push_back(hit2);
		hit2.cluster.push_back(hit1);
	}
* 
* 
*/
