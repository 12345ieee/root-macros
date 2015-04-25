#include <TCanvas.h>
//#include <TLegend.h>
//#include <TRandom3.h>
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "TSystemDirectory.h"

const int chnum = 3;

const double xmin_scint=   5;
const double xmax_scint=  18;
const double ymin_scint= -45;
const double ymax_scint= +45;

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
	bool is_inside_scint() {
		return (x>=xmin_scint && x<=xmax_scint && y >= ymin_scint && y<=ymax_scint);
	}
};
// eee
// calib

const int z0t = 145;
const int z1t = 100;
const int z2t = 23;
const int z0  = 149;
const int z1  = 95;
const int z2  = 49;
 
const double xmin = -0.5;
const double xmax = 23.5;

const double ymax = 84;

const int dx=1;
const double dy[chnum]={1.5,1.6,1.5};
const double dy_corr[chnum]={7.8, 3.8, 6.6}; // 2 sigma

const double yacc =0;  // to include the offset
const double yres = (dy[0]+dy[1]+dy[2])/3;
const double ydraw= yres/2*4; // to account for resolution

class eeevent
{
	int num1;
	
	int num_trig;
	int num_acq;
	
	public:
	
	int secs;
	int nsecs;
	
	int clk1;
	int clk2;
	
	double hv0;
	double hv1;
	double hv2;
	
	vector<hit> hits [3];  // array of vector of hit
	
	eeevent(int num1, int num_trig, int num_acq, int secs, int nsecs, int clk1, int clk2, double hv0, double hv1, double hv2):
		num1(num1), num_trig(num_trig), num_acq(num_acq), secs(secs), nsecs(nsecs),
		clk1(clk1), clk2(clk2), hv0(hv0), hv1(hv1), hv2(hv2) {}
	
	void gethit(double x, double y, int z) {
		hit newhit(x, y);
		switch (z) {
			case z0:
			case z0t:
				this->hits[0].push_back(newhit);
				break;
			case z1:
			case z1t:
				this->hits[1].push_back(newhit);
				break;
			case z2:
			case z2t:
				if (x!=18)  // kill dead channel on calib data
				this->hits[2].push_back(newhit);
				break;
			default:
				cout << "Unrecognized chamber" << endl;
		}
	}
	
	double time() {
		return secs+1e-9*nsecs;
	}
	double delta_time_ms(eeevent& ev1) {
		int dsecs = this->secs - ev1.secs;
		int dnsecs = this->nsecs - ev1.nsecs;
		return dsecs*1000+1e-6*dnsecs;
	}
};

typedef vector<eeevent> eeevector;
typedef vector<hit> hitvector;

void correct_coords (double* x, double* y, string mode) 
{
	double xmin = -60;
	double xstep =  5;
	double yratio = 0.64;
	
	if (mode=="calib") {
		xmin = -46.5;
		xstep=  3.875;
		yratio*= 15./17.;
	}
	*x = round(((*x)-xmin)/xstep);
	*y = (*y)*yratio;
}

void parse_file(string path, eeevector& evector, string mode)
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
	
	double hv0, hv1, hv2;                  // variables for high voltage values
	
	double dsecs=0;
	double dnsecs=0;
	if (evector.size()>0) {
		dsecs=evector.back().secs;
		dnsecs=evector.back().nsecs;
	}
	
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
			sscanf(next, "HV %lf %*f %lf %*f %lf %*f", &hv0, &hv1, &hv2); 
			//cout << "hv: " << hv0 << " - " << hv1 << " - " << hv2 << endl;
		}
		else if (descriptor=="EVENT") {
			
			piece = strtok (NULL, " \t");    // get num_trig
			int num_trig = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get num_acq
			int num_acq = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get seconds
			int secs = dsecs + atoi(piece);
			
			piece = strtok (NULL, " \t");    // get nseconds
			int nsecs = dnsecs + atoi(piece);
			
			secs+= nsecs/(int)1e9;
			nsecs= nsecs%(int)1e9;
			
			piece = strtok (NULL, " \t");    // get clk1
			int clk1 = atoi(piece);
			
			piece = strtok (NULL, " \t");    // get clk2
			int clk2 = atoi(piece);
			
			eeevent neweeevent(num1, num_trig, num_acq, secs, nsecs, clk1, clk2, hv0, hv1, hv2);
			
			piece = strtok (NULL, " \t");    // get events
			while (piece != NULL) {
				double x = atof(piece);
				piece = strtok (NULL, " \t");    // get more pieces
				
				if (piece!=NULL) {
					double y = atof(piece);
					piece = strtok (NULL, " \t");    // get more pieces
					
					int z = atoi(piece);
					piece = strtok (NULL, " \t");    // get more pieces
					
					correct_coords(&x, &y, mode);
					
					neweeevent.gethit(x, y, z);
					//cout << x << " - " << y << " - " << z << endl;
				}
			}
			
			evector.push_back(neweeevent);
		}
	}
	cout << "Parsing completed: " << evector.size() << " events captured" << endl << endl;
}

void parse_dir(const char* dirname, eeevector& evector, string mode)
{
	TSystemDirectory dir(dirname, dirname);
	string dirnames(dirname);
	TList *files = dir.GetListOfFiles();
	files->Sort();
	if (files) {
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.BeginsWith("PISA")) {
				stringstream s;
				s << fname;
				parse_file(dirnames + s.str(), evector, mode);
			}
		}
	}
}


void plot_events(int ch, eeevector evector, double* medians, double* means)
{
	char name[256];
	sprintf(name, "canv_ch%d", ch);
	
	char title[256];
	sprintf(title, "Eventi ch%d", ch+1); //IT
	
	TCanvas *canv = new TCanvas(name,"canvas ev", 800, 600);
	TH2D *hist = new TH2D("", title, 24, xmin, xmax, ymax/ydraw, -ymax-yacc, ymax+yacc);
	
	
	for (size_t i=0; i<evector.size(); ++i) {              // for every event
		eeevent event1 = evector[i];
		for (size_t j=0; j<event1.hits[ch].size(); ++j) {   // for every hit in ch
			hit hit1 = event1.hits[ch][j];
			hist->Fill(hit1.x, hit1.y);
		}
	}
	
	canv->cd();
	hist->SetStats(kFALSE);
	hist->SetXTitle("Coordinata x (# strip)");
	hist->SetYTitle("Coordinata y (cm)");
	hist->Draw("colz");
	
	
	double prob[1]={0.5};
	printf("Bin\tMean\t\tMedian\n");
	for (int bin=0; bin<24; bin++) {   // compute median and mean
		TH1D *hist1 = hist->ProjectionY(strcat(title, " bin"), bin+1, bin+1);
		if (hist1->GetEntries() > 0) {
			hist1->GetQuantiles(1, medians+bin, prob);
			means[bin] = hist1->GetMean();
		}
		else {
			medians[bin]=0;
			means[bin]=0;
		}
		cout << bin << "\t" << means[bin] << "\t\t" << medians[bin] << endl;
	}
	//hist1->Draw();
	
	cout << endl << "Chamber " << ch << ": " << hist->Integral() << " events" << endl << endl;
}

void plot_offset(double offsets[3][24])
{
	const double strips[24]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}; // yes, I know...
	for (int ch=0; ch<3 ; ch++) {
		char name[256];
		sprintf(name, "canv_off_ch%d", ch);
		
		char title[256];
		sprintf(title, "Offset ch%d", ch);
		
		TCanvas *canv = new TCanvas(name,"canvas off", 800, 600);
		TGraph *g_res = new TGraph(24, strips, offsets[ch]);
		g_res->SetTitle(title);
		g_res->GetXaxis()->SetLimits(xmin, xmax);
		g_res->GetXaxis()->SetTitle("x (# strip)");
		g_res->GetYaxis()->SetTitle("Offset (cm)");
		canv->cd();
		g_res->Draw("AB");
	}
}

void plot_offsets(double medians[3][24], double means[3][24], double resolutions[3][24])
{
	const double strips[24]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}; // yes, I know...
	for (int ch=0; ch<3 ; ch++) {
		char name_med[256];
		sprintf(name_med, "canv_med_ch%d", ch);
		
		char title_med[256];
		sprintf(title_med, "Medians ch%d", ch);
		
		TCanvas *canv_med = new TCanvas(name_med,"canvas med", 800, 600);
		TGraph *g_med = new TGraph(24, strips, medians[ch]);
		g_med->SetTitle(title_med);
		g_med->GetXaxis()->SetLimits(xmin, xmax);  
		canv_med->cd();
		g_med->Draw("AB");
		
		char name_mea[256];
		sprintf(name_mea, "canv_mea_ch%d", ch);
		
		char title_mea[256];
		sprintf(title_mea, "Means ch%d", ch);
		
		TCanvas *canv_mea = new TCanvas(name_mea,"canvas mea", 800, 600);
		TGraph *g_mea = new TGraph(24, strips, means[ch]);
		g_mea->SetTitle(title_mea);
		g_mea->GetXaxis()->SetLimits(xmin, xmax); 
		canv_mea->cd();
		g_mea->Draw("AB");
		
		char name_res[256];
		sprintf(name_res, "canv_res_ch%d", ch);
		
		char title_res[256];
		sprintf(title_res, "Resolution ch%d", ch);
		
		TCanvas *canv_res = new TCanvas(name_res,"canvas res", 800, 600);
		TGraph *g_res = new TGraph(24, strips, resolutions[ch]);
		g_res->SetTitle(title_res);
		g_res->GetXaxis()->SetLimits(xmin, xmax); 
		canv_res->cd();
		g_res->Draw("AB");
	}
}


void rescale_offset(eeevector& evector, double offsets[3][24])
{
	for (size_t evn=0; evn < evector.size(); ++evn) {
		for (int ch=0; ch < 3; ++ch) {
			hitvector& hits = evector[evn].hits[ch];
			for (size_t hitn=0; hitn < hits.size(); ++hitn) {
				int bin = hits[hitn].x;
				hits[hitn].y -= offsets[ch][bin];
			}
		}
	}
}

void cut_y(eeevector& evector)
{
	for (int ch=0; ch < 3; ++ch) {
		char title[256];
		sprintf(title, "Riscalamento lungo y ch%d", ch+1); //IT
		
		TH2D *hist = new TH2D("", title, 24, xmin, xmax, ymax, -ymax-yacc, ymax+yacc);
		
		for (size_t i=0; i<evector.size(); ++i) {              // for every event
			eeevent event1 = evector[i];
			for (size_t j=0; j<event1.hits[ch].size(); ++j) {   // for every hit in ch
				hit hit1 = event1.hits[ch][j];
				hist->Fill(hit1.x, hit1.y);
			}
		}
		
		
		char name[256];
		sprintf(name, "canv_cut_ch%d", ch);
		
		char titley[256];
		sprintf(titley, "Proiezione Y ch%d", ch+1); //IT
		
		TCanvas *canv = new TCanvas(name, "canvas pr", 800, 600);	
		canv->cd();
		
		TH1D *histy = hist->ProjectionY(titley);
		histy->SetXTitle("Coordinata Y (cm)");
		histy->SetYTitle("Eventi");
		histy->Draw();
	}
}

void average_hit_number(eeevector evector)
{
	for (int ch=0; ch<3; ++ch) {
		char name[256];
		sprintf(name, "canv_hits_ch%d", ch);
		
		char title[256];
		sprintf(title, "Numero medio di hit ch%d", ch+1); // IT
		
		TCanvas *canv = new TCanvas(name,"canvas", 800, 600);
		int nbins=20;
		TH1I *hist = new TH1I("", title, nbins, 0, nbins/2);
		
		for (size_t i=0; i != evector.size(); ++i) {  // for every event
			int nhits=evector[i].hits[ch].size();
			hist->Fill(nhits);
		}
		
		canv->cd();
		hist->SetXTitle("Numero di hit");
		hist->SetYTitle("Eventi");
		hist->Draw();
	}
}

void voltage(eeevector evector)
{
	TCanvas *canv = new TCanvas("canv_volt", "canvas", 800, 600);
	int nbins=evector.size();
	double* voltvec0 = new double[nbins];
	double* voltvec1 = new double[nbins];
	double* voltvec2 = new double[nbins];
	double* times    = new double[nbins];
	
	for (size_t i=0; i != evector.size(); ++i) {  // for every event
		voltvec0[i] = evector[i].hv0;
		voltvec1[i] = evector[i].hv1;
		voltvec2[i] = evector[i].hv2;
		times[i]    = evector[i].time();
	}
	TGraph *g0 = new TGraph(nbins, times, voltvec0);
	g0->SetMinimum(8700);
	TGraph *g1 = new TGraph(nbins, times, voltvec1);
	g1->SetLineColor(4);
	TGraph *g2 = new TGraph(nbins, times, voltvec2);
	g2->SetLineColor(2);
	
	canv->cd();
	g0->Draw("AL");
	g1->Draw("same");
	g2->Draw("same");
}

hit compute_hitpoint(int chtest, hit hit1, hit hit2)
{
	double z, za, zb;
	switch (chtest) {
		case 0: z = z0;     za = z1+0.5; zb = z2;     break;
		case 1: z = z1+0.5; za = z2;     zb = z0;     break;
		case 2: z = z2;     za = z0;     zb = z1+0.5; break;
		default: cout << "Unrecognized chamber" << endl; exit(1);
	}
	// line equation is: (x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)=A
	
	double x1 = hit1.x;
	double x2 = hit2.x;
	double y1 = hit1.y;
	double y2 = hit2.y;
	
	double A = (z-za)/(zb-za);
	double x = (x2-x1)*A +x1;
	double y = (y2-y1)*A +y1;
	return {x, y};
} 

void efficiency_calculator(eeevector evector)
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
	
	TCanvas* canv_noi[3];
	TH2D* hist_noi[3];
	
	for (int ch=0; ch < chnum; ++ch) {
		char name[256];
		sprintf(name, "canv_noi_ch%d", ch);
		
		char title[256];
		sprintf(title, "Noise ch%d", ch);
		
		canv_noi[ch] = new TCanvas(name,"canvas noi", 800, 600);
		hist_noi[ch] = new TH2D("", title, 24, xmin, xmax, ymax/ydraw, -ymax-yacc, ymax+yacc);
	}
	
	TCanvas* canv_diffx[3];
	TH1D* hist_diffx[3];
	
	for (int ch=0; ch < chnum; ++ch) {
		char name[256];
		sprintf(name, "canv_diffx_ch%d", ch);
		
		char title[256];
		sprintf(title, "Distribuzione dei punti ricostruiti in x ch%d", ch+1); //IT
		
		canv_diffx[ch] = new TCanvas(name,"canvas diffx", 800, 600);
		hist_diffx[ch] = new TH1D("", title, 2*24, -xmax, xmax);
		hist_diffx[ch]->SetXTitle("Distanza dal punto misurato (# strip)");
		hist_diffx[ch]->SetYTitle("Eventi");
	}
	
	TCanvas* canv_diffy[3];
	TH1D* hist_diffy[3];
	
	for (int ch=0; ch < chnum; ++ch) {
		char name[256];
		sprintf(name, "canv_diffy_ch%d", ch);
		
		char title[256];
		sprintf(title, "Distribuzione dei punti ricostruiti in y ch%d", ch+1); //IT
		
		canv_diffy[ch] = new TCanvas(name,"canvas diffy", 800, 600);
		hist_diffy[ch] = new TH1D("", title, ymax, -ymax-yacc, ymax+yacc);
		hist_diffy[ch]->SetXTitle("Distanza dal punto misurato (cm)");
		hist_diffy[ch]->SetYTitle("Eventi");
	}
	
	for (int ch=0; ch < 3; ++ch) {
		
		char name_exp[256];
		sprintf(name_exp, "canv_exp_ch%d", ch);
		char name_det[256];
		sprintf(name_det, "canv_det_ch%d", ch);
		char name_eff[256];
		sprintf(name_eff, "canv_eff_ch%d", ch);
		
		char title_exp[256];
		sprintf(title_exp, "Expected ch%d", ch);
		char title_det[256];
		sprintf(title_det, "Detected ch%d", ch);
		char title_eff[256];
		sprintf(title_eff, "Efficienza ch%d", ch+1); //IT
		
		//TCanvas *canv_exp = new TCanvas(name_exp,"canvas exp", 800, 600);
		TH2D *hist_exp = new TH2D("", title_exp, 24, xmin, xmax, ymax/ydraw, -ymax-yacc, ymax+yacc);
		//TCanvas *canv_det = new TCanvas(name_det,"canvas det", 800, 600);
		TH2D *hist_det = new TH2D("", title_det, 24, xmin, xmax, ymax/ydraw, -ymax-yacc, ymax+yacc);
		TCanvas *canv_eff = new TCanvas(name_eff,"canvas eff", 800, 600);
		TH2D *hist_eff = new TH2D("", title_eff, 24, xmin, xmax, ymax/ydraw, -ymax-yacc, ymax+yacc);
		
		int seen=0;
		
		for (size_t i=0; i<goodevents[ch].size(); ++i) {              // for every event
			eeevent event1 = goodevents[ch][i];
			hit hit1 = event1.hits[(ch+1)%3][0];
			hit hit2 = event1.hits[(ch+2)%3][0];
			if (hit1.is_inside_scint() && hit2.is_inside_scint()) {
				hit hit_test = compute_hitpoint(ch, hit1, hit2);
				if (hit_test.is_inside_scint()) {
					hist_exp->Fill(hit_test.x, hit_test.y);
					for (size_t h=0; h<event1.hits[ch].size(); ++h) {   // for every hit in chtest
						hit hitf = event1.hits[ch][h];
						hist_diffx[ch]->Fill(hitf.x - hit_test.x);
						hist_diffy[ch]->Fill(hitf.y - hit_test.y);
						if (abs(hitf.x - hit_test.x) <= dx && abs(hitf.y - hit_test.y) <= dy_corr[ch]) {
							hist_det->Fill(hitf.x, hitf.y);
							seen++;
							break;
						}
					}
				}
			}
			else {
				if (!hit1.is_inside_scint()) hist_noi[(ch+1)%3]->Fill(hit1.x, hit1.y);
				if (!hit2.is_inside_scint()) hist_noi[(ch+2)%3]->Fill(hit2.x, hit2.y);
			}
		}
		//~ canv_exp->cd();
		//~ hist_exp->Draw();
		//~ 
		//~ canv_det->cd();
		//~ hist_det->Draw();
		
		for (int binx=0; binx < hist_det->GetNbinsX(); ++binx) {
			for (int biny=0; biny < hist_det->GetNbinsY(); ++biny) {
				if (hist_det->GetBinContent(binx+1,biny+1) > hist_exp->GetBinContent(binx+1,biny+1)) {
					hist_det->SetBinContent(binx+1, biny+1, hist_exp->GetBinContent(binx+1,biny+1)*0.6);
				}
			}
		}
		
		cout << "Events seen in chamber " << ch << ": " << seen << endl;
		
		double eff=0;
		int nhit=0;
		hist_eff->Divide(hist_det, hist_exp);
		
		for (int binx=0; binx < hist_eff->GetNbinsX(); ++binx) {
			for (int biny=0; biny < hist_eff->GetNbinsY(); ++biny) {
				if (hist_eff->GetBinContent(binx+1,biny+1) > 0) {
					eff+=hist_eff->GetBinContent(binx+1,biny+1);
					nhit++;
				}
			}
		}
		
		eff/=nhit;
		cout << "Efficiency of chamber " << ch << ": " << eff << endl;
		
		for (int binx=0; binx < hist_eff->GetNbinsX(); ++binx) {
			for (int biny=0; biny < hist_eff->GetNbinsY(); ++biny) {
				//if (hist_eff->GetBinContent(binx+1, biny+1)==0) hist_eff->SetBinContent(binx+1, biny+1, eff);
				if (ch==0) if (binx==8 || binx==18 || binx==22) hist_eff->SetBinContent(binx+1, biny+1, 0);
				if (ch==1) if (binx==9) hist_eff->SetBinContent(binx+1, biny+1, 0);
				if (ch==2) if (binx==16 || binx==18 || binx==20) hist_eff->SetBinContent(binx+1, biny+1, 0);
			}
		}
		
		canv_eff->cd();
		hist_eff->SetStats(kFALSE);
		hist_eff->SetXTitle("Coordinata x (# strip)");
		hist_eff->SetYTitle("Coordinata y (cm)");
		hist_eff->Draw("colz");
		
		//~ cout << endl << "{";
		//~ for (int binx=0; binx < hist_eff->GetNbinsX(); ++binx) {
			//~ cout << "{";
			//~ for (int biny=0; biny < hist_eff->GetNbinsY(); ++biny) {
				//~ cout << hist_eff->GetBinContent(binx+1,biny+1) << ", ";
			//~ }
			//~ cout << "\b\b}, ";
		//~ }
		//~ cout << "\b\b}" << endl;
	}
	
	for (int ch=0; ch < chnum; ++ch) {
		canv_noi[ch]->cd();
		hist_noi[ch]->SetXTitle("Coordinata x (# strip)");
		hist_noi[ch]->SetYTitle("Coordinata y (cm)");
		hist_noi[ch]->Draw("colz");
	}
	for (int ch=0; ch < chnum; ++ch) {
		canv_diffx[ch]->cd();
		cout << "\n\nx - ch" << ch << ":" << endl;
		hist_diffx[ch]->Fit("gaus");
		hist_diffx[ch]->Draw();
	}
	for (int ch=0; ch < chnum; ++ch) {
		canv_diffy[ch]->cd();
		cout << "\n\ny - ch" << ch << ":" << endl;
		hist_diffy[ch]->Fit("gaus");
		hist_diffy[ch]->Draw();
	}
}

void print_offset(double offsets[chnum][24])
{
	cout << endl << "{";
	for (int ch=0; ch < chnum; ++ch) {
		cout << "{";
		for (int strip=0; strip < 24; ++strip) {
			cout << offsets[ch][strip] << ", ";
		}
		cout << "\b\b}, ";
	}
	cout << "\b\b}" << endl;
}

void thphi_distribution(eeevector evector)
{
	eeevector perfevents;     // events with a single hit/ch
	for (size_t i=0; i != evector.size(); ++i) {  // for every event
		eeevent event1=evector[i];
		if (event1.hits[0].size()==1
		&&  event1.hits[1].size()==1
		&&  event1.hits[2].size()==1) perfevents.push_back(event1);
	}
	cout << "Perfect events: " << perfevents.size() << endl;
	
	TCanvas *canv_th = new TCanvas("canv_th","canvas th", 800, 600);
	TH1D *hist_th = new TH1D("", "Distribuzione in #theta", 101, 0, 90); //IT
	hist_th->SetXTitle("#theta (°)");
	hist_th->SetYTitle("Eventi");
	TCanvas *canv_phi = new TCanvas("canv_phi","canvas eff", 800, 600);
	TH1D *hist_phi = new TH1D("", "Distribuzione in #phi", 78, 0, 360); // IT
	hist_phi->SetXTitle("#phi [deg]");
	hist_phi->SetYTitle("Eventi");
	hist_phi->SetMinimum(0);
	
	int Nev=0;
	
	for (size_t i=0; i != perfevents.size(); ++i) {  // for every event
		eeevent event1=perfevents[i];
		hit hit0 = event1.hits[0][0];
		hit hit1 = event1.hits[1][0];
		hit hit2 = event1.hits[2][0];
		
		hit hit_test = compute_hitpoint(1, hit0, hit2);
		if (abs(hit1.x - hit_test.x) <= dx && abs(hit1.y - hit_test.y) <= dy_corr[1]) {
			double yp = hit0.y-hit2.y;
			double xp = (hit0.x-hit2.x)*3.2;
			hist_phi->Fill(atan2(yp, xp)*180/TMath::Pi()+180);
			
			double yt = sqrt((hit0.x-hit2.x)*(hit0.x-hit2.x)*3.2*3.2 + (hit0.y-hit2.y)*(hit0.y-hit2.y));
			double xt = z0-z2;
			
			hist_th->Fill(atan2(yt, xt)*180/TMath::Pi());
			
			Nev++;
		}
	}
	cout << "Valid events: " << Nev << endl << endl;
	
	canv_phi->cd();
	hist_phi->Draw();
	canv_th->cd();
	hist_th->Draw();
	
}

double fitfun(double* arg, double* par)
{
	double tm=par[0];
	double N=par[1];
	double cexp=par[2];
	double x=*arg;
	
	if (x < tm) return 0;
	return N*exp(-cexp*x);
}

void time_distrib(eeevector evector)
{
	TCanvas *canv_time = new TCanvas("canv_time","canvas time", 800, 600);
	TH1D *hist_time = new TH1D("", "Distribuzione in tempo", 1000, 0, 500);
	
	if (evector.size()<2) return;
	for (size_t ev=1; ev != evector.size(); ++ev) {  // for every event
		eeevent ev1=evector[ev-1];
		if (ev1.hits[0].size()!=0 && ev1.hits[1].size()!=0 && ev1.hits[2].size()!=0) {
			eeevent ev2=evector[ev];
			 if (ev2.hits[0].size()!=0 && ev2.hits[1].size()!=0 && ev2.hits[2].size()!=0) { 
				double dt = ev2.delta_time_ms(ev1);
				//cout << dt << endl;
				hist_time->Fill(dt);
			}
		}
	}
	hist_time->GetXaxis()->SetTitle("Distanza in tempo (ms)");
	hist_time->GetYaxis()->SetTitle("Eventi");
	canv_time->SetLogy();
	canv_time->cd();
	
	TF1 *fexp = new TF1("Fexp", fitfun, 0, 500, 3);
	fexp->SetParameters(10, 1000, 0.02);
	fexp->SetParNames("tm", "Norm", "cexp");
	fexp->SetNpx(1000); // more points
	hist_time->Fit(fexp);
	hist_time->Draw();
}

int parser(string path, string mode="eee")
{
	
	eeevector evector;
	if (mode == "calib") parse_file(path, evector, mode);
	else if (mode == "eee") parse_dir(path.c_str(), evector, mode);
	else {
		cout << "Mode unrecognized" << endl;
		exit(1);
	}
	
	double medians[3][24];
	double means[3][24];
	double reported_medians[chnum][24]={  // medians from all the EEE data of 08-03
		{3.32419, 5.29597, 5.33823, 5.12025, 5.73931, 7.21682, 7.19375, 10.0718, 0, 4.64609, 5.05202, 5.67403, 1.46907, 1.71985, 1.83848, 2.79605, 1.37313, 3.24093, 0, 2.77568, 3.03059, 3.46353, 0, 3.40497},
		{-4.19656, -4.73269, 13.4538, -4.94234, -3.01792, -1.6836, -1.7194, -1.05278, 3.65286, 0, -2.85605, 1.21787, -4.9753, -2.798, -2.71347, -6.35292, 3.20602, 2.11266, -0.457665, 2.8667, 5.13112, 3.02974, 1.93134, 1.49208},
		{-3.75076, 1.01684, -2.0115, 0.215216, -0.798211, -2.71105, -1.02473, -4.51589, 4.16105, 1.97378, 1.69023, -0.589424, -3.37038, -3.05634, -5.50312, -5.30694, 0, -3.74178, 0, -3.56896, 0, -4.82906, -5.9285, -5.58281}
	};
	
	//for (int i = 0; i < 3; i++) plot_events(i, evector, medians[i], means[i]);
	
	//plot_offset(medians, means, resolutions);
	
	rescale_offset(evector, reported_medians);
	//plot_offset(reported_medians);
	
	//for (int i = 0; i < 3; i++) plot_events(i, evector, medians[i], means[i]);
	//plot_offsets(medians, means, resolutions);
	
	//efficiency_calculator(evector);
	
	thphi_distribution(evector);
	
	//time_distrib(evector);
	
	//cut_y(evector);                // needed just one time
	//average_hit_number(evector);   // needed just one time
	//voltage(evector);              // needed just one time
	 
	return 0;
}

int main(int, char** argv)
{
	string path = argv[1];
	string mode = argv[2];
	parser(path, mode);
	
	return 0;
}

/*
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
		case 1: *ch1 = 2; *ch2 = 0; break;
		case 2: *ch1 = 0; *ch2 = 1; break;
		default: cout << "Unrecognized chamber" << endl;
	}
}

*/
/* Code graveyard
*
* 
* 
void build_resolutions(double resolutions[3][24], double medians[3][24], double means[3][24])
{
	for (int ch=0; ch<3; ++ch) {
		for (int strip=0; strip<24; ++strip) {
			resolutions[ch][strip] = abs(medians[ch][strip] - means[ch][strip]);
		}
	}
}
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
