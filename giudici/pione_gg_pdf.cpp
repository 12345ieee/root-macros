#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
//#include <fstream>
#include <iostream>

using namespace std;


double pdf(double* arg, double* par)
{
	// arg is for the argument
	// par is for parameters
	double x = (*arg)*(*arg);     // electron pair invariant mass
	double me = par[0];  // electron mass
	double mp = par[1];  // pion mass
	
	double frac = 2*me*me/x;
	double p1 = sqrt(1-2*frac);
	double p2 = 1 + frac;
	double p3 = pow(1 - x/(mp*mp), 3);

	return p1*p2*p3/x;
} 

double h(double* arg, double* par)
{
	double x = (*arg);
	double me = par[0];  // electon mass
	double mp = par[1];  // pion mass
	
	double norm = log(mp/(2*me));
	return 2/(norm*x);
}


void pione_gg_pdf()
{
	const double mp=134.977; // MeV
	const double me=0.511; // MeV
	
	// get canvas
	TCanvas *canv = new TCanvas("Canvas", "Canvas Title", 800, 600);
	canv->cd();
	//canv->SetGrid(1,1)
	
	TF1 *f1 = new TF1("pdf", pdf, 2*me, mp, 2);
	f1->SetParameters(me, mp);
	f1->SetParNames("Electron mass", "Pion mass");
	f1->SetNpx(1000);

	TH1 *base_hist = f1->GetHistogram();
	base_hist->SetTitle("PI-> gamma gamma");
	base_hist->GetXaxis()->SetTitle("Boring X axis (Gbored/s)");
	base_hist->GetYaxis()->SetTitle("Exciting Y axis (Gfun/s)");
	base_hist->GetYaxis()->SetTitleOffset(1.4);
	
	
	// Double extraction method
	TRandom3 *rng = new TRandom3(0);
	const int goal=100000;
	const int nbins=10000;
	double ymax=f1->GetMaximum(2*me, mp);
	
	TH1D *hist = new TH1D("", "", nbins, 2*me, mp);
	
	double par[2]={me, mp};
	
	TStopwatch *wtc = new TStopwatch();
	
	//~ int pts=0;
	//~ while (pts<goal) {
		//~ double x = rng->Uniform(2*me, mp);
		//~ double y = rng->Uniform(0, ymax);
		//~ if (y < pdf(&x, par)) {
			//~ pts++;
			//~ hist->Fill(x);
		//~ }
	//~ }
	
	int pts=0;
	while (pts<goal) {
		double csi = rng->Uniform(0, 1);
		double x = 2*me*pow(mp/(2*me), csi);
		double y = rng->Uniform(0, h(&x, par));
		if (y < pdf(&x, par)) {
			pts++;
			hist->Fill(x);
		}
	}
	
	cout << wtc->RealTime() << endl;
	
	// renormalize
	hist->Scale(1./goal);
	hist->Scale(nbins/mp);
	
	hist->Draw();
	f1->Draw("same");
	
	
	TF1 *f2 = new TF1("h", h, 2*me, mp, 2);
	f2->SetParameters(me, mp);
	
	f2->Draw("same");
	
	
}
