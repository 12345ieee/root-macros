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
	double x = *arg;     // electron pair invariant mass
	double me = par[0];  // electon mass
	double mp = par[1];  // pion mass
	
	double frac = 2*me*me/x;
	double p1 = sqrt(1-2*frac);
	double p2 = 1 + frac;
	double p3 = pow(1 - x/(mp*mp), 3);
	double I = log(mp/me)-3.5;

	return (1./I)*p1*p2*p3/x;
} 

void pione_gg_pdf()
{
	const double mp=134.977; // MeV
	const double me=0.511; // MeV
	
	// get canvas
	TCanvas *canv = new TCanvas("Canvas", "Canvas Title", 800, 600);
	canv->cd();
	//canv->SetGrid(1,1)
	
	TF1 *f1 = new TF1("pdf", pdf, 4*me*me, mp*mp, 2);
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
	const int goal=10000;
	const int nbins=1000;
	double ymax=f1->GetMaximum(4*me*me, mp*mp);
	
	TH1D *hist = new TH1D("", "", nbins, 4*me*me, mp*mp);
	hist->GetXaxis()->SetRange(4*me*me,mp*mp);
	
	double par[2]={me, mp};
	
	int pts=0;
	while (pts<goal) {
		double x = rng->Uniform(4*me*me, mp*mp);
		double y = rng->Uniform(0, ymax);
		if (y < pdf(&x, par)) {
			pts++;
			hist->Fill(x);
		}
	}
	// renormalize
	hist->Scale(1./goal);
	hist->Scale(nbins/(mp*mp));
	
	hist->Draw();
	f1->Draw("same");
}
