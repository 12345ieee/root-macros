#include <TCanvas.h>
#include <TRandom3.h>
#include "TH1.h"
#include "TMath.h"

const long int N=10000000;
const double sigma=1;

void gausstest()
{
	// get canvas
	TCanvas *canv = new TCanvas("nome", "Gaussian test canvas", 800, 600);
	
	// make this canvas the active one
	canv->cd();
	
	// get histogram (empty)
	TH1D *hist = new TH1D("Gaussiana?", "Gaussian test", 1000, -5, 5);
	
	// get randgen
	TRandom3 rng(0);
	
	// make data
	// and fill the histogram
	int i;
	double* gnum= (double*)malloc(N*sizeof(double));
	for (i=0; i<N; ++i) {
		double csi1=rng.Uniform(0,1);
		double csi2=rng.Uniform(0,1);
		gnum[i]=sigma*sqrt(-2*log(csi1))*cos(2*TMath::Pi()*csi2);
		//gnum[i+1]=sigma*sqrt(-log(csi1))*sen(2*TMath::Pi()*csi2);
		hist->Fill(gnum[i]);
	}
	
	hist->Draw();
	
	
}
