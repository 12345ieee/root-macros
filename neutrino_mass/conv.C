#include <Riostream.h>
#include <cstdio>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TProfile.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"

//definition of the interval
double xlow=18.5504;
const double xupp=18.6176;
const int np=100.0;//number of subintervals
const double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)

double conv( double *x, double *par){

	xlow=18.52;
	// Range of convolution integral

	double step = (18.6-par[2]-xlow) / np;
	double sum=0;
	// Convolution integral of Spectrum and Gaussian by sum
	for(int i=1.0; i<=np; i++) {
		double xx = xlow + (i-.5) * step;
		double fland = par[0]*xx*(par[1]-xx)*TMath::Sqrt((par[1]-xx)*(par[1]-xx)-(par[2]*par[2]));
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
/*
		xx = xupp - (i-.5) * step;
		fland = par[0]*xx*(par[1]-xx)*TMath::Sqrt((par[1]-xx)*(par[1]-xx)-(par[2]*par[2]));
		sum += fland * TMath::Gaus(x[0],xx,par[3]);
*/
	}

	return (step * sum* invsq2pi / par[3]);
}

/*

double conv( double *x, double *par)
{

	return ((par[0]/par[2])*TMath::Exp((x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2])*par[3]*(par[6]-x[0])*(par[4]-(par[6]-x[0]))*TMath::Sqrt((par[4]-(par[6]-x[0]))*(par[4]-(par[6]-x[0]))-par[5]*par[5]))*(b-a)/N;

}

TF1 *conv_fun=new TF1("conv", conv, a, b, 7);

double conv_int( double *x, double *par)
{
	conv_fun->SetParameter(6, x[0]);

	double sum=0;
	for(int i=0;i<N;i++){

		sum=sum+((par[0]/par[2])*TMath::Exp((x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2])*par[3]*(par[6]-x[0])*(par[4]-(par[6]-x[0]))*TMath::Sqrt((par[4]-(par[6]-x[0]))*(par[4]-(par[6]-x[0]))-par[5]*par[5]))*(b-a)/N;
	}
	return sum;

}

 */

void conv() {

	ifstream inputfile_1, inputfile_2;
	inputfile_1.open("gun.dat");
	inputfile_2.open("kurie.dat");
	cout.precision(9);

	TH1D *h=new TH1D("h", "h", 100, xlow, xupp);
	TH1D *h1=new TH1D("h1", "h1", 85, xlow-0.0004, xupp+0.0004);//binning?
	h->Sumw2();
	h1->Sumw2();

	double nentries, energy, sum=0;

	while (!inputfile_1.eof()) {
		inputfile_1 >> energy >> nentries;

		for(int i=0; i<nentries; i++){
			h->Fill(energy);
		}
	}
	while (!inputfile_2.eof()) {
		inputfile_2 >> energy >> nentries;

		for(int i=0; i<nentries; i++){
			h1->Fill(energy);
		}
		sum++;

	}
	cout<<sum<<endl;
	// Start machinery for convolution

	TF1 *convol=new TF1("conv", conv, xlow-0.0004, xupp+0.0004, 4);
	//TF1 *convol=new TF1("conv", conv, xlow-0.0004, 18.600, 4);
	convol->SetParNames("Norm", "Q", "mass", "sigma");
	convol->SetParameters(1e7, 18.60, 0.003, 0.01);
	//convol->SetParLimits(0, 1e6, 1e9);
	//convol->FixParameter(1, 18.60);
	convol->SetParLimits(1, 18.55, 18.65);
	convol->SetParLimits(2, 0, 1e-1);
	//convol->SetParLimits(4, 1e-5, 1);
	//convol->FixParameter(3,1.01353e-02);

	//h1->Fit(convol, "VRL");
	//h1->Draw("");

	TFitResultPtr r = h1->Fit(convol, "VRSL");
		TMatrixDSym cov = r->GetCovarianceMatrix();
		r->Print("V");
		h1->Draw("");
	//convol->Draw("sames");

}
