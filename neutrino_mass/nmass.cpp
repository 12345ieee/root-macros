#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFitResult.h"

using namespace std;

void gun(string filename="gun.dat")
{
	const double emin = 18.5;
	const double emax = 18.7;
	
	ifstream gunfile;
	gunfile.open(filename.c_str());
	if (!gunfile) {
		cout << "File does not exist" << endl;
		return;
	}
	
	TCanvas* canv = new TCanvas("cgun");
	TH1D* hgun = new TH1D("hgun", "Calibration plot", 200, emin, emax);
	
	double energy;
	int nevents;
	
	while (!gunfile.eof()) {
		gunfile >> energy >> nevents;
		//cout << energy << " - "<< nevents << endl;
		for (int i=0; i<nevents; ++i)
			hgun->Fill(energy);
	}
	
	TF1* f = new TF1("3gauss", "gaus(0)+gaus(3)+gaus(6)", emin, emax);
	f->SetParameters(
		13300, 18.55, 0.01,
		17700, 18.60, 0.01,
		 8900, 18.65, 0.01);
	
	TFitResultPtr pResult = hgun->Fit(f, "S"); // "S" gets a hold of the result data
	
	canv->cd();
	hgun->Draw();
	
	double sigma[3];
	double errSigma[3];
	double num=0;
	double den=0;
	
	cout << endl;
	for (int i=0; i<3; ++i) {
		sigma[i] = pResult->Value(3*i+2)*1000;
		errSigma[i] = pResult->ParError(3*i+2)*1000;
			cout << "Resolution" << i+1 << " = (" << sigma[i] << " \\pm " << errSigma[i] << ") eV" << endl;
		num += sigma[i]/errSigma[i]/errSigma[i];
		den += 1/errSigma[i]/errSigma[i];
	}
	
	
	double resolution = num/den;
	double errResolution = 1/sqrt(den);
	cout << endl << "Resolution = (" << resolution << " \\pm " << errResolution << ") eV" << endl;
}

double kurie(double x, double A, double Q, double m)
{
	if (x > abs(Q-m)) return 0;
	return A*x*(Q-x)*sqrt((Q-x)*(Q-x)-m*m);
}

double cgaus(double x, double s)
{
	return TMath::Gaus(x, 0, s, kTRUE);
}

const int nsteps  = 100;
const int nsigmas = 4;

const double emin = 18.5504;  // min of file
const double emax = 18.6176;  // max of file
const double estep=  0.0008;  // step for file

const int    nbins= 85;       // should be (emax-emin)/estep + 1

const double pmin = 18.55;   // min of plot: emin-estep/2
const double pmax = 18.618;  // max of plot: emax+estep/2

double kconv(double* arg, double* par)
{
	double x=arg[0];
	double A=par[0];
	double Q=par[1];
	double m=par[2];
	double s=par[3];
	
	double min = pmin - nsigmas*s;
	
	double interval = Q - m - min;
	double step = interval/nsteps;
	
	double hsum=0;
	double y = min + step/2; 
	
	for (int i=0; i<nsteps; ++i, y += step) {
		double Gy = cgaus(x-y, s);
		double Ky = kurie(y, A, Q, m);
		//cout << y << " - " << Gy << " - " << Ky << endl;
		hsum+=Gy*Ky;
	}
	return hsum*step;
}

void nmass(string filename="kurie.dat")
{
	// res from before: (10.009 \pm 0.008) eV
	double sigma = 0.01; // in KeV, keep it simple
	// Q from Giudici's assignment: 18.600 +- 0.005 KeV
	double Q = 18.600;   // in KeV
	
	ifstream file;
	file.open(filename.c_str());
	if (!file) {
		cout << "File does not exist" << endl;
		return;
	}
	
	TCanvas* canv = new TCanvas("ckurie");
	TH1D* hk = new TH1D("hk", "Kurie plot", nbins, pmin, pmax);
	
	double energy;
	int nevents;
	
	while (!file.eof()) {
		file >> energy >> nevents;
		for (int i=0; i<nevents; ++i)
			hk->Fill(energy);
	}
	
	TF1* f = new TF1("kurie", kconv, pmin, pmax, 4);
	f->SetParNames("A", "Q", "m", "sigma");
	//                 A    Q  m     sigma
	f->SetParameters(1.5e7, Q, 0.03, sigma);
	f->FixParameter(1, Q);
	f->FixParameter(3, sigma);
	
	TFitResultPtr pResult = hk->Fit(f, "LS"); // "S" gets a hold of the result data
	
	canv->cd();
	canv->SetLogy();
	hk->Draw("PE");
}
