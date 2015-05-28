#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
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
			cout << "Resolution" << i+1 << " = (" << sigma[i] << " \\pm " << errSigma[i] << ") meV" << endl;
		num += sigma[i]/errSigma[i]/errSigma[i];
		den += 1/errSigma[i]/errSigma[i];
	}
	
	
	double resolution = num/den;
	double errResolution = 1/sqrt(den);
	cout << endl << "Resolution = (" << resolution << " \\pm " << errResolution << ") meV" << endl;
}


