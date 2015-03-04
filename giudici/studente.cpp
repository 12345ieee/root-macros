#include <TCanvas.h>
#include <TRandom3.h>
#include "TH1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"

const int Nstud=10000;
const int Npoints=10;
const int Nmeasures=3;
const double sigma_smearing=1;
const double m=1;
const int q=1;

double f(double x)
{
	return m*x+q;
}

void studente()
{
	// get canvas
	TCanvas *canv = new TCanvas("canvas1", "studente", 800, 600);
	
	// make this canvas the active one
	canv->cd();
	
	// get randgen
	TRandom3 rng(0);
	
	// make data
	// and fill the plots
	
	double data[Nstud][Npoints][Nmeasures];  // here I'll get the students' results
	double xpoints[Npoints]={0,1,2,3,4,5,6,7,8,9};  // x points
	
	int student;
	int point;
	int measure;
	for (student=0; student<Nstud; ++student) {
		for (point=0; point<Npoints; ++point) {
			for (measure=0; measure<Nmeasures; ++measure) {  // for every possible measure
				double y=f(xpoints[point]);                  // get the theoretical result
				double smearing=rng.Gaus(0, sigma_smearing); // we smear it
				data[student][point][measure]=y+smearing;
			}
		}
	}
	
	
	TH1D *scarto = new TH1D("Scarto", "Scarto", 200, -2, 10);
	
	double sigma_mis[Nstud][Npoints];
	
	for (student=0; student<Nstud; ++student) {
		for (point=0; point<Npoints; ++point) {
			double s=0;
			for (measure=0; measure<Nmeasures; ++measure) {  // for every possible measure
				s+= pow(data[student][point][measure] - f(xpoints[point]), 2);
			}
			s/=(Nmeasures-1);
			sigma_mis[student][point]=sqrt(s);
			scarto->Fill(sigma_mis[student][point]);
		}
	}
	
	scarto->Draw();
	
	
	//~ //Create a TGraphErrors
	//~ student=0;
	//~ TGraph *grafico = new TGraph(Npoints*Nmeasures, xpoints, data[student]);
	//~ grafico->SetMarkerStyle(20);
	//~ grafico->SetMinimum(0);
	//~ // Drawing options in http://root.cern.ch/root/html/THistPainter.html
	//~ grafico->Draw("AP");
	//~ 
	//~ // Create Function
	//~ TF1 *f_fit = new TF1("f_fit",straight_line, 0, 10, 2);
	//~ // Set function parameter
	//~ f_fit->SetParameter(0, 1.);//q
	//~ f_fit->SetParameter(1, 0.3);//m
	//~ grafico-> Fit(f_fit, "", "", 0, 10);

	
	
}
