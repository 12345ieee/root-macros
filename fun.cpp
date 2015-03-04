#include <TCanvas.h>
//#include <TRandom3.h>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
//#include <string>
//#include <cstring>
//#include <fstream>
#include <iostream>

using namespace std;

double fun1(double* arg, double* par)
{
	// arg is for the argument
	// par is for parameters
	double x= *arg;
	
	return x/(1+x*x);
} 

void fun()
{
	// get canvas
	TCanvas *canv = new TCanvas("Canvas", "Canvas Title", 800, 600);
	canv->cd();
	//canv->SetGrid(1,1)
	
	//                   name        fun  m  M
	//TF1 *f1 = new TF1("Funzione", "x", 0, 10);
	
	//                 name        fun*  m  M   npar
	TF1 *f1 = new TF1("Function", fun1, 0, 10, 0);
	//fun->SetParameters(0, 1, 2);
	//fun->SetParNames("c0", "c1", "c2");

	// use this to set properies of the plot
	// base_hist->GetWhatever()->SetWhatever();
	TH1 *base_hist = f1->GetHistogram();
	base_hist->SetTitle("Function plot");
	base_hist->GetXaxis()->SetTitle("Boring X axis (Gbored/s)");
	base_hist->GetYaxis()->SetTitle("Exciting Y axis (Gfun/s)");
	base_hist->GetYaxis()->SetTitleOffset(1.2);
	
	f1->Draw();
	
	// add a legend
	//                         lx    by   rx   ty
	TLegend *leg = new TLegend(0.7, 0.75, 0.99, 0.9);
	leg->SetHeader("Legend Title");
	leg->AddEntry(base_hist, "Function","l");
	leg->AddEntry("a","Random string","l");
	leg->Draw();
	
}
