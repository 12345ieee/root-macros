#include <TCanvas.h>
#include <TLegend.h>
//#include <TRandom3.h>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
//#include <fstream>
#include <iostream>
#include "TGraphErrors.h"

using namespace std;

void effplot()
{
	// get canvas
	TCanvas *canv = new TCanvas("Canvas", "Canvas Title", 800, 600);
	canv->cd();
	
	double xv[6]={7800, 8200, 8600, 9000, 9400, 9800};
	double dxv[6]={0};
	double y1v[6]={0.194948, 0.21719, 0.275192, 0.473395, 0.585464, 0.558046};
	double dy1v[6]={0.03, 0.01, 0.008, 0.009, 0.01, 0.02};
	double y2v[6]={0.423499, 0.413033, 0.536609, 0.654653, 0.688038, 0.659541};
	double dy2v[6]={0.06, 0.02, 0.01, 0.01, 0.01, 0.02};
	double y3v[6]={0.532275, 0.479036, 0.546016, 0.653466, 0.677822, 0.661247};
	double dy3v[6]={0.08, 0.03, 0.02, 0.01, 0.01, 0.02};
	
	TGraphErrors* tg1 = new TGraphErrors(6, xv, y1v, dxv, dy1v);
	TGraphErrors* tg2 = new TGraphErrors(6, xv, y2v, dxv, dy2v);
	TGraphErrors* tg3 = new TGraphErrors(6, xv, y3v, dxv, dy3v);
	//f1->SetParameters(0, 1, 2);
	//f1->SetParNames("c0", "c1", "c2");
	//f1->SetNpx(1000); // more points
	tg1->SetMarkerStyle(21);
	tg2->SetMarkerStyle(21);
	tg3->SetMarkerStyle(21);
	tg2->SetMarkerColor(2);
	tg3->SetMarkerColor(4);
	tg2->SetLineColor(2);
	tg3->SetLineColor(4);

	// use this to set properies of the plot
	// base_hist->GetWhatever()->SetWhatever();
	TH1 *base_hist = tg1->GetHistogram();
	base_hist->SetTitle("Efficienza");
	base_hist->GetXaxis()->SetTitle("Tensione (V)");
	base_hist->GetYaxis()->SetTitle("Efficienza");
	base_hist->GetYaxis()->SetTitleOffset(1.2);
	
	tg1->Draw("ALP");
	tg2->Draw("SAMELP");
	tg3->Draw("SAMELP");
	
	// add a legend
	//                           lx    by   rx   ty
	TLegend *leg = new TLegend(0.7, 0.75, 0.99, 0.9);
	//leg->SetHeader("Legenda");
	leg->AddEntry(tg1, "Camera 1","p");
	leg->AddEntry(tg2, "Camera 2","p");
	leg->AddEntry(tg3, "Camera 3","p");
	leg->Draw();
	
}
