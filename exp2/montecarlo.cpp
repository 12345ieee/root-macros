#include <TLegend.h>
#include <TRandom3.h>
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>

const double xmax= 76.1;
const double ymax= 168;

inline bool is_inside_x(double x)
{
	return 0<x && x<xmax;
}

inline bool is_inside_y(double y)
{
	return 0<y && y<ymax;
}

double Traccia(int N_volte){

	TRandom *rand= new TRandom ();
	TCanvas *c1 = new TCanvas("c1","multipads",900,700);
	TCanvas *c2 = new TCanvas("c2","multipads",900,700);
	TF1 *f1= new TF1 ("f1", "sin(x)*pow(cos(x),2)", 0, TMath::Pi()/2);
	TH1D *h1= new TH1D ("Phi", "Distribuzione in Phi", 500, 0, 360);
	TH1D *h2= new TH1D ("Theta", "Distribuzione in Theta",500 , 0, 90);

	int N;
	int Na=0;
	const double z12 = 149   - 95.5;
	const double z23 =  95.5 - 49;

	for(N=0; N<N_volte; N++){
		
		Double_t x1=rand->Uniform(0,xmax);
		Double_t y1=rand->Uniform(0,ymax);
		Double_t theta1=f1->GetRandom();					//Coordinate della prima camera generate
		Double_t phi1=rand->Uniform(0, 2*TMath::Pi());
		
		Double_t d2=z12*tan(theta1);		//Coordinate polari della seconda camera
		Double_t x2=d2*cos(phi1)+x1;
		if (is_inside_x(x2)) {
			
			Double_t y2=d2*sin(phi1)+y1;
			if (is_inside_y(y2)) {
				
				Double_t d3=z23*tan(theta1);		//Coordinate polari della terza camera
				Double_t x3=d3*cos(phi1)+x2;
				if (is_inside_x(x3)) {
					
					Double_t y3=d3*sin(phi1)+y2;
					if (is_inside_y(y3)) {
						
						//~ Double_t x_step=3.2;
						//~ Double_t x_1=floor(x1/x_step);
						//~ Double_t x_2=floor(x2/x_step);
						//~ Double_t x_3=floor(x3/x_step);

						Na++;
						h2->Fill(theta1*180/TMath::Pi());
						h1->Fill(phi1*180/TMath::Pi());

						/*if((x_1==9) || (x_1==19) || (x_1==23) || (x_2==10) || (x_3==17) || (x_3==21))
						
						else {Na++;h2->Fill(theta1*360/6.28);h1->Fill(phi1*360/6.28);}*/
					}
				}
			}
		}
	}
	double acc=(double)Na/N;
	printf("\nN:\t%d\nNa:\t%d\nAcc:\t%f\n", N, Na, acc);

	c1->cd();
	h1->SetXTitle("Phi[deg]");
	h1->SetYTitle("Eventi");
	h1->Draw();


	c2->cd();
	h2->SetXTitle("Theta[deg]");
	h2->SetYTitle("Eventi");
	h2->Draw();

	return acc*1.23648*2*TMath::Pi();

}


		//printf("%f %f %f\n", x1, x2, x3);

		/*
			int pos_x=x/xstep;
			pos_y=;
			double eff_cella=eff[pos_x][pos_y];
			if (eff==1) {x}
			else if (eff==0) {y}
			else {
			double e = rand->Uniform(0,1);
			if (e>eff_cella) {x}
			else {y}
			}
		*/


void montecarlo()
{
	int N_volte=10000000;
	printf("L'accettanza è %f\n", Traccia(N_volte));
}

int main()
{
	montecarlo();
}
