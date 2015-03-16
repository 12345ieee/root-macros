#include <iostream>
#include "TGraph.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMath.h"
#include <math.h>

double func(double *x, double *par){

  double x1=x[0];
  double M=par[0]; 
  double cosa=(74-27*M*M*x1);
  double p1= 22*pow(2, 1./3)/(pow((cosa + sqrt(5324+cosa*cosa)), 1./3));
  double p2= pow(2, 2./3)*pow((cosa + sqrt(5324+cosa*cosa)), 1./3);
  //return 1./6*(-2+p1-p2);
  return cosa;
}

int main()
{
	double M=80;
	double* pM=&M;
	double x=0.0001;
	std::cout << func(&x, pM) << std::endl;
	return 18;
}


void MT_studies(){

	const int N=17;
	const int M=12;

	double p_mu[N]={0.018, 0.019, 0.020, 0.021, 0.022,0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031, 0.032, 0.033, 0.034};
	double recoil[N]={0.07841,0.11, 0.2159, 0.2727, 0.2495, 0.2064, 0.1542, 0.0988, 0.05343, 0.02059, -0.00205, -0.01847, -0.0293, -0.03589, -0.03426, -0.02853, -0.0253};

	double p_mu_eta[M]={0.020, 0.021, 0.022,0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031};
	double recoil_eta[M]={0.3627, 0.3119, 0.2862, 0.1708, 0.1040, 0.05081, -0.04301, -0.1049, -0.1676, -0.2293, -0.1631, -0.1926};

	for(int j=0; j<N; j++){
		p_mu[j]=p_mu[j]*p_mu[j];
	}

	for(int i=0; i<M; i++){
		p_mu_eta[i]=p_mu_eta[i]*p_mu_eta[i];
	}

	TGraph *g1= new TGraph(M, p_mu_eta, recoil_eta);
	TGraph *g= new TGraph(N, p_mu, recoil);

	g->GetXaxis()->SetTitle("1/p_{T}(#mu)^2");
	g->GetYaxis()->SetTitle("h*p_{T}(#mu)/p_{T}(#mu)^2");
	g->SetMarkerStyle(20);
	g->SetMarkerColor(2);
	g1->SetMarkerStyle(20);
	//g1->Draw("AP");
	//g->Draw("P");

	TF1 *f1 =new TF1("first_order", "[0]*x+1",0, 1);
	f1->SetParameter(0,-1600);
	//f1->Draw("same");

	TF1 *f2 =new TF1("second_order", "2-[0]*TMath::Sqrt(x)",0, 1);
	f2->SetParameter(0, 80);
	f2->SetLineColor(3);
	//f2->Draw("same");

	TF1 *f3 =new TF1("third_order","(1./6)*(2 - (22 * pow(2, (1./3))/pow(74 - 27*[0]*[0]*x + sqrt(5324 + pow((74 - 27*[0]*[0]*x), 2)), 1./3)) - pow(2, 2./3)*pow((74 - 27*[0]*[0]*x + sqrt(5324 + pow((74 - 27*[0]*[0]*x), 2))), 1./3))", 0.0004, 0.002);
        f3->SetParameter(0, 80);
        f3->SetLineColor(5);
        f3->Draw();


}
