#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TH1F.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

void charge() {

  const int N_events = 9999;

    // Constants
    const double tenp9 = pow(10, 9);
    const double resistance = 909.09; // ohms (1k||10k)
    const double scnd_gain = 1.37;    // (1/2)*(1+RH4/RG4); RH4=270, RG4=470

    // Histograms
    TH1F* hist_charge = new TH1F("hist_charge", "hist_charge", 1200, -10., 110.);    
    TH1F* hist_efficiency = new TH1F("hist_efficiency", "Efficiency vs Threshold", 1200, -10., 110.);
 
    // Open input ROOT file
    TFile* infile = TFile::Open("trig_1_2_4_and_charge3.root", "read");
    infile->cd();

    for (int ij = 0; ij < N_events; ij++) {
        // Retrieve raw signal histogram for each event
        char name[300];
        sprintf(name, "raw_%i", ij);
        TH1F* chnl_raw = dynamic_cast<TH1F*>(infile->Get(name));
	cout<<"chnl_raw binwidth "<<chnl_raw->GetMaximum()<<" "<<chnl_raw->GetMaximumBin()<<endl;//0.1 ns
	 double inv_widthx = 1.0 / chnl_raw->GetBinWidth(2);
	 
	 if (!chnl_raw) {
            std::cerr << "Error: Could not retrieve histogram for event " << ij << std::endl;
            continue;
        }
	 //chnl_raw ->Draw();
        // Analysis variables
        double pedestal = 0.0;
	double pedchi = 0.0;
      

	// Fit pedestal
          
	TFitResultPtr ptr = chnl_raw->Fit("pol0", "QSR", "linear_fit", chnl_raw->GetBinCenter(1.), chnl_raw->GetBinCenter(100.*inv_widthx));
	Int_t fitStatus = ptr;
	if (fitStatus == 0) {
	  pedestal = ptr->Parameter(0);
	  pedchi = ptr->Chi2();
	  cout<<"pedestal "<<pedestal<<endl;
	  //pedestal =  chnl_raw->Integral(1, 100*inv_widthx, "width");

	  if(ptr) {ptr=0;}
	  //   cout<<"pedestal "<<pedestal<<endl;
	} else {
	  std::cerr << "Error: Failed pedestal fit for event " << ij << std::endl;
	  delete chnl_raw; chnl_raw=0;
	  if (ptr) {  ptr=0;}
	  continue;
	}

	
	int peak_bin = chnl_raw->GetMaximumBin();

	int tstart = peak_bin -  20.*inv_widthx;
	int tend = peak_bin + 80.*inv_widthx; 
	
	// TLine* line1 = new TLine(chnl_raw->GetBinCenter(tstart),0.0,chnl_raw->GetBinCenter(tstart),500);
	// line1->Draw();
	// line1->SetLineColor(kRed);
	// line1->SetLineWidth(2);

	
	// TLine* line2 = new TLine(chnl_raw->GetBinCenter(tend),0.0,chnl_raw->GetBinCenter(tend),500);
	// line2->Draw();
	// line2->SetLineColor(kRed);
	// line2->SetLineWidth(2);



	double signal = chnl_raw->Integral(tstart, tend, "width");// * 100.0 * inv_widthx;
	//cout<<"signal "<<signal<<endl;
	double charge = (signal - pedestal) / (resistance * scnd_gain);
	//hist_peak->Fill( (signal-pedestal) / (resistance * scnd_gain));
	cout<<charge<<endl;
	hist_charge->Fill(charge );
		
    

     
    // Clean up raw signal histogram
	delete chnl_raw; chnl_raw=0;
}

    double total_integral = hist_charge->Integral();
    for (int bin = 1; bin <= hist_charge->GetNbinsX(); ++bin) {
      double threshold_integral = hist_charge->Integral(bin, hist_charge->GetNbinsX());
      double efficiency = 100. * threshold_integral / total_integral;
      hist_efficiency->SetBinContent(bin, efficiency);
    }

    // Create output ROOT file
    TFile* outfile = new TFile("Charge.root", "RECREATE");
    outfile->cd();

    // Write histograms to output file
    hist_charge->Write(0,TObject::kOverwrite);
    hist_efficiency->Write(0,TObject::kOverwrite);
    // Clean up histograms and close files
    
    outfile->Close();
    delete outfile;
    delete hist_charge;
    infile->Close();
}


//fit_scint_birks.C
Double_t gaussianfun(Double_t *x, Double_t *par) {
  return par[0]*TMath::Gaus(x[0], par[1], par[2]);
}

Double_t langaufun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]*par[1]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  // /*
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]*par[1]; 
  double scale=1; // par[1];  
  double scale2=1; //anormglb; // for notmalisation this is one, otehrwise use the normalisation;
  //  if (scale2 <.1) scale2 = 0.1;
  //  double scale=par[1];
  // Range of convolution integral
  xlow = x[0] - sc * scale*par[3];
  xupp = x[0] + sc * scale*par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(double ij=1.0; ij<=np/2; ij++) {
    xx = xlow + (ij-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]*par[1], kTRUE); // / par[0];
    if (xx>par[1]) { fland *=exp(-(xx-par[1])/par[4]);}
    sum += fland * TMath::Gaus(x[0],xx,scale*par[3]);
    xx = xupp - (ij-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]*par[1], kTRUE); // / par[0];
    if (xx>par[1]) { fland *=exp(-(xx-par[1])/par[4]);}
    sum += fland * TMath::Gaus(x[0],xx,scale*par[3]);
  }
  return (par[2] * step * sum * invsq2pi / (scale2*par[3]));
}

Double_t totalfun(Double_t* x, Double_t* par){
  return langaufun(x,par) + gaussianfun(x, &par[5]);
}









//Plot pulse,  charge spectrum, efficiency 

void plot(){
  gROOT->ForceStyle();
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(0);    
  gStyle->SetOptTitle(1);
  gStyle->SetOptLogy(0);
  
  gStyle->SetStatW(.40); //40);
    gStyle->SetStatH(.20); //30);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetPadLeftMargin(0.19);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.09); //(0.03);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetStatFont(42);        // Times New Roman
    gStyle->SetTextFont(42);        // Times New Roman
    gStyle->SetTitleFont(42,"XYZ"); // Times New Roman
    gStyle->SetLabelFont(42,"XYZ"); // Times New Roman
    gStyle->SetLabelSize(0.06, "XYZ"); // Times New Roman  
    gStyle->SetNdivisions(606, "XYZ");
    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetTitleX(0.55); //40);
    gStyle->SetTitleY(1.0); //40);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.06);
    latex.SetTextFont(42);
    latex.SetTextAlign(1); //(31); // align right

  // Open the ROOT file
  TFile* file = TFile::Open("Charge.root");

  // Check if the file is open successfully
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening ROOT file!" << std::endl;
    return;
  }

  // Retrieve the histograms
  TH1F* hist_charge = (TH1F*)file->Get("hist_charge");
  TH1F* hist_efficiency = (TH1F*)file->Get("hist_efficiency");

  // Check if histograms were retrieved successfully
  if (!hist_charge || !hist_efficiency) {
    std::cerr << "Error retrieving histograms!" << std::endl;
    return;
  }

  // Create a canvas to plot the histograms
  TCanvas* c1 = new TCanvas("c1", "Charge Spectrum and Efficiency", 1200, 500);

  // Divide the canvas into two pads (for side-by-side plotting)
  c1->Divide(2, 1);

  // Plot hist_charge on the left
  c1->cd(1);
  hist_charge->SetTitle("Charge Spectrum;Charge [pC];Counts");
  hist_charge->Rebin(2);
  //hist_charge->Scale(1./hist_charge->Integral());
  //hist_charge->Draw();


  TF1* fitfunc = new TF1("fitfunc",totalfun,-10.,100.0, 8);//8
  fitfunc->SetParameter(0, 0.2); //histcharge[ijk][ij]->GetRMS());
  fitfunc->SetParameter(1, 4.5); //histcharge[ijk][ij]->GetMean());
  fitfunc->SetParameter(2, 5500); //histcharge[ijk][ij]->Integral());
  fitfunc->SetParameter(3, 1.4);
  fitfunc->SetParameter(4, 1.0);

  int zerobin = hist_charge->FindBin(0);
  double binwidth =  hist_charge->GetBinWidth(zerobin);     
  double pedheight = hist_charge->GetBinContent(zerobin);


  fitfunc->SetParameter(5, pedheight);    
  //if(ijk==3 || ijk==5)fitfunc->SetParLimits(5, 0.*pedheight, 5000);
  fitfunc->SetParLimits(5, 0.*pedheight, 2.0*pedheight);
  fitfunc->SetParameter(6, 0.0);
  fitfunc->SetParLimits(6, -0.25, 0.25);
  fitfunc->SetParameter(7, 0.2);
  fitfunc->SetParLimits(7, -0.5, 0.5);

  hist_charge->Fit(fitfunc, "LBMR");
  hist_charge->GetXaxis()->SetRangeUser(0.,60.0);
  
  latex.DrawLatex(0.48,0.8,Form("MPV: %.2f #pm %.2f pC",fitfunc->GetParameter(1), fitfunc->GetParError(1)));

  // Plot hist_efficiency on the right
  c1->cd(2);
  gStyle->SetOptLogy(1);
  hist_efficiency->SetTitle("Efficiency vs Threshold;Threshold [pC];Efficiency (%)");
  hist_efficiency->Draw();
  hist_efficiency->GetXaxis()->SetRangeUser(0.,3.0);
 



}
