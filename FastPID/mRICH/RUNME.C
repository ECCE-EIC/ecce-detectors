#include "mRICH.C"
#include <vector>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>

std::vector<PID*> Detectors;
std::vector<TH1*> Plots;

PID::type   myType = PID::pi_k;  // This means I want to plot the pi_K separation;
double    numSigma = 3.001;        // This means I want to plot teh contour for three sigma separation
double mom=5.;//GeV/c

void RUNME()
{

  gStyle->SetOptStat(0);

  //  Add as many detectors as your heart desires....
  Detectors.push_back( new mRICH(0.00175, 1, 1, mom) ); //  5 psec @ 150 cm
  Detectors.push_back( new mRICH(0.00175, 1, 2, mom) ); // 10 psec @ 100 cm
  Detectors.push_back( new mRICH(0.00175, 1, 3, mom) ); // 20 psec @ 100 cm
  
  //  Select the booking criterion for the performance plot.
  TH1* Performance = new TH1D("Performance","Performance",100,-6,6);  // Format for detector Performance.
  
  //----------------------------------------------------------------------------------
  //  At least in the beginning, you might not want to change the lines below here...
  //----------------------------------------------------------------------------------

  // Now we call upon all detectors to present their descriptions and assumptions.
  for (int i=0; i<Detectors.size(); i++)
    {
      Detectors[i]->description();
    }

  // Now we clone the original performance plot once for every known detector.
  for (int i=0; i<Detectors.size(); i++)
    {
      Plots.push_back( (TH1*)Performance->Clone(Detectors[i]->name().c_str()) );
      Plots[i]->SetLineColor(i+1);  // Works well only for the first 9 colors...
      Plots[i]->SetLineWidth(4);
   }
  
  // Now we fill all the performance plots:
  for (int i=1; i<Plots[0]->GetNbinsX(); i++) // Remember root is weird...bin0 is the underflow bin
    {
      double eta = Plots[0]->GetBinCenter(i);
      for (int j=0; j<Detectors.size(); j++)
	{
	  if (Detectors[j]->valid(eta,mom))
	    {
	      Plots[j]->SetBinContent(i, Detectors[j]->maxP(eta, mom, numSigma, myType) ); 
	      //cout<<numSigma<<"\t"<<Detectors[j]->maxP(eta, numSigma, myType)<<endl;
	    }
	}
    }
  
  TCanvas *c1 = new TCanvas("c1","c1",700,500);
  // Now we display the performance plots
  for (int i=0; i<Plots.size(); i++)
    {
      Plots[i]->Draw("same");
    }

  //  Now we put colored names on top of the plot (too fancy?)...
  TPaveText *pt = new TPaveText(0.15,0.7,0.4,0.90,"NDC");
  for (int i=0; i<Detectors.size(); i++)
    {
      pt->AddText( Detectors[i]->name().c_str() ); 
      ((TText*)pt->GetListOfLines()->Last())->SetTextColor(i+1);
    }
  pt->Draw();

}
