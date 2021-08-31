#include "dualRICH.h"

int particleAngle = 15; // [deg]
bool thresholdMode = false;
double nsigma = 3.;

void
dRICH()
{

  std::vector<PID *> pidDetectors;
  std::vector<int> fillStyle;
  
  /** dualRICH aerogel **/
  auto aerogel = new dualRICH_aerogel();
  aerogel->setThresholdMode(thresholdMode);
  pidDetectors.push_back(aerogel);
  fillStyle.push_back(3004);

  /** dualRICH C2F6 **/
  auto gas = new dualRICH_C2F6();
  gas->setThresholdMode(thresholdMode);
  pidDetectors.push_back(gas);
  fillStyle.push_back(3004);

  /** drawing n-sigma vs. p **/

  auto hNSigmaKpi_aerogel = new TH1F("hNSigmaKpi_aerogel", ";p (GeV/c);", 2000, 0., 200.);
  hNSigmaKpi_aerogel->SetLineColor(kAzure-3);
  hNSigmaKpi_aerogel->SetLineStyle(kSolid);
  hNSigmaKpi_aerogel->SetLineWidth(2);

  auto hNSigmaKpi_gas = new TH1F("hNSigmaKpi_gas", ";p (GeV/c);", 2000, 0., 200.);
  hNSigmaKpi_gas->SetLineColor(kAzure-3);
  hNSigmaKpi_gas->SetLineStyle(kDashed);
  hNSigmaKpi_gas->SetLineWidth(2);

  auto hNSigmaKp_aerogel = new TH1F("hNSigmaKp_aerogel", ";p (GeV/c);", 2000, 0., 200.);
  hNSigmaKp_aerogel->SetLineColor(kRed+1);
  hNSigmaKp_aerogel->SetLineStyle(kSolid);
  hNSigmaKp_aerogel->SetLineWidth(2);

  auto hNSigmaKp_gas = new TH1F("hNSigmaKp_gas", ";p (GeV/c);", 2000, 0., 200.);
  hNSigmaKp_gas->SetLineColor(kRed+1);
  hNSigmaKp_gas->SetLineStyle(kDashed);
  hNSigmaKp_gas->SetLineWidth(2);

  auto particleEta = -log(tan(particleAngle * 0.017453293 * 0.5));
  for (int i = 0; i < hNSigmaKpi_aerogel->GetNbinsX(); ++i) {
    auto p = hNSigmaKpi_aerogel->GetBinCenter(i + 1);
    hNSigmaKpi_aerogel->SetBinContent(i + 1, aerogel->numSigma(particleEta, p, PID::pi_k));
    hNSigmaKpi_gas->SetBinContent(i + 1, gas->numSigma(particleEta, p, PID::pi_k));
    hNSigmaKp_aerogel->SetBinContent(i + 1, aerogel->numSigma(particleEta, p, PID::k_p));
    hNSigmaKp_gas->SetBinContent(i + 1, gas->numSigma(particleEta, p, PID::k_p));
  }

  TLine line;
  line.SetLineColor(kBlack);
  line.SetLineWidth(1);
  line.SetLineStyle(7);

  TLatex latex;
  latex.SetNDC(kTRUE);
  latex.SetTextSize(0.05);
  latex.SetTextAlign(31);

  auto s = new TCanvas("s", "s", 1600, 800);
  s->cd(1)->Divide(2, 1);  

  s->cd(1)->DrawFrame(0., 0.2, 100., 200., ";#it{p} (GeV/#it{c});n#sigma_{K#pi}");
  line.DrawLine(0.0, 3., 100., 3.);
  s->cd(1)->SetLogy();
  hNSigmaKpi_aerogel->Draw("same,l");
  hNSigmaKpi_gas->Draw("same,l");
  latex.DrawLatex(0.88, 0.85, Form("#vartheta_{track} = %d#circ", particleAngle));
  
  s->cd(2)->DrawFrame(0., 0.2, 100., 200., ";#it{p} (GeV/#it{c});n#sigma_{pK}");
  line.DrawLine(0.0, 3., 100., 3.);
  s->cd(2)->SetLogy();
  hNSigmaKp_aerogel->Draw("same,l");
  hNSigmaKp_gas->Draw("same,l");
  latex.DrawLatex(0.88, 0.85, Form("#vartheta_{track} = %d#circ", particleAngle));
  
  s->SaveAs(Form("dRICH.separation.%ddeg.png", particleAngle));
  
  /** drawing 3-sigma separation in eta-p plane **/
  
  std::vector<double> etaf, etab;
  for (double eta = -6.; eta <= 6.; eta += 0.1)
    etaf.push_back(eta);
  etab = etaf;
  std::reverse(std::begin(etab), std::end(etab));

  auto c = new TCanvas("c", "c", 1600, 800);
  c->Divide(2, 1);  
  
  c->cd(1)->DrawFrame(-2.0, 0.1, 5.0, 150.0, ";#eta;#it{p} (GeV/#it{c})");
  c->cd(1)->SetLogy();
  auto kapi_legend = new TLegend(0.2, 0.2, 0.5, 0.4);
  kapi_legend->SetBorderSize(0);
  for (int idet = 0; idet < pidDetectors.size(); ++idet) {

    auto poly = new TPolyLine();
    poly->SetLineColor(kAzure-3);
    poly->SetLineStyle(kSolid);
    poly->SetLineWidth(1);
    poly->SetFillColor(kAzure-3);
    poly->SetFillStyle(fillStyle[idet]);

    for (auto eta : etaf) {
      if (!pidDetectors[idet]->valid(eta, 1.e6)) continue;
      poly->SetNextPoint(eta, pidDetectors[idet]->minP(eta, nsigma, PID::pi_k));
    }
    for (auto eta : etab) {
      if (!pidDetectors[idet]->valid(eta, 1.e6)) continue;
      poly->SetNextPoint(eta, pidDetectors[idet]->maxP(eta, nsigma, PID::pi_k));
    }
    poly->SetNextPoint(poly->GetX()[0], poly->GetY()[0]);
    c->cd(1);
    poly->Draw();
    poly->Draw("f");

    kapi_legend->AddEntry(poly, pidDetectors[idet]->name().c_str(), "f");
  }
  kapi_legend->SetHeader(Form("K/#pi %.0f#sigma separation", nsigma));
  kapi_legend->Draw("same");

  
  c->cd(2)->DrawFrame(-2.0, 0.1, 5.0, 150.0, ";#eta;#it{p} (GeV/#it{c})");
  c->cd(2)->SetLogy();
  auto prka_legend = new TLegend(0.2, 0.2, 0.5, 0.4);
  prka_legend->SetBorderSize(0);
  for (int idet = 0; idet < pidDetectors.size(); ++idet) {
    
    auto poly = new TPolyLine();
    poly->SetLineColor(kRed+1);
    poly->SetLineStyle(kSolid);
    poly->SetLineWidth(1);
    poly->SetFillColor(kRed+1);
    poly->SetFillStyle(fillStyle[idet]);

    for (auto eta : etaf) {
      if (!pidDetectors[idet]->valid(eta, 1.e6)) continue;
      poly->SetNextPoint(eta, pidDetectors[idet]->minP(eta, nsigma, PID::k_p));
    }
    for (auto eta : etab) {
      if (!pidDetectors[idet]->valid(eta, 1.e6)) continue;
      poly->SetNextPoint(eta, pidDetectors[idet]->maxP(eta, nsigma, PID::k_p));
    }
    poly->SetNextPoint(poly->GetX()[0], poly->GetY()[0]);
    c->cd(2);
    poly->Draw();
    poly->Draw("f");

    prka_legend->AddEntry(poly, pidDetectors[idet]->name().c_str(), "f");
  }
  prka_legend->SetHeader(Form("p/K %.0f#sigma separation", nsigma));
  prka_legend->Draw("same");
  
  c->SaveAs("dRICH.png");

}
