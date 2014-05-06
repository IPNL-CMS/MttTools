//Breit-Wigner function
Double_t mybw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4);
}

#include "tdrstyle.C"

void plotMttResolution() {

  setTDRStyle();

  TCanvas *c = new TCanvas("c", "c", 800, 800);

  TH1F *mtt = (TH1F*) _file0->Get("mtt_resolution_mva");
  mtt->SetTitle("");
  mtt->SetName("m_{t#bar{t}}");
  mtt->GetXaxis()->SetTitle("m_{t#bar{t}} - m_{t#bar{t}}^{gen} (GeV)");
  mtt->SetLineColor(TColor::GetColor("#542437"));
  mtt->SetLineColor(TColor::GetColor("#8A9B0F"));
  mtt->SetLineColor(TColor::GetColor("#C02942"));
  mtt->SetFillColor(mtt->GetLineColor());
  mtt->SetFillStyle(3004);
  mtt->GetXaxis()->SetTitleOffset(1.2);

  TF1 *voigt = new TF1("voigt", "gaus", -100, 100);
	voigt->SetParName(0, "amp");
	voigt->SetParName(2, "mean");
	voigt->SetParName(1, "sigma");

  voigt->SetParameter(0, mtt->GetMaximum());
	voigt->SetParameter(2, 0);
	voigt->SetParameter(1, 100);

  mtt->Fit(voigt, "RVN");

  voigt->SetLineColor(TColor::GetColor("#8A9B0F"));
  voigt->SetLineWidth(2);

  mtt->Draw();
  voigt->Draw("same");

  TLatex* latex = new TLatex();
  latex->SetTextFont(42);
  latex->SetTextSize(0.033);
  latex->DrawLatexNDC(0.25, 0.84, TString::Format("mean = %.2f", voigt->GetParameter(1)));
  latex->DrawLatexNDC(0.25, 0.8, TString::Format("#sigma = %.2f", voigt->GetParameter(2)));

  gPad->Modified();
  gPad->Update();

  c->Print("mtt_resolution.pdf");
}
