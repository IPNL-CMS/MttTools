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

  TH1F *mtt = (TH1F*) _file0->Get("mtt_resolution");
  mtt->SetTitle("");
  mtt->SetName("m_{t#bar{t}}");
  mtt->GetXaxis()->SetTitle("m_{t#bar{t}} - m_{t#bar{t}}^{gen} (GeV)");
  mtt->SetLineColor(TColor::GetColor("#542437"));
  mtt->SetLineColor(TColor::GetColor("#8A9B0F"));
  mtt->SetLineColor(TColor::GetColor("#C02942"));
  mtt->SetFillColor(mtt->GetLineColor());
  mtt->SetFillStyle(3004);
  mtt->GetXaxis()->SetTitleOffset(1.2);

  //TF1* func = new TF1("g", "[0] * TMath::BreitWigner(x, [1], [2])", -100, 100);
  TF1* func = new TF1("g", "gaus", -100, 100);
  g->SetParameter(0, mtt->GetMaximum());
  g->SetParameter(1, 0);
  g->SetParameter(2, 50);
  mtt->Fit(func, "RQN");

  func->SetLineColor(TColor::GetColor("#8A9B0F"));
  func->SetLineWidth(2);

  mtt->Draw("hist");
  func->Draw("same");

  TLatex* latex = new TLatex();
  latex->SetTextFont(42);
  latex->SetTextSize(0.033);
  latex->DrawLatexNDC(0.25, 0.84, TString::Format("mean = %.2f", func->GetParameter(1)));
  latex->DrawLatexNDC(0.25, 0.8, TString::Format("#sigma = %.2f", func->GetParameter(2)));

  gPad->Modified();
  gPad->Update();

  c->Print("mtt_resolution.pdf");
}
