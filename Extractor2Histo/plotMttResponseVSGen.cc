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

void plotMttResponse() {

  setTDRStyle();

  TCanvas *c = new TCanvas("c", "c", 800, 800);

  TGraphErrors * mtt_response = (TGraphErrors*) _file0->Get("mtt_gen_vs_mtt_reco_chi2sel/mtt_gen_vs_mtt_reco_chi2sel_response");
  mtt_response->SetTitle("");
  mtt_response->GetYaxis()->SetTitle("m_{t#bar{t}} (GeV)");
  mtt_response->GetXaxis()->SetTitle("m_{t#bar{t}}^{gen} (GeV)");
  mtt_response->SetLineColor(TColor::GetColor("#542437"));
  mtt_response->SetLineColor(TColor::GetColor("#8A9B0F"));
  mtt_response->SetLineColor(TColor::GetColor("#C02942"));
  mtt_response->SetMarkerColor(mtt_response->GetLineColor());
  mtt_response->SetMarkerStyle(20);
  mtt_response->SetMarkerSize(1.5);
  //mtt_response->SetFillColor(mtt->GetLineColor());
  //mtt_response->SetFillStyle(3004);
  //mtt_response->GetXaxis()->SetTitleOffset(1.2);

  TF1 *line = new TF1("line", "pol1", 0, 2000);
  mtt_response->Fit(line, "QN");

  line->SetLineColor(TColor::GetColor("#8A9B0F"));
  line->SetLineWidth(2);

  mtt_response->Draw("AP");
  line->Draw("same");
  mtt_response->Draw("P same");

  TLatex* latex = new TLatex();
  latex->SetTextFont(42);
  latex->SetTextSize(0.033);
  latex->DrawLatexNDC(0.25, 0.84, TString::Format("a = %.2f #pm %.2f", line->GetParameter(0), line->GetParError(0)));
  latex->DrawLatexNDC(0.25, 0.8, TString::Format("b = %.2f #pm %.2e", line->GetParameter(1), line->GetParError(1)));

  gPad->Modified();
  gPad->Update();

  c->Print("mtt_response.pdf");
}
