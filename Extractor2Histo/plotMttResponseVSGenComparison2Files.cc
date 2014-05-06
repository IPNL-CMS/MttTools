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

void plotMttResponseVSGenComparison2Files() {

  setTDRStyle();

  TCanvas *c = new TCanvas("c", "c", 800, 800);

  TGraphErrors * mtt_response = (TGraphErrors*) _file0->Get("mtt_gen_vs_mtt_reco_chi2sel/mtt_gen_vs_mtt_reco_chi2sel_response");
  mtt_response->SetTitle("");
  mtt_response->GetYaxis()->SetTitle("m_{t#bar{t}} (GeV)");
  mtt_response->GetXaxis()->SetTitle("m_{t#bar{t}}^{gen} (GeV)");
  //mtt_response->SetLineColor(TColor::GetColor("#542437"));
  //mtt_response->SetLineColor(TColor::GetColor("#8A9B0F"));
  mtt_response->SetLineColor(TColor::GetColor("#C02942"));
  mtt_response->SetMarkerColor(mtt_response->GetLineColor());
  mtt_response->SetMarkerStyle(20);
  mtt_response->SetMarkerSize(1.5);
  //mtt_response->SetFillColor(mtt->GetLineColor());
  //mtt_response->SetFillStyle(3004);
  //mtt_response->GetXaxis()->SetTitleOffset(1.2);

  TGraphErrors * mtt_response_mva = (TGraphErrors*) _file1->Get("mtt_gen_vs_mtt_reco_chi2sel/mtt_gen_vs_mtt_reco_chi2sel_response");
  mtt_response_mva->SetTitle("");
  mtt_response_mva->GetYaxis()->SetTitle("m_{t#bar{t}} (GeV)");
  mtt_response_mva->GetXaxis()->SetTitle("m_{t#bar{t}}^{gen} (GeV)");
  //mtt_response_mva->SetLineColor(TColor::GetColor("#542437"));
  mtt_response_mva->SetLineColor(TColor::GetColor("#8A9B0F"));
  //mtt_response_mva->SetLineColor(TColor::GetColor("#C02942"));
  mtt_response_mva->SetMarkerColor(mtt_response_mva->GetLineColor());
  mtt_response_mva->SetMarkerStyle(20);
  mtt_response_mva->SetMarkerSize(1.5);

  TF1 *line_a = new TF1("line_a", "pol1", 0, 2000);
  mtt_response->Fit(line_a, "QN");
  line_a->SetLineColor(mtt_response->GetLineColor());
  line_a->SetLineWidth(2);

  TF1 *line_b = new TF1("line_b", "pol1", 0, 2000);
  mtt_response_mva->Fit(line_b, "QN");
  line_b->SetLineColor(mtt_response_mva->GetLineColor());
  line_b->SetLineWidth(2);
  
  TMultiGraph* mg = new TMultiGraph();

  mg->Add(mtt_response, "P");
  mg->Add(mtt_response_mva, "P");

  mg->Draw("a");

  mg->GetYaxis()->SetTitle("m_{t#{bar}} (GeV)");
  mg->GetXaxis()->SetTitle("m_{t#bar{t}}^{gen} (GeV)");

  line_a->Draw("same");
  line_b->Draw("same");

  TLatex* latex = new TLatex();
  latex->SetTextFont(42);
  latex->SetTextSize(0.033);
  latex->SetTextColor(line_a->GetLineColor());
  latex->DrawLatexNDC(0.25, 0.84, TString::Format("a = %.2f #pm %.2f", line_a->GetParameter(0), line_a->GetParError(0)));
  latex->DrawLatexNDC(0.25, 0.8, TString::Format("b = %.2f #pm %.2e", line_a->GetParameter(1), line_a->GetParError(1)));

  latex->SetTextColor(line_b->GetLineColor());
  latex->DrawLatexNDC(0.25, 0.74, TString::Format("a = %.2f #pm %.2f", line_b->GetParameter(0), line_b->GetParError(0)));
  latex->DrawLatexNDC(0.25, 0.7, TString::Format("b = %.2f #pm %.2e", line_b->GetParameter(1), line_b->GetParError(1)));

  TLegend l(0.6, 0.15, 0.9, 0.3);
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.SetTextFont(42);
  l.AddEntry(mtt_response, "#chi^{2} algorithm", "P");
  l.AddEntry(mtt_response_mva, "BDT", "P");
  l.Draw("same");

  gPad->Modified();
  gPad->Update();

  c->Print("mtt_response_vs_gen_comparison_chi2_bdt.pdf");
}
