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

void plotS0PositiveNegativeGen() {

  setTDRStyle();

  TCanvas *c = new TCanvas("c", "c", 800, 800);

  TH1 * positive = (TH1*) _file0->Get("mtt_gen_nosel_positive");
  TH1 * negative = (TH1*) _file0->Get("mtt_gen_nosel_negative");
  TH1 * sum = (TH1*) _file0->Get("mtt_gen_nosel");
  positive->SetTitle("");
  positive->GetYaxis()->SetTitle("");
  positive->GetXaxis()->SetTitle("m_{t#bar{t}} (GeV)");
  positive->SetLineColor(TColor::GetColor("#542437"));
  positive->SetLineColor(TColor::GetColor("#8A9B0F"));
  positive->SetLineColor(TColor::GetColor("#C02942"));
  positive->SetMarkerColor(positive->GetLineColor());
  positive->SetMarkerStyle(20);
  positive->SetMarkerSize(1.5);
  positive->GetXaxis()->SetTitleOffset(1.3);

  negative->SetTitle("");
  negative->SetLineColor(TColor::GetColor("#542437"));
  negative->SetLineColor(TColor::GetColor("#8A9B0F"));
  //negative->SetLineColor(TColor::GetColor("#C02942"));
  negative->SetMarkerColor(negative->GetLineColor());
  negative->SetMarkerStyle(20);
  negative->SetMarkerSize(1.5);
  //mtt_resolution->SetFillColor(mtt->GetLineColor());
  //mtt_resolution->SetFillStyle(3004);
  //mtt_resolution->GetXaxis()->SetTitleOffset(1.2);

  //TF1 *line = new TF1("line", "pol1", 0, 2000);
  //mtt_resolution->Fit(line, "QN");

  //line->SetLineColor(TColor::GetColor("#8A9B0F"));
  //line->SetLineWidth(2);
  
  sum->SetLineColor(TColor::GetColor("#53777A"));
  sum->SetLineStyle(5);
  
  float max = positive->GetMaximum();
  float min = negative->GetMinimum();

  positive->Draw("hist");

  positive->SetMaximum(max * 1.2);
  positive->SetMinimum(min * 1.2);

  negative->Draw("hist same");

  sum->Draw("hist same");

  //TMultiGraph* mg = new TMultiGraph();
  //mg->Add(mtt_resolution, "P");
  //mg->Add(mtt_resolution_four_jets, "P");

  //mg->Draw("a");
  //mg->GetYaxis()->SetTitle("Resolution (%)");
  //mg->GetXaxis()->SetTitle("m_{t#bar{t}}^{gen} (GeV)");
  //line->Draw("same");
  //mtt_resolution->Draw("P same");

  //TLatex* latex = new TLatex();
  //latex->SetTextFont(42);
  //latex->SetTextSize(0.033);
  //latex->DrawLatexNDC(0.25, 0.84, TString::Format("a = %.2f #pm %.2f", line->GetParameter(0), line->GetParError(0)));
  //latex->DrawLatexNDC(0.25, 0.8, TString::Format("b = %.2f #pm %.2e", line->GetParameter(1), line->GetParError(1)));
  
  TLegend l(0.6, 0.15, 0.9, 0.3);
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.SetTextFont(42);
  l.SetTextSize(0.025);
  l.AddEntry(positive, "Positive weight events", "L");
  l.AddEntry(negative, "Negative weight events", "L");
  l.AddEntry(sum, "All events", "L");
  l.Draw("same");

  gPad->Modified();
  gPad->Update();

  c->Print("s0_positive_negative.pdf");
}
