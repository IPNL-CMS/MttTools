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

void plotMttResolutionFromSkim() {

  setTDRStyle();
  tdrStyle->SetLabelOffset(0.02, "XYZ");

  TCanvas *c = new TCanvas("c", "c", 800, 800);

  TTree *mtt_tree = (TTree*) _file0->Get("Mtt");
  Mtt->Draw("mttResolution_AfterChi2/MC_mtt*100>>mtt_resolution(70,-150,150)", "");
  Mtt->Draw("mttResolution_AfterChi2/MC_mtt*100>>mtt_resolution_good(70,-150,150)", "recoJetsAssociatedWithChi2 && eventIsAssociable");
  Mtt->Draw("mttResolution_AfterChi2/MC_mtt*100>>mtt_resolution_wrong(70,-150,150)", "recoJetsAssociatedWithChi2 == 0 && eventIsAssociable");
  Mtt->Draw("mttResolution_AfterChi2/MC_mtt*100>>mtt_resolution_unmatched(70,-150,150)", "eventIsAssociable == 0");

  TH1F *mtt = (TH1F*) gDirectory->Get("mtt_resolution");
  mtt->SetTitle("");
  mtt->GetXaxis()->SetTitle("Resolution (%)");
  mtt->SetLineColor(TColor::GetColor("#542437"));
  mtt->SetLineColor(TColor::GetColor("#8A9B0F"));
  mtt->SetLineColor(TColor::GetColor("#C02942"));
  mtt->SetFillColor(mtt->GetLineColor());
  mtt->SetFillStyle(3004);
  mtt->GetXaxis()->SetTitleOffset(1.2);

  TH1F *mtt_good = (TH1F*) gDirectory->Get("mtt_resolution_good");
  mtt_good->SetTitle("");
  mtt_good->GetXaxis()->SetTitle("Resolution (%)");
  mtt_good->SetLineColor(TColor::GetColor("#542437"));
  mtt_good->SetLineColor(TColor::GetColor("#8A9B0F"));
  mtt_good->SetLineColor(TColor::GetColor("#C02942"));
  mtt_good->SetFillColor(mtt_good->GetLineColor());
  mtt_good->SetFillStyle(1001);
  mtt_good->GetXaxis()->SetTitleOffset(1.2);

  TH1F *mtt_wrong = (TH1F*) gDirectory->Get("mtt_resolution_wrong");
  mtt_wrong->SetTitle("");
  mtt_wrong->GetXaxis()->SetTitle("Resolution (%)");
  mtt_wrong->SetLineColor(TColor::GetColor("#556270"));
  //mtt_wrong->SetLineColor(TColor::GetColor("#8A9B0F"));
  //mtt_wrong->SetLineColor(TColor::GetColor("#C02942"));
  mtt_wrong->SetFillColor(mtt_wrong->GetLineColor());
  mtt_wrong->SetFillStyle(1001);
  mtt_wrong->GetXaxis()->SetTitleOffset(1.2);

  TH1F *mtt_unmatched = (TH1F*) gDirectory->Get("mtt_resolution_unmatched");
  mtt_unmatched->SetTitle("");
  mtt_unmatched->GetXaxis()->SetTitle("Resolution (%)");
  //mtt_unmatched->SetLineColor(TColor::GetColor("#556270"));
  mtt_unmatched->SetLineColor(TColor::GetColor("#ECD078"));
  //mtt_unmatched->SetLineColor(TColor::GetColor("#C02942"));
  mtt_unmatched->SetFillColor(mtt_unmatched->GetLineColor());
  mtt_unmatched->SetFillStyle(1001);
  mtt_unmatched->GetXaxis()->SetTitleOffset(1.2);

  THStack* stack = new THStack();
  stack->Add(mtt_unmatched, "hist");
  stack->Add(mtt_wrong, "hist");
  stack->Add(mtt_good, "hist");

  //TF1* func = new TF1("g", "[0] * TMath::BreitWigner(x, [1], [2])", -mtt->GetRMS(), mtt->GetRMS());
  TF1* func = new TF1("g", "gaus", -mtt->GetRMS() / 2., mtt->GetRMS() / 2);
  g->SetParameter(0, mtt->GetMaximum());
  g->SetParameter(1, 0);
  g->SetParameter(2, 50);
  mtt->Fit(func, "RQN");

  func->SetLineColor(TColor::GetColor("#D95B43"));
  func->SetLineWidth(2);

  stack->Draw("hist");
  func->Draw("same");

  stack->GetXaxis()->SetTitle("Resolution (%)");

  TLatex* latex = new TLatex();
  latex->SetTextFont(42);
  latex->SetTextSize(0.030);
  latex->DrawLatexNDC(0.20, 0.84, TString::Format("mean = %.2f #pm %.2e", func->GetParameter(1), func->GetParError(1)));
  latex->DrawLatexNDC(0.20, 0.8, TString::Format("#sigma = %.2f #pm %.2e", func->GetParameter(2), func->GetParError(2)));
  latex->DrawLatexNDC(0.20, 0.76, TString::Format("RMS = %.2f #pm %.2e", mtt->GetRMS(), mtt->GetRMSError()));

  gPad->Modified();
  gPad->Update();

  TLegend l(0.63, 0.75, 0.9, 0.85);
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.SetTextFont(42);
  l.SetTextSize(0.025);
  l.AddEntry(mtt_good, "Good combinations", "F");
  l.AddEntry(mtt_wrong, "Wrong combinations", "F");
  l.AddEntry(mtt_unmatched, "Non-associable events", "F");
  l.Draw("same");

  c->Print("mtt_resolution.pdf");
}
