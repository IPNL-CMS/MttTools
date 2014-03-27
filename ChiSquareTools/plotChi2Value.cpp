void histoToGraph(TH1* histo, TGraphErrors* graph) {
  for (int i = 1; i < histo->GetNbinsX(); i++) {
    graph->SetPoint(i - 1, histo->GetBinCenter(i), histo->GetBinContent(i));
    graph->SetPointError(i - 1, 0, histo->GetBinError(i));
  }

  static int markerStyle = 20;

  graph->SetMarkerSize(1);
  graph->SetLineStyle(0);
  graph->SetLineWidth(0);
  graph->SetFillStyle(0);
  graph->SetMarkerStyle(markerStyle++);
  graph->GetXaxis()->SetTitle("Generated m_{tt}");
}

#include "tdrstyle.C"

void plotChi2Value(TFile* f) {

  setTDRStyle();

  //TFile* no_btag = TFile::Open("selection_efficiency_no_htfrac_nobtag.root");
  //TFile* with_btag = TFile::Open("efficiency_no_htfrac_with_btag.root");

  //TH1* eff_btag = (TH1*) with_btag->Get("chi2_efficiency_vs_mtt");

  TTree* mtt = (TTree*) f->Get("Mtt");

  mtt->Draw("solChi2>>chi2_all(100,0,50)");
  TH1* chi2_all = (TH1*) gDirectory->Get("chi2_all");

  mtt->Draw("bestSolChi2>>chi2_best(100,0,50)");
  TH1* chi2_best = (TH1*) gDirectory->Get("chi2_best");

  TCanvas* c1 = new TCanvas("c", "c", 800, 800);

  chi2_all->Scale(chi2_best->Integral() / chi2_all->Integral());

  chi2_best->SetLineColor(TColor::GetColor("#C44D58"));
  chi2_all->SetLineColor(TColor::GetColor("#556270"));

  chi2_best->SetXTitle("#chi^{2} value");
  chi2_best->SetYTitle("Events");

  chi2_best->Draw();
  chi2_all->Draw("same");

  float height = 0.88 - 0.67;
  TLegend* l = new TLegend(0.5, 0.70, 0.88, height + 0.70);
  l->SetTextFont(42);
  l->SetBorderSize(0);
  l->SetFillColor(kWhite);
  l->AddEntry(chi2_all, "All combinations", "L");
  l->AddEntry(chi2_best, "Good combination", "L");
  l->Draw();

  c1->SetLogy();

  c1->Print("chi2_distribution.pdf");
}
