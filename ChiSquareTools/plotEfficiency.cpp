void histoToGraph(TH1* histo, TGraphErrors* graph) {
  for (int i = 1; i < histo->GetNbinsX(); i++) {
    graph->SetPoint(i - 1, histo->GetBinCenter(i), histo->GetBinContent(i));
    graph->SetPointError(i - 1, 0, histo->GetBinError(i));
  }

  static int markerStyle = 21;

  graph->SetMarkerSize(1);
  graph->SetLineStyle(0);
  graph->SetLineWidth(0);
  graph->SetFillStyle(0);
  graph->SetMarkerStyle(markerStyle++);
  graph->GetXaxis()->SetTitle("Generated m_{tt}");
}

void plotEfficiency() {

  TFile* no_btag = TFile::Open("efficiency_no_htfrac_nobtag.root");
  TFile* with_btag = TFile::Open("efficiency_no_htfrac_with_btag.root");

  TH1* eff_btag = (TH1*) with_btag->Get("chi2_efficiency_vs_mtt");
  TH1* eff_nobtag = (TH1*) no_btag->Get("chi2_efficiency_vs_mtt");
  TH1* eff_jets = (TH1*) no_btag->Get("jets_efficiency_vs_mtt");

  TMultiGraph* mg = new TMultiGraph();

  TGraphErrors* g_eff_btag = new TGraphErrors();
  TGraphErrors* g_eff_nobtag = new TGraphErrors();
  TGraphErrors* g_eff_jets = new TGraphErrors();

  histoToGraph(eff_btag, g_eff_btag);
  histoToGraph(eff_nobtag, g_eff_nobtag);
  histoToGraph(eff_jets, g_eff_jets);

  g_eff_btag->SetMarkerColor(kRed);
  g_eff_btag->SetTitle("#chi^{2} method with b-tags");
  g_eff_nobtag->SetTitle("#chi^{2} method");
  g_eff_jets->SetMarkerColor(kBlue - 3);
  g_eff_jets->SetTitle("taking the four leading jets");

  mg->Add(g_eff_btag/*, "#chi^{2} method with b-tags"*/);
  mg->Add(g_eff_nobtag/*, "#chi^{2} method"*/);
  mg->Add(g_eff_jets/*, "taking four leading jets"*/);

  mg->SetMinimum(0.);
  mg->SetMaximum(1.);

  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("Generated m_{tt}");
  mg->GetYaxis()->SetTitle("Jet selection efficiency");

  float height = 0.88 - 0.67;
  c1->BuildLegend(0.5, 0.15, 0.88, height + 0.15);
}
