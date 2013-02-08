void plotOverlaidKeysPdf(int mass, int btag = 2) {
  RooWorkspace* w = _file0->Get("w");

  TParameter<int> *n = static_cast<TParameter<int>*>(w->obj("iterations"));

  RooRealVar* mtt = static_cast<RooRealVar*>(w->arg("mtt"));
  mtt->SetTitle("m_{t#bar{t}}");

  TString btagTitle = "";
  if (btag == 1) {
    btagTitle = "one";
  } else {
    btagTitle = "#geq 2";
  }

  TString title = TString::Format("m_{t#bar{t}} - semi-%%s channel, %s b-tagged jets (GeV/c^{2})", btagTitle.Data());

  RooPlot* plot_mu = mtt->frame();
  plot_mu->SetTitle("");
  plot_mu->SetTitleOffset(1.30, "y");
  plot_mu->SetXTitle(TString::Format(title.Data(), "muonic"));
  plot_mu->SetYTitle("Projection of Z' signal PDF");
  RooPlot* plot_e = mtt->frame();
  plot_e->SetTitle("");
  plot_e->SetTitleOffset(1.30, "y");
  plot_e->SetXTitle(TString::Format(title.Data(), "electronic"));
  plot_e->SetYTitle("Projection of Z' signal PDF");

  for (int i = 0; i < n->GetVal(); i++) {
    TString mu = TString::Format("signal_muon_%d_%d", btag, i);
    TString e = TString::Format("signal_electron_%d_%d", btag, i);

    RooAbsPdf* mu_pdf = static_cast<RooAbsPdf*>(w->pdf(mu));
    RooAbsPdf* e_pdf = static_cast<RooAbsPdf*>(w->pdf(e));

    mu_pdf->plotOn(plot_mu, RooFit::LineWidth(1), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::Name("toys"));
    e_pdf->plotOn(plot_e, RooFit::LineWidth(1), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::Name("toys"));

  }

  RooAbsPdf* mu_good_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format("old_signal_muon_%d", btag)));
  RooAbsPdf* e_good_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format("old_signal_electron_%d", btag)));

  mu_good_pdf->plotOn(plot_mu, RooFit::Name("signal"));
  e_good_pdf->plotOn(plot_e, RooFit::Name("signal"));

  plot_mu->Draw();

  TString legendEntry = TString::Format("Z' signal (%d GeV)", mass);
  TString legendEntryToys = TString::Format("Generated Z' signal (%d GeV)", mass);

  double x1 = 0.55, x2 = 0.85;
  double y1 = 0.65, y2 = 0.87;
  if (mass >= 1500) {
    x1 = 0.15;
    x2 = 0.45;
    y1 = 0.15;
    y2 = y1 + 0.22;
  }

  TString output = TString::Format("keyspdf_overlaid_%d_%d_btag_%%s.pdf", mass, btag);

  TLegend *leg1 = new TLegend(x1, y1, x2, y2);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->SetTextFont(42);
  leg1->AddEntry("signal", legendEntry, "L");
  leg1->AddEntry("toys", legendEntryToys, "L");
  leg1->Draw();

  TString f = TString::Format(output.Data(), "muon");
  c1->SaveAs(f);

  plot_e->Draw();
  leg1->Draw();

  f = TString::Format(output.Data(), "electron");
  c1->SaveAs(f);
}
