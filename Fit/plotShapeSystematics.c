/*
 * syst =
 *   1: JEC
 *   2: JER
 *   3: PU
 */
void doPlot(int mass, int btag, int syst) {

  RooWorkspace* w = _file0->Get("w");

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

  TString systStr;
  if (syst == 1)
    systStr = "jec";
  else if (syst == 2)
    systStr = "jer";
  else if (syst == 3)
    systStr = "pu";
  else if (syst == 4)
    systStr = "pdf";

  std::cout << systStr << std::endl;

  TString mu = TString::Format("signal_mu_%db_%s%%s", btag, systStr.Data());
  TString e = TString::Format("signal_e_%db_%s%%s", btag, systStr.Data());

  RooAbsPdf* mu_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format(mu.Data(), "Up")));
  RooAbsPdf* e_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format(e.Data(), "Up")));

  mu_pdf->plotOn(plot_mu, /*RooFit::LineWidth(1), */RooFit::LineColor(kRed),/*, RooFit::LineStyle(kDashed),*/ RooFit::Name("up"));
  e_pdf->plotOn(plot_e, /*RooFit::LineWidth(1), */RooFit::LineColor(kRed),/*, RooFit::LineStyle(kDashed),*/ RooFit::Name("up"));

  mu_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format(mu.Data(), "Down")));
  e_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format(e.Data(), "Down")));

  mu_pdf->plotOn(plot_mu, /*RooFit::LineWidth(1), */RooFit::LineColor(kGreen),/*, RooFit::LineStyle(kDashed),*/ RooFit::Name("down"));
  e_pdf->plotOn(plot_e, /*RooFit::LineWidth(1), */RooFit::LineColor(kGreen),/*, RooFit::LineStyle(kDashed),*/ RooFit::Name("down"));

  RooAbsPdf* mu_good_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format("signal_mu_%db", btag)));
  RooAbsPdf* e_good_pdf = static_cast<RooAbsPdf*>(w->pdf(TString::Format("signal_e_%db", btag)));

  mu_good_pdf->plotOn(plot_mu, RooFit::Name("nominal"));
  e_good_pdf->plotOn(plot_e, RooFit::Name("nominal"));

  //TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);

  plot_mu->Draw();

  TString systStrPretty;
  if (syst == 1)
    systStrPretty = "JEC";
  else if (syst == 2)
    systStrPretty = "JER";
  else if (syst == 3)
    systStrPretty = "PU";
  else if (syst == 4)
    systStrPretty = "PDF";

  TString legendEntry = TString::Format("Nominal Z' signal (%d GeV)", mass);
  TString legendEntryUp = TString::Format("Z' signal with %s syst. shifted up  (%d GeV)", systStrPretty.Data(), mass);
  TString legendEntryDown = TString::Format("Z' signal with %s syst. shifted down  (%d GeV)", systStrPretty.Data(), mass);

  double x1 = 0.50, x2 = 0.88;
  double y1 = 0.60, y2 = 0.87;
  if (mass >= 1500) {
    x1 = 0.15;
    x2 = 0.45;
    y1 = 0.15;
    y2 = y1 + 0.22;
  }

  TString output = TString::Format("shape_systematic_%s_%d_%d_btag_%%s.pdf", systStr.Data(), mass, btag);

  TLegend *leg1 = new TLegend(x1, y1, x2, y2);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->SetTextFont(42);
  leg1->AddEntry("signal", legendEntry, "L");
  leg1->AddEntry("up", legendEntryUp, "L");
  leg1->AddEntry("down", legendEntryDown, "L");
  leg1->Draw();

  TString f = TString::Format(output.Data(), "muon");
  c1->SaveAs(f);

  plot_e->Draw();
  leg1->Draw();

  f = TString::Format(output.Data(), "electron");
  c1->SaveAs(f);
}

void plotShapeSystematics(int mass) {
  doPlot(mass, 1, 1);
  doPlot(mass, 2, 1);

  doPlot(mass, 1, 2);
  doPlot(mass, 2, 2);

  doPlot(mass, 1, 3);
  doPlot(mass, 2, 3);

  doPlot(mass, 1, 4);
  doPlot(mass, 2, 4);
}
