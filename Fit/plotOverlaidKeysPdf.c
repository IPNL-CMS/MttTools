{
  RooWorkspace* w = _file0->Get("w");

  TParameter<int> *n = static_cast<TParameter<int>*>(w->obj("iterations"));

  RooRealVar* mtt = static_cast<RooRealVar*>(w->arg("mtt"));
  RooPlot* plot = mtt->frame();

  for (int i = 0; i < n->GetVal(); i++) {
    TString mu = TString::Format("signal_muon_1_%d", i);
    TString e = TString::Format("signal_electron_1_%d", i);

    RooAbsPdf* mu_pdf = static_cast<RooAbsPdf*>(w->pdf(mu));
    RooAbsPdf* e_pdf = static_cast<RooAbsPdf*>(w->pdf(e));

    mu_pdf->plotOn(plot, RooFit::LineWidth(1), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

    //TString _mu = TString::Format("dataset_muon_%d", i);
    //TString _e = TString::Format("dataset_electron_%d", i);

    //RooAbsData* _mu_pdf = w->data(_mu);
    //RooAbsData* _e_pdf = w->data(_e);

//    _mu_pdf->plotOn(plot);
  }

  RooAbsPdf* mu_good_pdf = static_cast<RooAbsPdf*>(w->pdf("old_signal_muon_1"));
  mu_good_pdf->plotOn(plot);

  plot->Draw();
}
