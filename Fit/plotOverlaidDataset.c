{
  RooWorkspace* w = _file0->Get("w");

  TParameter<int> *n = static_cast<TParameter<int>*>(w->obj("iterations"));

  RooRealVar* mtt = static_cast<RooRealVar*>(w->arg("mtt"));
  RooPlot* plot = mtt->frame();

  for (int i = 0; i < n->GetVal(); i++) {
    TString mu = TString::Format("dataset_muon_%d", i);
    TString e = TString::Format("dataset_electron_%d", i);

    RooAbsData* mu_pdf = w->data(mu);
    RooAbsData* e_pdf = w->data(e);

    std::cout << "Dataset events: " << mu_pdf->numEntries() << std::endl;

    mu_pdf->plotOn(plot, RooFit::LineWidth(1), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
  }

  RooAbsData* mu_good_pdf = w->data("binned_dataset_muon");
  mu_good_pdf->plotOn(plot);

  plot->Draw();
}
