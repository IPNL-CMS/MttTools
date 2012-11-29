{
  RooWorkspace* w = _file0->Get("w");

  TParameter<int> *n = static_cast<TParameter<int>*>(w->obj("iterations"));

  RooRealVar* mtt = static_cast<RooRealVar*>(w->arg("mtt"));
  RooPlot* plot = mtt->frame();

  RooAbsData* dataset = (RooAbsData*) w->data("binned_dataset_muon");
  dataset->plotOn(plot);

  plot->Draw();
}
