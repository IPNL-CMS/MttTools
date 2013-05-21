{
#include "tdrstyle.C"
  setTDRStyle();
  RooWorkspace* w = (RooWorkspace*) _file0->Get("w");

  // Datasets*
  
  RooAbsData* mu_1b = w->data("data_obs_mu_1b");
  RooAbsData* mu_2b = w->data("data_obs_mu_2b");
  RooAbsData* e_1b = w->data("data_obs_e_1b");
  RooAbsData* e_2b = w->data("data_obs_e_2b");

  double totalEntries = mu_1b->sumEntries() + mu_2b->sumEntries() + e_1b->sumEntries() + e_2b->sumEntries();

  RooConstVar mu_1b_frac("mu_1b_frac", "", mu_1b->sumEntries() / totalEntries);
  RooConstVar mu_2b_frac("mu_2b_frac", "", mu_2b->sumEntries() / totalEntries);
  RooConstVar e_1b_frac("e_1b_frac", "", e_1b->sumEntries() / totalEntries);
  RooConstVar e_2b_frac("e_2b_frac", "", e_2b->sumEntries() / totalEntries);

  RooRealVar* mtt = w->var("mtt");
  const float resolution = 50.;
  const int nBinsForHisto = (mtt.getMax() - mtt.getMin() + 0.5) / resolution;
  mtt->setBins(nBinsForHisto);

  RooDataHist merged_dataset("merged_dataset", "merged_dataset", RooArgSet(*mtt), *mu_1b);
  merged_dataset.add(*mu_2b);
  merged_dataset.add(*e_1b);
  merged_dataset.add(*e_2b);

  // PDFs
  RooAbsPdf* global_pdf_mu_1b = w->pdf("global_pdf_mu_1b");
  RooAbsPdf* global_pdf_mu_2b = w->pdf("global_pdf_mu_2b");
  RooAbsPdf* global_pdf_e_1b = w->pdf("global_pdf_e_1b");
  RooAbsPdf* global_pdf_e_2b = w->pdf("global_pdf_e_2b");

  RooAddPdf merged_pdf("merged_pdf", "", RooArgList(*global_pdf_mu_1b, *global_pdf_mu_2b, *global_pdf_e_1b, *global_pdf_e_2b), RooArgList(mu_1b_frac, mu_2b_frac, e_1b_frac, e_2b_frac));

  RooPlot* plot = mtt->frame(nBinsForHisto);
  merged_dataset.plotOn(plot, RooFit::Name("data"));

  // Signal + background
  merged_pdf->plotOn(plot, RooFit::LineColor(kBlue), RooFit::LineWidth(1), RooFit::Name("sig_bkg"));

  // Background only
  merged_pdf->plotOn(plot, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::Components("background_{electron;1-btag},background_{electron;2-btag},background_{muon;1-btag},background_{muon;2-btag}"), RooFit::Name("bkg"));

  // Draw signal only
  double totalSignalEvents = 288.8890 + 275.0030 + 259.0330 + 247.6470;
  RooConstVar mu_1b_signal("mu_1b_signal", "", 288.8890);
  RooConstVar mu_2b_signal("mu_2b_signal", "", 259.0330);
  RooConstVar e_1b_signal("e_1b_signal", "", 275.0030);
  RooConstVar e_2b_signal("e_2b_signal", "", 247.6470);

  // PDFs
  RooAbsPdf* signal_pdf_mu_1b = w->pdf("signal_{muon;1-btag}");
  RooAbsPdf* signal_pdf_mu_2b = w->pdf("signal_{muon;2-btag}");
  RooAbsPdf* signal_pdf_e_1b = w->pdf("signal_{electron;1-btag}");
  RooAbsPdf* signal_pdf_e_2b = w->pdf("signal_{electron;2-btag}");

  RooExtendPdf ext_signal_pdf_mu_1b("ext_signal_pdf_mu_1b", "", *signal_pdf_mu_1b, mu_1b_signal);
  RooExtendPdf ext_signal_pdf_mu_2b("ext_signal_pdf_mu_2b", "", *signal_pdf_mu_2b, mu_2b_signal);
  RooExtendPdf ext_signal_pdf_e_1b("ext_signal_pdf_e_1b", "", *signal_pdf_e_1b, e_1b_signal);
  RooExtendPdf ext_signal_pdf_e_2b("ext_signal_pdf_e_2b", "", *signal_pdf_e_2b, e_2b_signal);

  RooAddPdf merged_ext_signal_pdf("merged_ext_signal_pdf", "", RooArgList(ext_signal_pdf_mu_1b, ext_signal_pdf_mu_2b, ext_signal_pdf_e_1b, ext_signal_pdf_e_2b));

  merged_ext_signal_pdf.plotOn(plot, RooFit::LineColor(39), RooFit::LineStyle(8), RooFit::LineWidth(2), RooFit::Normalization(totalSignalEvents, RooAbsReal::NumEvent), RooFit::Name("sig"));

  float binningSize = (mtt.getBinning().highBound() - mtt.getBinning().lowBound()) / (float) nBinsForHisto;

  plot->SetTitleOffset(1.3, "Y");
  plot->SetTitleOffset(1, "X");

  plot->SetXTitle("#font[62]{M_{t#bar{t}} [GeV/c^{2}]}");

  //plot->SetYTitle("#font[62]{event yield}");
  plot->SetYTitle("#font[62]{Events}");

  // Y title with bin size
  //{
    //plot->SetTitleOffset(1.24, "Y");
    //plot->SetYTitle(TString::Format("#font[62]{Events / (%.0f GeV/c^{2})}", binningSize));
  //}

  plot->SetTitleSize(0.06, "XZ");

  plot->GetXaxis()->SetNdivisions(5, true);
  plot->GetYaxis()->SetRangeUser(13, 1e5);

  plot->Draw();

  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.10);

  gPad->SetGridx(0);
  gPad->SetGridy(0);

  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.05);

  //t.DrawLatex(0.552013, 0.845, "#font[62]{m_{Z'} = 750 GeV, 1.2% width}");
  t.DrawLatex(0.199664, 0.82, "#font[62]{e+#mu, N_{b-tag} #geq 1}");

  TLegend* title = new TLegend(0.15, 0.91, 0.95, 0.96);
  title->SetBorderSize(0);
  title->SetTextAlign(31);
  title->SetTextFont(62);
  title->SetFillColor(10);
  title->SetLineColor(1);
  title->SetMargin(0.12);
  title->SetHeader("CMS, 19.6 fb^{-1}, #sqrt{s} = 8 TeV");
  title->Draw();

  TLegend *leg1 = new TLegend(0.54,0.597902,0.927852,0.818182);
  leg1->SetFillColor(kWhite);
  leg1->SetLineColor(kWhite);
  leg1->SetBorderSize(0);
  leg1->SetTextFont(62);
  leg1->SetMargin(0.20);
  leg1->AddEntry("data","CMS data 2012", "P");
  leg1->AddEntry("sig_bkg", "Signal + background","L");
  leg1->AddEntry("bkg", "Background only", "L");
  leg1->AddEntry("sig", "Z' 750 Gev/c^{2} (1%)", "L");
  leg1->SetTextSize(0.036);
  leg1->Draw();

  gPad->SetLogy();
  gPad->Modified();

  c1->SaveAs("likelihood_projection_with_signal_zprime_750.pdf");
  //c1->SaveAs("likelihood_projection_zprime_750.pdf");
}
