{
  gROOT->SetBatch();
  gSystem->Load("PlotIt_cc");

  PlotIt p("semimu2012.list", 1);

/*  p.plotstack("hNGoodJets_beforesel");
  c1.SaveAs("plots/nJets_beforesel_semimu.pdf");*/

  p.plotstack("hNGoodJets");
  c1.SaveAs("plots/nJets_semimu.pdf");

/*  p.plotstack("hNBtaggedJets_beforesel");
  c1.SaveAs("plots/nBJets_beforesel_semimu.pdf");*/

  p.plotstack("hNBtaggedJets", 1, 1, "Number of b-tagged jets");
  c1.SaveAs("plots/nBJets_semimu.pdf");

  p.plotstack("hmttSelected2b", 5);
  c1.SaveAs("plots/mtt_2btag_semimu.pdf");

/*  p.plotstack("hLeptonPt_beforesel", 4);
  c1.SaveAs("plots/MuonPt_beforesel_semimu.pdf");
*/
  p.plotstack("hLeptonPt", 4);
  c1.SaveAs("plots/MuonPt_semimu.pdf");

  p.plotstack("h1stjetpt", 4);
  c1.SaveAs("plots/firstjet_semimu.pdf");

/*  p.plotstack("h1stjetpt_beforesel", 4);
  c1.SaveAs("plots/firstjet_beforesel_semimu.pdf");*/

  p.plotstack("h2ndjetpt", 4);
  c1.SaveAs("plots/secondjet_semimu.pdf");

/*  p.plotstack("h2ndjetpt_beforesel", 4);
  c1.SaveAs("plots/secondjet_beforesel_semimu.pdf");*/

  p.plotstack("h3rdjetpt", 4);
  c1.SaveAs("plots/thirdjet_semimu.pdf");

/*  p.plotstack("h3rdjetpt_beforesel", 4);
  c1.SaveAs("plots/thirdjet_beforesel_semimu.pdf");*/

  p.plotstack("h4thjetpt", 4, 1, "4^{th} jet p_{T} [GeV/c]");
  c1.SaveAs("plots/fourthjet_semimu.pdf");

/*  p.plotstack("h4thjetpt_beforesel", 4);
  c1.SaveAs("plots/fourthjet_beforesel_semimu.pdf");*/

  p.plotstack("hMET", 4);
  c1.SaveAs("plots/MET_semimu.pdf");

/*  p.plotstack("hMET_beforesel", 4);
  c1.SaveAs("plots/MET_beforesel_semimu.pdf");*/

  p.plotstack("hNVtx", 2);
  c1.SaveAs("plots/nvertex_semimu.pdf");

  p.plotstack_ratio("hNVtx", 2);
  c1.SaveAs("plots/nvertex_ratio_semimu.pdf");

/*  p.plotstack("hNVtx_beforesel", 2);
  c1.SaveAs("plots/nvertex_beforesel_semimu.pdf");*/

  exit(0);
}