TLegend* createLegend(TH1* a, TH1* b, bool larger = false) {
  TLegend* l = NULL;
  l = new TLegend(0.40, 0.57, 0.88, 0.88);

  l->SetTextFont(42);
  if (larger)
    l->SetTextSize(0.040);
  l->SetFillColor(kWhite);
  l->AddEntry(a, "", "p");
  l->AddEntry(b, "", "f");

  return l;
}

void plotGoodWrongCombinaison() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(42);
  gStyle->SetTitleFont(42);
  gStyle->SetLegendFont(42);

  TFile* f = TFile::Open("good_wrong_combinaisons.root");
  
  TH1* w_mass = (TH1*) f->Get("w_mass_large");
  TH1* w_mass_wrong = (TH1*) f->Get("w_mass_wrong");

  TH1* h_top_mass = (TH1*) f->Get("hadronic_top_mass_large");
  TH1* h_top_mass_wrong = (TH1*) f->Get("hadronic_top_mass_wrong");

  TH1* l_e_top_mass = (TH1*) f->Get("leptonic_top_mass_e_large");
  TH1* l_e_top_mass_wrong = (TH1*) f->Get("leptonic_top_mass_e_wrong");

  TH1* l_mu_top_mass = (TH1*) f->Get("leptonic_top_mass_mu_large");
  TH1* l_mu_top_mass_wrong = (TH1*) f->Get("leptonic_top_mass_mu_wrong");

  TH1* pt_syst = (TH1*) f->Get("pt_system_large");
  TH1* pt_syst_wrong = (TH1*) f->Get("pt_system_wrong");

  TH1* ht_frac = (TH1*) f->Get("ht_frac");
  TH1* ht_frac_wrong = (TH1*) f->Get("ht_frac_wrong");

   TIter next(gDirectory->GetList());
   TObject *obj;
   while ((obj=next())) {
      if (obj->InheritsFrom("TH1")) {
         TH1 *h = (TH1*)obj;
         std::cout << h << std::endl;
         //h->SetLabelSize(0.03);
         h->SetLabelFont(42, "xy");
         h->SetTitleFont(42, "xy");
      }
   }

  w_mass->SetTitle("Hadronic W mass");
  w_mass->GetXaxis()->SetTitle("Hadronic W mass (GeV)");
  w_mass->SetMarkerStyle(20);
  w_mass->SetMarkerSize(0.8);
  w_mass->Rebin(5);
  w_mass_wrong->SetTitle("Hadronic W mass (all comb.)");

  h_top_mass->SetTitle("Hadronic top mass");
  h_top_mass->GetXaxis()->SetTitle("Hadronic top mass (GeV)");
  h_top_mass->SetMarkerStyle(20);
  h_top_mass->SetMarkerSize(0.8);
  h_top_mass->Rebin(4);
  h_top_mass_wrong->SetTitle("Hadronic top mass (all comb.)");

  l_e_top_mass->SetTitle("Leptonic top mass");
  l_e_top_mass->GetXaxis()->SetTitle("Leptonic top mass (semi-e channel) (GeV)");
  l_e_top_mass->SetMarkerStyle(20);
  l_e_top_mass->SetMarkerSize(0.8);
  l_e_top_mass->Rebin(4);
  l_e_top_mass_wrong->SetTitle("Leptonic top mass (all comb.)");

  l_mu_top_mass->SetTitle("Leptonic top mass");
  l_mu_top_mass->GetXaxis()->SetTitle("Leptonic top mass (semi-mu channel) (GeV)");
  l_mu_top_mass->SetMarkerStyle(20);
  l_mu_top_mass->SetMarkerSize(0.8);
  l_mu_top_mass->Rebin(4);
  l_mu_top_mass_wrong->SetTitle("Leptonic top mass (all comb.)");

  pt_syst->SetTitle("t#bar{t} system p_{t}");
  pt_syst->GetXaxis()->SetTitleSize(0.04);
  pt_syst->GetXaxis()->SetTitle("t#bar{t} system p_{t} (GeV)");
  pt_syst->SetMarkerStyle(20);
  pt_syst->SetMarkerSize(0.8);
  pt_syst->Rebin(4);
  pt_syst_wrong->SetTitle("t#bar{t} system p_{t} (all comb.)");

  ht_frac->SetTitle("H_{t} fraction");
  ht_frac->GetXaxis()->SetTitle("H_{t} fraction");
  ht_frac->SetMarkerStyle(20);
  ht_frac->SetMarkerSize(0.8);
  ht_frac->Rebin(2);

  w_mass_wrong->SetLineColor(8);
  w_mass_wrong->SetFillColor(8);
  w_mass_wrong->SetFillStyle(1001);
  w_mass_wrong->Rebin(4);

  h_top_mass_wrong->SetLineColor(8);
  h_top_mass_wrong->SetFillColor(8);
  h_top_mass_wrong->SetFillStyle(1001);
  h_top_mass_wrong->Rebin(4);

  l_e_top_mass_wrong->SetLineColor(8);
  l_e_top_mass_wrong->SetFillColor(8);
  l_e_top_mass_wrong->SetFillStyle(1001);
  l_e_top_mass_wrong->Rebin(4);

  l_mu_top_mass_wrong->SetLineColor(8);
  l_mu_top_mass_wrong->SetFillColor(8);
  l_mu_top_mass_wrong->SetFillStyle(1001);
  l_mu_top_mass_wrong->Rebin(4);

  pt_syst_wrong->SetLineColor(8);
  pt_syst_wrong->SetFillColor(8);
  pt_syst_wrong->SetFillStyle(1001);
  pt_syst_wrong->Rebin(4);

  ht_frac_wrong->SetLineColor(8);
  ht_frac_wrong->SetFillColor(8);
  ht_frac_wrong->SetFillStyle(1001);
  ht_frac_wrong->Rebin(2);

  TCanvas *c1 = new TCanvas("c1", "c1", 400, 400);

  w_mass->DrawNormalized("hist P");
  w_mass_wrong->DrawNormalized("hist same");
  w_mass->DrawNormalized("P same");
  createLegend(w_mass, w_mass_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_w_mass.pdf");

  h_top_mass->DrawNormalized("hist P");
  h_top_mass_wrong->DrawNormalized("hist same");
  h_top_mass->DrawNormalized("P same");
  createLegend(h_top_mass, h_top_mass_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_hadronic_top_mass.pdf");

  l_e_top_mass->DrawNormalized("hist P");
  l_e_top_mass_wrong->DrawNormalized("hist same");
  l_e_top_mass->DrawNormalized("P same");
  createLegend(l_e_top_mass, l_e_top_mass_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_leptonic_top_mass_e.pdf");

  l_mu_top_mass->DrawNormalized("hist P");
  l_mu_top_mass_wrong->DrawNormalized("hist same");
  l_mu_top_mass->DrawNormalized("P same");
  createLegend(l_mu_top_mass, l_mu_top_mass_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_leptonic_top_mass_mu.pdf");

  pt_syst->DrawNormalized("hist P");
  pt_syst_wrong->DrawNormalized("hist same");
  pt_syst->DrawNormalized("P same");
  createLegend(pt_syst, pt_syst_wrong, true)->Draw();

  c1->Print("good_wrong_combinaisons_pt_syst.pdf");

  ht_frac->DrawNormalized("hist P");
  ht_frac_wrong->DrawNormalized("hist same");
  ht_frac->DrawNormalized("P same");
  createLegend(ht_frac, ht_frac_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_ht_frac.pdf");
}
