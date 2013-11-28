TLegend* createLegend(TH1* a, TH1* b, bool left = false, bool larger = false) {
  TLegend* l = NULL;
  if (! left)
    l = new TLegend(0.40, 0.57, 0.88, 0.88);
  else
    l = new TLegend(0.12, 0.57, 0.45, 0.88);

  l->SetTextFont(42);
  if (larger)
    l->SetTextSize(0.040);
  l->SetFillColor(kWhite);
  l->AddEntry(a, "", "p");
  l->AddEntry(b, "", "f");

  return l;
}

//Breit-Wigner function
Double_t bw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4);
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

  int binmax = pt_syst->GetMaximumBin();
  double ptsyst_x = pt_syst->GetXaxis()->GetBinCenter(binmax);

  binmax = ht_frac->GetMaximumBin();
  double htfrac_x = ht_frac->GetXaxis()->GetBinCenter(binmax);

   TIter next(gDirectory->GetList());
   TObject *obj;
   while ((obj=next())) {
      if (obj->InheritsFrom("TH1")) {
         TH1 *h = (TH1*)obj;
         //h->SetLabelSize(0.03);
         h->SetLabelFont(42, "xy");
         h->SetTitleFont(42, "xy");
      }
   }

  w_mass->SetTitle("Hadronic W mass");
  w_mass->GetXaxis()->SetTitle("Hadronic W mass (GeV)");
  w_mass->SetMarkerStyle(20);
  w_mass->SetMarkerSize(0.8);
  w_mass->Rebin(2);
  w_mass_wrong->SetTitle("Hadronic W mass (all comb.)");

  h_top_mass->SetTitle("Hadronic top mass");
  h_top_mass->GetXaxis()->SetTitle("Hadronic top mass (GeV)");
  h_top_mass->SetMarkerStyle(20);
  h_top_mass->SetMarkerSize(0.8);
  h_top_mass->Rebin(2);
  h_top_mass_wrong->SetTitle("Hadronic top mass (all comb.)");

  l_e_top_mass->SetTitle("Leptonic top mass");
  l_e_top_mass->GetXaxis()->SetTitle("Leptonic top mass (semi-e channel) (GeV)");
  l_e_top_mass->SetMarkerStyle(20);
  l_e_top_mass->SetMarkerSize(0.8);
  l_e_top_mass->Rebin(2);
  l_e_top_mass_wrong->SetTitle("Leptonic top mass (all comb.)");

  l_mu_top_mass->SetTitle("Leptonic top mass");
  l_mu_top_mass->GetXaxis()->SetTitle("Leptonic top mass (semi-mu channel) (GeV)");
  l_mu_top_mass->SetMarkerStyle(20);
  l_mu_top_mass->SetMarkerSize(0.8);
  l_mu_top_mass->Rebin(2);
  l_mu_top_mass_wrong->SetTitle("Leptonic top mass (all comb.)");

  pt_syst->SetTitle("t#bar{t} system p_{t}");
  pt_syst->GetXaxis()->SetTitleSize(0.04);
  pt_syst->GetXaxis()->SetTitle("t#bar{t} system p_{t} (GeV)");
  pt_syst->SetMarkerStyle(20);
  pt_syst->SetMarkerSize(0.8);
  pt_syst->Rebin(2);
  pt_syst_wrong->SetTitle("t#bar{t} system p_{t} (all comb.)");

  ht_frac->SetTitle("H_{t} fraction");
  ht_frac->GetXaxis()->SetTitle("H_{t} fraction");
  ht_frac->SetMarkerStyle(20);
  ht_frac->SetMarkerSize(0.8);
  //ht_frac->Rebin(2);

  w_mass_wrong->SetLineColor(8);
  w_mass_wrong->SetFillColor(8);
  w_mass_wrong->SetFillStyle(1001);
  w_mass_wrong->Rebin(2);

  h_top_mass_wrong->SetLineColor(8);
  h_top_mass_wrong->SetFillColor(8);
  h_top_mass_wrong->SetFillStyle(1001);
  h_top_mass_wrong->Rebin(2);

  l_e_top_mass_wrong->SetLineColor(8);
  l_e_top_mass_wrong->SetFillColor(8);
  l_e_top_mass_wrong->SetFillStyle(1001);
  l_e_top_mass_wrong->Rebin(2);

  l_mu_top_mass_wrong->SetLineColor(8);
  l_mu_top_mass_wrong->SetFillColor(8);
  l_mu_top_mass_wrong->SetFillStyle(1001);
  l_mu_top_mass_wrong->Rebin(2);

  pt_syst_wrong->SetLineColor(8);
  pt_syst_wrong->SetFillColor(8);
  pt_syst_wrong->SetFillStyle(1001);
  pt_syst_wrong->Rebin(2);

  ht_frac_wrong->SetLineColor(8);
  ht_frac_wrong->SetFillColor(8);
  ht_frac_wrong->SetFillStyle(1001);
  //ht_frac_wrong->Rebin(2);

  // Fit plots
  TF1* top_gaussian = new TF1("f1", "gaus", 150, 210);
  //TF1* top_gaussian = new TF1("f1", "gaus", 0, 2000);
  //TF1* top_gaussian = new TF1("f1", bw, 0, 2000, 3);
  top_gaussian->SetParameter(0, 500);
  top_gaussian->SetParameter(1, 173);
  top_gaussian->SetParameter(2, 40);
  top_gaussian->SetLineColor(kRed);
  top_gaussian->SetLineWidth(2);

  TF1* top_gaussian_lept_e = new TF1("f3", "gaus", 150, 190);
  //TF1* top_gaussian_lept_e = new TF1("f3", bw, 120, 240, 3);
  top_gaussian_lept_e->SetParameter(0, 500);
  top_gaussian_lept_e->SetParameter(1, 173);
  top_gaussian_lept_e->SetParameter(2, 40);
  top_gaussian_lept_e->SetLineColor(kRed);
  top_gaussian_lept_e->SetLineWidth(2);

  TF1* top_gaussian_lept_mu = new TF1("f4", "gaus", 150, 190);
  top_gaussian_lept_mu->SetParameter(1, 173);
  top_gaussian_lept_mu->SetLineColor(kRed);
  top_gaussian_lept_mu->SetLineWidth(2);

  TF1* w_gaussian = new TF1("f2", "gaus", 70, 100);
  //TF1* w_gaussian = new TF1("f2", "gaus", 0, 2000);
  //TF1* w_gaussian = new TF1("f2", bw, 70, 100, 3);
  w_gaussian->SetParameter(0, 500);
  w_gaussian->SetParameter(1, 80);
  w_gaussian->SetParameter(2, 10);
  w_gaussian->SetLineColor(kRed);
  w_gaussian->SetLineWidth(2);

  w_mass->Fit(w_gaussian, "R");
  h_top_mass->Fit(top_gaussian, "QR");
  l_mu_top_mass->Fit(top_gaussian_lept_mu, "RQ");
  l_e_top_mass->Fit(top_gaussian_lept_e, "RQ");

  std::cout << "Report:" << std::endl;
  std::cout <<"_______________" << std::endl << std::endl;
  //if (! CSV_MODE) {
    std::cout << "hadronic W mass: " << w_gaussian->GetParameter(1) << " +/- " << w_gaussian->GetParameter(2) << std::endl;
    std::cout << "hadronic top mass: " << top_gaussian->GetParameter(1) << " +/- " << top_gaussian->GetParameter(2) << std::endl;
    std::cout << "leptonic top mass (semi-mu): " << top_gaussian_lept_mu->GetParameter(1) << " +/- " << top_gaussian_lept_mu->GetParameter(2) << std::endl;
    std::cout << "leptonic top mass (semi-e): " << top_gaussian_lept_e->GetParameter(1) << " +/- " << top_gaussian_lept_e->GetParameter(2) << std::endl;

    std::cout << "TT system pt: " << ptsyst_x << " +/- " << pt_syst->GetRMS() << std::endl;
    std::cout << "HT frac: " << htfrac_x << " +/- " << ht_frac->GetRMS() << std::endl;
  //} else {
    //std::cout << w_gaussian->GetParameter(1) << "\t" << w_gaussian->GetParameter(2) << std::endl;
    //std::cout << top_gaussian->GetParameter(1) << "\t" << top_gaussian->GetParameter(2) << std::endl;
    //std::cout << top_gaussian_lept_mu->GetParameter(1) << "\t" << top_gaussian_lept_mu->GetParameter(2) << std::endl;
    //std::cout << top_gaussian_lept_e->GetParameter(1) << "\t" << top_gaussian_lept_e->GetParameter(2) << std::endl;
  //}

  TCanvas *c1 = new TCanvas("c1", "c1", 400, 400);

  w_mass_wrong->Scale(w_mass->Integral() / w_mass_wrong->Integral());
  w_mass->Draw("hist P");
  w_mass_wrong->Draw("hist same");
  w_mass->Draw("P same");
  w_gaussian->Draw("same");
  createLegend(w_mass, w_mass_wrong)->Draw("same");

  c1->Print("good_wrong_combinaisons_w_mass.pdf");

  h_top_mass_wrong->Scale(h_top_mass->Integral() / h_top_mass_wrong->Integral());
  h_top_mass->Draw("hist P");
  h_top_mass_wrong->Draw("hist same");
  h_top_mass->Draw("P same");
  top_gaussian->Draw("same");
  createLegend(h_top_mass, h_top_mass_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_hadronic_top_mass.pdf");

  l_e_top_mass_wrong->Scale(l_e_top_mass->Integral() / l_e_top_mass_wrong->Integral());
  l_e_top_mass->Draw("hist P");
  l_e_top_mass_wrong->Draw("hist same");
  l_e_top_mass->Draw("P same");
  top_gaussian_lept_e->Draw("same");
  createLegend(l_e_top_mass, l_e_top_mass_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_leptonic_top_mass_e.pdf");

  l_mu_top_mass_wrong->Scale(l_mu_top_mass->Integral() / l_mu_top_mass_wrong->Integral());
  l_mu_top_mass->Draw("hist P");
  l_mu_top_mass_wrong->Draw("hist same");
  l_mu_top_mass->Draw("P same");
  top_gaussian_lept_mu->Draw("same");
  createLegend(l_mu_top_mass, l_mu_top_mass_wrong)->Draw();

  c1->Print("good_wrong_combinaisons_leptonic_top_mass_mu.pdf");

  pt_syst->DrawNormalized("hist P");
  pt_syst_wrong->DrawNormalized("hist same");
  pt_syst->DrawNormalized("P same");
  createLegend(pt_syst, pt_syst_wrong, false, true)->Draw();

  c1->Print("good_wrong_combinaisons_pt_syst.pdf");

  ht_frac->DrawNormalized("hist P");
  ht_frac_wrong->DrawNormalized("hist same");
  ht_frac->DrawNormalized("P same");
  createLegend(ht_frac, ht_frac_wrong, true)->Draw();

  c1->Print("good_wrong_combinaisons_ht_frac.pdf");
}
