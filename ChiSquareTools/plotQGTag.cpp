TLegend* createLegend(TH1* a, TH1* b, bool left = false) {
  TLegend* l = NULL;
  if (! left)
    l = new TLegend(0.40, 0.75, 0.88, 0.88);
  else
    l = new TLegend(0.12, 0.75, 0.45, 0.88);

  l->SetTextFont(42);
  l->SetFillColor(kWhite);
  l->SetBorderSize(0);
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

void plotQGTag(TFile* f) {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelFont(42);
  gStyle->SetTitleFont(42);
  gStyle->SetLegendFont(42);

  TH1F* h_qgLikelihood_hadronicBJets = (TH1F*) f->Get("h_qgLikelihood_hadronicBJets");
  TH1F* h_qgLikelihood_leptonicBJets = (TH1F*) f->Get("h_qgLikelihood_leptonicBJets");
  TH1F* h_qgLikelihood_lightJets = (TH1F*) f->Get("h_qgLikelihood_lightJets");
  TH1F* h_qgLikelihood_allJets = (TH1F*) f->Get("h_qgLikelihood_allJets");
  TH1F* h_qgLikelihood_allTTJets = (TH1F*) f->Get("h_qgLikelihood_allTTJets");
  TH1F* h_qgLikelihood_allNoTTJets = (TH1F*) f->Get("h_qgLikelihood_allNoTTJets");

  int binmax = h_qgLikelihood_allJets->GetMaximumBin();
  double qgLikelihood_x = h_qgLikelihood_allJets->GetXaxis()->GetBinCenter(binmax);

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

  int red_color = TColor::GetColor("#C44D58");
  int bkg_color = TColor::GetColor("#ECD078");
  //int bkg_color = TColor::GetColor("#D95B43");

  h_qgLikelihood_allJets->SetTitle("QG Likelihood (all jets)");
  h_qgLikelihood_allJets->SetLineColor(bkg_color);
  h_qgLikelihood_allJets->SetFillColor(bkg_color);
  h_qgLikelihood_allJets->SetFillStyle(1001);
  //h_qgLikelihood_allJets->Rebin(2);

  h_qgLikelihood_allNoTTJets->SetTitle("QG Likelihood (all no TT jets)");
  h_qgLikelihood_allNoTTJets->SetLineColor(bkg_color);
  h_qgLikelihood_allNoTTJets->SetFillColor(bkg_color);
  h_qgLikelihood_allNoTTJets->SetFillStyle(1001);
  //h_qgLikelihood_allNoTTJets->Rebin(2);

  h_qgLikelihood_hadronicBJets->SetTitle("QG Likelihood, hadr. B jets");
  h_qgLikelihood_hadronicBJets->GetXaxis()->SetTitle("QG Likelihood");
  h_qgLikelihood_hadronicBJets->SetMarkerStyle(20);
  h_qgLikelihood_hadronicBJets->SetMarkerSize(0.8);
  h_qgLikelihood_hadronicBJets->SetMarkerColor(red_color);
  //h_qgLikelihood_hadronicBJets->Rebin(2);
  
  h_qgLikelihood_leptonicBJets->SetTitle("QG Likelihood, lept. B jets");
  h_qgLikelihood_leptonicBJets->GetXaxis()->SetTitle("QG Likelihood");
  h_qgLikelihood_leptonicBJets->SetMarkerStyle(20);
  h_qgLikelihood_leptonicBJets->SetMarkerSize(0.8);
  h_qgLikelihood_leptonicBJets->SetMarkerColor(red_color);  
  //h_qgLikelihood_leptonicBJets->Rebin(2);

  h_qgLikelihood_lightJets->SetTitle("QG Likelihood, light jets");
  h_qgLikelihood_lightJets->GetXaxis()->SetTitle("QG Likelihood");
  h_qgLikelihood_lightJets->SetMarkerStyle(20);
  h_qgLikelihood_lightJets->SetMarkerSize(0.8);
  h_qgLikelihood_lightJets->SetMarkerColor(red_color); 
  //h_qgLikelihood_lightJets->Rebin(2);

  h_qgLikelihood_allTTJets->SetTitle("QG Likelihood, all TT jets");
  h_qgLikelihood_allTTJets->GetXaxis()->SetTitle("QG Likelihood");
  h_qgLikelihood_allTTJets->SetMarkerStyle(20);
  h_qgLikelihood_allTTJets->SetMarkerSize(0.8);
  h_qgLikelihood_allTTJets->SetMarkerColor(red_color); 
  //h_qgLikelihood_allTTJets->Rebin(2);


  std::cout << "Report:" << std::endl;
  std::cout <<"_______________" << std::endl << std::endl;
  std::cout << "QG Lieklihood: " << qgLikelihood_x << " +/- " << h_qgLikelihood_allJets->GetRMS() << std::endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);

  h_qgLikelihood_allJets->Scale(h_qgLikelihood_hadronicBJets->Integral() / h_qgLikelihood_allJets->Integral());
  h_qgLikelihood_hadronicBJets->Draw("hist P");
  h_qgLikelihood_allJets->Draw("hist same");
  h_qgLikelihood_hadronicBJets->Draw("P same");
  createLegend(h_qgLikelihood_hadronicBJets, h_qgLikelihood_allJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_hadronicBJets_comparisonAllJets.pdf");

  h_qgLikelihood_allJets->Scale(h_qgLikelihood_leptonicBJets->Integral() / h_qgLikelihood_allJets->Integral());
  h_qgLikelihood_leptonicBJets->Draw("hist P");
  h_qgLikelihood_allJets->Draw("hist same");
  h_qgLikelihood_leptonicBJets->Draw("P same");
  createLegend(h_qgLikelihood_leptonicBJets, h_qgLikelihood_allJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_leptonicBJets_comparisonAllJets.pdf");

  h_qgLikelihood_allJets->Scale(h_qgLikelihood_lightJets->Integral() / h_qgLikelihood_allJets->Integral());
  h_qgLikelihood_lightJets->Draw("hist P");
  h_qgLikelihood_allJets->Draw("hist same");
  h_qgLikelihood_lightJets->Draw("P same");
  createLegend(h_qgLikelihood_lightJets, h_qgLikelihood_allJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_lightJets_comparisonAllJets.pdf");

  h_qgLikelihood_allJets->Scale(h_qgLikelihood_allTTJets->Integral() / h_qgLikelihood_allJets->Integral());
  h_qgLikelihood_allTTJets->Draw("hist P");
  h_qgLikelihood_allJets->Draw("hist same");
  h_qgLikelihood_allTTJets->Draw("P same");
  createLegend(h_qgLikelihood_allTTJets, h_qgLikelihood_allJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_allTTJets_comparisonAllJets.pdf");


  /////////////

  h_qgLikelihood_allNoTTJets->Scale(h_qgLikelihood_hadronicBJets->Integral() / h_qgLikelihood_allNoTTJets->Integral());
  h_qgLikelihood_hadronicBJets->Draw("hist P");
  h_qgLikelihood_allNoTTJets->Draw("hist same");
  h_qgLikelihood_hadronicBJets->Draw("P same");
  createLegend(h_qgLikelihood_hadronicBJets, h_qgLikelihood_allNoTTJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_hadronicBJets_comparisonAllNoTTJets.pdf");

  h_qgLikelihood_allNoTTJets->Scale(h_qgLikelihood_leptonicBJets->Integral() / h_qgLikelihood_allNoTTJets->Integral());
  h_qgLikelihood_leptonicBJets->Draw("hist P");
  h_qgLikelihood_allNoTTJets->Draw("hist same");
  h_qgLikelihood_leptonicBJets->Draw("P same");
  createLegend(h_qgLikelihood_leptonicBJets, h_qgLikelihood_allNoTTJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_leptonicBJets_comparisonAllNoTTJets.pdf");

  h_qgLikelihood_allNoTTJets->Scale(h_qgLikelihood_lightJets->Integral() / h_qgLikelihood_allNoTTJets->Integral());
  h_qgLikelihood_lightJets->Draw("hist P");
  h_qgLikelihood_allNoTTJets->Draw("hist same");
  h_qgLikelihood_lightJets->Draw("P same");
  createLegend(h_qgLikelihood_lightJets, h_qgLikelihood_allNoTTJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_lightJets_comparisonAllNoTTJets.pdf");

  h_qgLikelihood_allNoTTJets->Scale(h_qgLikelihood_allTTJets->Integral() / h_qgLikelihood_allNoTTJets->Integral());
  h_qgLikelihood_allTTJets->Draw("hist P");
  h_qgLikelihood_allNoTTJets->Draw("hist same");
  h_qgLikelihood_allTTJets->Draw("P same");
  createLegend(h_qgLikelihood_allTTJets, h_qgLikelihood_allNoTTJets, true)->Draw();

  c1->Print("plots/QGTag/qgLikelihood_allTTJets_comparisonAllNoTTJets.pdf");
}
