void make_plots()
{
    gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C");
    TFile * f = TFile::Open("Bs2Jpsiphi_hists.root");
    TH1F * Bs_M = (TH1F*)f->Get("Bs_M");
    TH1F * Bs_PT = (TH1F*)f->Get("Bs_PT");
    TH1F * Bs_M_Kp2pip = (TH1F*)f->Get("Bs_M_Kp2pip");
    TH1F * Bs_M_Kp2pip_Km2pim = (TH1F*)f->Get("Bs_M_Kp2pip_Km2pim");
    TH1F * Jpsi_M = (TH1F*)f->Get("Jpsi_M");
    TH1F * phi_M = (TH1F*)f->Get("phi_M");

    TCanvas * c0 = new TCanvas();
    Bs_M->Draw();
    Bs_M->GetYaxis()->SetTitle("Events / ( 0.002 GeV/#it{c}^{2} )");
    Bs_M->GetXaxis()->SetTitle("#it{m}(#it{B}_{#it{s}}^{0}) [GeV/#it{c}^{2}]");
    c0->SaveAs("Bs2JpsiPhi_Bs_M.pdf");

    TCanvas * c1 = new TCanvas();
    Bs_PT->Draw();
    Bs_PT->GetYaxis()->SetTitle("Events / ( 0.002 GeV/#it{c}^{2} )");
    Bs_PT->GetXaxis()->SetTitle("#it{p_{#rm T}}(#it{B}_{#it{s}}^{0}) [GeV/#it{c}^{2}]");
    c1->SaveAs("Bs2JpsiPhi_Bs_PT.pdf");

    TCanvas * c2 = new TCanvas();
    Jpsi_M->Draw();
    Jpsi_M->GetYaxis()->SetTitle("Events / ( 0.002 GeV/#it{c}^{2} )");
    Jpsi_M->GetXaxis()->SetTitle("#it{m}(#it{J}/#psi) [GeV/#it{c}^{2}]");
    c2->SaveAs("Bs2JpsiPhi_Jpsi_M.pdf");

    TCanvas * c3 = new TCanvas();
    phi_M->Draw();
    c3->SetLogy();
    phi_M->GetYaxis()->SetTitle("Events / ( 0.002 GeV/#it{c}^{2} )");
    phi_M->GetXaxis()->SetTitle("#it{m}(#phi) [GeV/#it{c}^{2}]");
    c3->SaveAs("Bs2JpsiPhi_phi_M.pdf");

    TCanvas * c4 = new TCanvas();
    Bs_M_Kp2pip->Draw();
    Bs_M_Kp2pip->GetYaxis()->SetTitle("Events / ( 0.002 GeV/#it{c}^{2} )");
    Bs_M_Kp2pip->GetXaxis()->SetTitle("#it{m}(#it{B}_{#it{s}}^{0} (K^{#plus}#rightarrow #pi^{#plus})) [GeV/#it{c}^{2}]");
    c4->SaveAs("Bs2JpsiPhi_Bs_M_Kp2pip.pdf");

    TCanvas * c5 = new TCanvas();
    Bs_M_Kp2pip_Km2pim->Draw();
    Bs_M_Kp2pip_Km2pim->GetYaxis()->SetTitle("Events / ( 0.002 GeV/#it{c}^{2} )");
    Bs_M_Kp2pip_Km2pim->GetXaxis()->SetTitle("#it{m}(#it{B}_{#it{s}}^{0} (K^{#plus}#rightarrow #pi^{#plus},K^{#minus}#rightarrow #pi^{#minus})) [GeV/#it{c}^{2}]");
    c5->SaveAs("Bs2JpsiPhi_Bs_M_Kp2pip_Km2pim.pdf");

    f->ls();
}
