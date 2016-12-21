{
	gROOT->ProcessLine(".x ~/lhcb/lhcbStyle.C"); 
  	RooRealVar x   (    "x",        "#sigma(p_{T})",     -500., 500., "MeV/#it{c}");
   	RooRealVar mean(    "mean",     "mean",   0.0);
	RooRealVar sigma(   "sigma",    "sigma",  39.9);
	RooRealVar alpha1(  "alpha1",   "alpha1", 0.64);
	RooRealVar alpha2(  "alpha2",   "alpha2", -0.93);
	RooRealVar n1(      "n1",       "n1",     2.2);
	RooRealVar n2(      "n2",       "n2",     4.2);
	RooRealVar frac(    "frac2",    "frac",   0.2);

	RooCBShape  cb1("cb1", "cb1",  x, mean,   sigma,  alpha1, n1); 
	RooCBShape  cb2("cb2", "cb2",  x, mean,   sigma,  alpha2, n2); 
	RooAddPdf   pdf("pdf", "pdf",  RooArgList(cb2, cb1), RooArgList( frac ));

	RooDataSet * data = pdf.generate(x, 100000);

	TCanvas * c = new TCanvas();
	RooPlot * frame = x.frame();
	data->plotOn(frame);
	pdf.plotOn(frame);
	frame->Draw();
    c->SaveAs("electronSmearingHistogram.pdf");

    TFile * file = TFile::Open("electronSmearingHistogram.root", "RECREATE");
    RooDataHist * hist_data = new RooDataHist("hist_data", "hist_data", RooArgSet(x), *data);
    TH1D * hist = (TH1D*) hist_data->createHistogram("histE", x);
    hist->Write();
    file->Close();
}
