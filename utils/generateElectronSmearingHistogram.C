{
  	RooRealVar x   (    "x",        "#sigma(p_{T})",     -0.500, 0.500, "GeV/#it{c}");
   	RooRealVar mean(    "mean",     "mean",  -0.005);
	RooRealVar sigma(   "sigma",    "sigma",  0.0199);
	RooRealVar alpha1(  "alpha1",   "alpha1", 0.30);
	RooRealVar alpha2(  "alpha2",   "alpha2", -0.93);
	RooRealVar n1(      "n1",       "n1",     4.0);
	RooRealVar n2(      "n2",       "n2",     10.);
	RooRealVar frac(    "frac2",    "frac",   0.05);

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
