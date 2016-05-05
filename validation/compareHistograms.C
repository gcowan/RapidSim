int compareHistograms(TString mode) {
	int status(0);

	TFile *f1 = TFile::Open("../validation/rootfiles/"+mode+"_hists.root");
	TFile *f2 = TFile::Open("../validation/"+mode+"_hists.root");

	if(!f1 || !f2) {
		std::cout << "ERROR in compareHistograms : Files do not exist for mode " << mode << std::endl;
		return 1;
	}

	TIter next(f1->GetListOfKeys());
	TKey *key;
	
	while ((key = (TKey*)next())) {
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom("TH1")) continue;
		TH1 *h1 = (TH1*)key->ReadObj();
		TH1* h2 = (TH1*)f2->Get(h1->GetName());

		if(!h2) {
			std::cout << "WARNING in compareHistograms : Histogram " << h1->GetName() << " not found." << std::endl;
			status=1;
			continue;
		}

		double chi2 = h1->Chi2Test(h2,"UU");
		std::cout << "INFO in compareHstograms : " << h1->GetName() << " :\tchi2 = " << chi2 << std::endl;
		if(chi2<1.0) {
			std::cout << "WARNING in compareHistograms : Histograms do not match for " << h1->GetName() << std::endl;
			status=1;
		}
	}

	return status;
}

