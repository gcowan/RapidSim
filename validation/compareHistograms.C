int compareHistograms(TString mode) {
	TCanvas c1;

	TH1D hpval("pvals","",10,0.,1.);
	hpval.SetMinimum(0.);

	double sumpval(0.), sum(0.);

	TFile *f1 = TFile::Open("${RAPIDSIM_ROOT}/validation/rootfiles/"+mode+"_hists.root");
	TFile *f2 = TFile::Open(mode+"_hists.root");

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
			std::cout << "ERROR in compareHistograms : Histogram " << h1->GetName() << " not found." << std::endl;
			return 1;
		}

		if(TMath::Abs(h1->GetXaxis()->GetXmax()-h2->GetXaxis()->GetXmax())>0.001) {
			std::cout << "ERROR in compareHistograms : Histogram " << h1->GetName() << " has different maxima " << h1->GetXaxis()->GetXmax() << "\t" << h2->GetXaxis()->GetXmax() << std::endl;
			return 1;
		}
		if(TMath::Abs(h1->GetXaxis()->GetXmin()-h2->GetXaxis()->GetXmin())>0.001) {
			std::cout << "ERROR in compareHistograms : Histogram " << h1->GetName() << " has different minima " << h1->GetXaxis()->GetXmin() << "\t" << h2->GetXaxis()->GetXmin() << std::endl;
			return 1;
		}


		double pval = h1->Chi2Test(h2,"UU");
		std::cout << "INFO in compareHistograms : " << h1->GetName() << " :\tp-value = " << pval << std::endl;
		if(pval<0.01) {
			std::cout << "WARNING in compareHistograms : Histograms do not match for " << h1->GetName() << " (p < 1%)" << std::endl;
		}
		hpval.Fill(pval);
		sumpval+=pval;
		++sum;

		h1->Draw();
		h2->SetLineColor(kRed);
		h2->Draw("same");
		TString pName = "plots/";
		pName += mode;
		pName += "_";
		pName += h1->GetName();
		pName += ".pdf";
		c1.SaveAs(pName);

	}
	hpval.Fit("pol0");
	hpval.Draw();
	TString pName = "plots/";
	pName += mode;
	pName += "_pval.pdf";
	c1.SaveAs(pName);

	if(sumpval/sum < 0.1) {
		std::cout << "WARNING in compareHistograms : sum of p-values is unusually low - check plots" << std::endl;
		return 1;
	}
	return 0;
}

