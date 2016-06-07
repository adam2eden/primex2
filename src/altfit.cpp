#include "primex.h"
#include "alt_model_fit.h"

void plotresult(vector<TH1F*>& h, UImanager& uimanager);

int main(int argc, char* argv[]) //input in this program only refers to the input initial parameters, and output means input data etc. and result
{
	UImanager uimanager(argc, argv, "altfit");

	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1);

	TString inrootdir("./fitroot/");

	//load input root file with rotated mass histograms and sidebands
	TFile *inroot = new TFile(inrootdir + "rotatedmass" + uimanager.output_filename("data") + ".root");
	if (!inroot->IsOpen()) exit(open_err(string(inrootdir + "rotatedmass" + uimanager.output_filename("data") + ".root")));
    cout << "use " << inrootdir + "rotatedmass" + uimanager.output_filename("data") + ".root" << endl;

	vector<TH1F*> h(uimanager.output_nbins(), NULL), hside(uimanager.output_nbins(), NULL), omega(uimanager.output_nbins(), NULL), mcpeak(uimanager.output_nbins(), NULL);

	TVectorD *ratios = (TVectorD*)inroot->Get("ratios");
	for (int i = 0; i<uimanager.output_nbins(); i++){
		TString hname(Form("hrotd_%d", i));
		h[i] = (TH1F*)inroot->Get(hname);
		h[i]->SetDirectory(0);
        
		hside[i] = (TH1F*)inroot->Get(Form("hside_%d", i));
		hside[i]->SetDirectory(0);
		hside[i]->Scale((*ratios)[0]);
	}
	inroot->Close();

	//load omega bkg from M.C.
	inroot = new TFile(inrootdir + "rotatedmass" + uimanager.output_filename("omega") + ".root");
	if (!inroot->IsOpen()) exit(open_err(string(inrootdir + "rotatedmass" + uimanager.output_filename("omega") + ".root")));

	for (int i = 0; i < uimanager.output_nbins(); i++){
		omega[i] = (TH1F*)inroot->Get(Form("homg_%d", i));
		omega[i]->SetDirectory(0);
	}
	inroot->Close();

	//load peak fitting model from M.C. if using M.C. to fit rotated mass distribution
	if (uimanager.get_method() == 1) {
		inroot = new TFile(inrootdir + "fitmodel_mc" + uimanager.output_filename("fit") + ".root"); //corrected M.C. peak model
		if (!inroot->IsOpen()) exit(open_err(inrootdir + "fitmodel_mc" + uimanager.output_filename("fit") + ".root"));
        cout << "use " << inrootdir + "fitmodel_mc" + uimanager.output_filename("fit") + ".root" << endl;
		for (int i = 0; i < uimanager.output_nbins(); i++){
			mcpeak[i] = (TH1F*)inroot->Get(Form("haltinvm_%d", i));
			mcpeak[i]->SetDirectory(0);
		}
		inroot->Close();
	}

	vector<vector<double> > init_par;
	fitrod fitting(uimanager, init_par);

	//create output root files
	int NDF, Npfits, index_angle, npar = fitting.npar();
	double chi2, Npi0, Npi0_err, angle;
	vector<double> parameters(npar, 0.), errors(npar, 0.);

	string infilename;
	if (uimanager.get_method() == 1) infilename = "alt_model_correct_parameters" + uimanager.input_filename("fit") + ".dat";
	else infilename = "offsets" + uimanager.input_filename("fit") + ".dat";
	ifstream inputfile(infilename);

	if (inputfile.is_open()) {
		for (int i = 0; i < uimanager.input_nbins(); i++) {
			vector<double> oneline(npar, 0.);
            double err, chi2;
			if (uimanager.input_nbins() == uimanager.output_nbins()) {
				for (int j = 0; j < npar; j++) inputfile >> oneline[j] >> err;
                inputfile >> chi2;
				init_par.push_back(oneline);
			}
			else {
				int start = int(uimanager.get_input_angles()[i] / (max_angle / nangle) + 0.5);
				int end = int(uimanager.get_input_angles()[i + 1] / (max_angle / nangle) + 0.5);
				for (int j = start; j < end; j++) {
					for (int k = 0; k < npar; k++) inputfile >> oneline[k] >> err;
                    inputfile >> chi2;
					init_par.push_back(oneline);
				}
			}
		}
        fitting.setinitpar(init_par);
        cout << "use parameters from " << infilename << endl;
	}
	else {
		cout << "can't find input parameters; fit without them" << endl;
	}

	TFile *outroot = new TFile(TString("altfit" + uimanager.output_filename("fit") + ".root"), "RECREATE");
	TTree *fitresult = new TTree("fitresult", "fitresult");
	fitresult->Branch("index_angle", &index_angle, "index_angle/I");
	fitresult->Branch("npar", &npar, "npar/I");
	fitresult->Branch("angles", &angle, "angles/D");
	fitresult->Branch("parameters", &parameters[0], "parameters[npar]/D");
	fitresult->Branch("Npi0", &Npi0, "Npi0/D");
	fitresult->Branch("Npi0_err", &Npi0_err, "Npi0_err/D");
	fitresult->Branch("errors", &errors[0], "errors[npar]/D");
	fitresult->Branch("chi2", &chi2, "chi2/D");
	fitresult->Branch("NDF", &NDF, "NDF/I");
	fitresult->Branch("Npfits", &Npfits, "Npfits/I");

	ofstream output("offsets" + uimanager.output_filename("fit") + ".dat");

	for (int i = 0; i < uimanager.output_nbins(); i++) {
		angle = (uimanager.get_output_angles()[i] + uimanager.get_output_angles()[i + 1]) / 2;
		index_angle = i + 1;

		fitting.sethist(h[i], hside[i], mcpeak[i], omega[i]);
        if (uimanager.get_method() == 1) parameters[1] = -100.;
		fitting.initialize(parameters, angle);
		TFitResultPtr r = fitting.fitting(parameters, angle, "Q");
		chi2 = r->Chi2();
		NDF = r->Ndf();
		Npfits = r->NFreeParameters();
		fitting.calc_Npi0(parameters, Npi0, Npi0_err);

		for (int j = 0; j < npar - 2; j++) errors[j] = r->Error(j);
		for (int j = 0; j < npar - 2; j++) cout << r->Parameter(j) << " "; cout << endl;

		for (int j = 0; j < npar - 1; j++) output << parameters[j] << " " << errors[j] << " ";
		output << parameters[npar - 1] << " " << errors[npar - 1] << " " << chi2 / NDF << endl;
		
		cout << "-------------------------------------------------" << endl;
		cout << "****************** " << index_angle << ": angle: " << angle << " *************************" << endl;
		cout << "-------------------------------------------------" << endl;

		fitresult->Fill();
	}

	output.close();

	outroot->cd();
	fitresult->Write();
	for (int i = 0; i < uimanager.output_nbins(); ++i) h[i]->Write();
	outroot->Close();

    if (uimanager.get_method() == 1) fitting.make_plot();
	plotresult(h, uimanager);

    return 0;
}

void plotresult(vector<TH1F*>& h, UImanager& uimanager) {
	TString pdfname("altfit" + uimanager.output_filename("fit") + ".pdf");

	TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
	c1->Divide(2, 2);
	TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
	c2->Divide(2, 2);
	c1->SaveAs(pdfname + "[");

	for (int i = 0; i < uimanager.output_nbins(); i++){
		if (i / 4 >= uimanager.output_nbins() / 4)c2->cd(i % 4 + 1);
		else c1->cd(i % 4 + 1);
		h[i]->SetTitle(Form("%d: ", i + 1) + (TString)h[i]->GetTitle());
		h[i]->GetXaxis()->SetRangeUser(uimanager.fit_low_limit() - 0.02, uimanager.fit_high_limit() + 0.02);
		h[i]->SetMinimum(0);
		h[i]->Draw();

		TLegend *leg = new TLegend(0.15, 0.6, 0.45, 0.9);
		leg->AddEntry(h[i], "data");
		leg->AddEntry((TF1*)h[i]->FindObject("fitmc"), "total fit");
		leg->AddEntry((TF1*)h[i]->FindObject("mcfcn"), "peak fit");
		leg->AddEntry((TF1*)h[i]->FindObject("bkg"), "bkg fit");
		leg->AddEntry((TF1*)h[i]->FindObject("omgfcn"), "#omega bkg");
		if (h[i]->FindObject("accfcn")) leg->AddEntry((TF1*)h[i]->FindObject("accfcn"), "#acc. bkg");
		if (h[i]->FindObject("best2")) leg->AddEntry((TF1*)h[i]->FindObject("best2"), "#acc2 bkg");
		leg->SetTextSize(0.03);
		leg->Draw();

		if (i == uimanager.output_nbins() - 1)c2->SaveAs(pdfname);
		else if (i % 4 == 3)c1->SaveAs(pdfname);
	}
	c1->SaveAs(pdfname + "]");
}
