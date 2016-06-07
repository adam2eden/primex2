#include "primex.h"
#include "alt_model_fit.h"

int main(int argc, char* argv[])
{
	UImanager uimanager(argc, argv, "alt_model_correct");

	//gStyle->SetOptFit(1);
	gStyle->SetOptStat(1);

	TString inrootdir("./fitroot/");

	//load input root file with rotated mass histograms and sidebands
	TFile *inroot = new TFile(inrootdir + "rotatedmass" + uimanager.output_filename("data") + ".root");
	if (!inroot->IsOpen()) exit(open_err(inrootdir + "rotatedmass" + uimanager.output_filename("data") + ".root"));

	vector<TH1F*> h(uimanager.output_nbins(), NULL), hside(uimanager.output_nbins(), NULL), omega(uimanager.output_nbins(), NULL);
	vector<vector<TH1F*> > mcpeak(uimanager.get_nsigma(), vector<TH1F*>(uimanager.output_nbins(), NULL));

	TVectorD *ratios = (TVectorD*)inroot->Get("ratios");
	for (int i = 0; i < uimanager.output_nbins(); i++){
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

	for (int i = 0; i < uimanager.output_nbins(); i++){
		omega[i] = (TH1F*)inroot->Get(Form("homg_%d", i));
		omega[i]->SetDirectory(0);
	}
	inroot->Close();

	vector<double> mc_sigma(uimanager.output_nbins(), 0.), mc_mean(uimanager.output_nbins(), 0.);
	//load peak fitting model from M.C. if using M.C. to fit rotated mass distribution
    inroot = new TFile(inrootdir + "fitmodel_mc" + uimanager.output_filename("mc") + ".root");
	if (!inroot->IsOpen()) exit(open_err(inrootdir + "fitmodel_mc" + uimanager.output_filename("mc") + ".root"));
    for (int i = 0; i < uimanager.get_nsigma(); i++){
        for (int j = 0; j < uimanager.output_nbins(); j++) {
			mcpeak[i][j] = (TH1F*)inroot->Get(Form("haltinvm_%d_%.2f", j, uimanager.get_sigma_start() + i*uimanager.get_sigma_step()));
            mcpeak[i][j]->SetDirectory(0);
			if (fabs(uimanager.get_sigma_start() + i*uimanager.get_sigma_step() - 1) < 0.005) {
                mc_sigma[j] = mcpeak[i][j]->GetRMS();
                mc_mean[j] = mcpeak[i][j]->GetMean();
            }
        }
    }
    inroot->Close();

	vector<vector<double> > init_par, chi2(uimanager.output_nbins(), vector<double>(uimanager.get_nsigma(), 0.));
	vector<double> sigma(uimanager.get_nsigma(), 0.);

	fitrod fitting(uimanager, init_par);

	vector<vector<vector<double> > > parameters(uimanager.output_nbins(), vector<vector<double> >(uimanager.get_nsigma(), vector<double>(fitting.npar(), 0.))), errors(uimanager.output_nbins(), vector<vector<double> >(uimanager.get_nsigma(), vector<double>(fitting.npar(), 0.)));

	for (int i = 0; i < uimanager.get_nsigma(); i++) {
		sigma[i] = uimanager.get_sigma_start() + i*uimanager.get_sigma_step();		
		for (int j = 0; j < uimanager.output_nbins(); j++){
			fitting.sethist(h[j], hside[j], mcpeak[i][j], omega[j]);
			double angle = (uimanager.get_output_angles()[i] + uimanager.get_output_angles()[i + 1]) / 2;
            parameters[j][i][1] = -100;
			fitting.initialize(parameters[j][i], angle);

			TFitResultPtr r = fitting.fitting(parameters[j][i], angle, "Q");
			for (int k = 0; k < fitting.npar(); k++) errors[j][i][k] = r->Error(k);
            
			chi2[j][i] = r->Chi2() / r->Ndf();
		}
        cout << "--------------------- sigma change: " << sigma[i] << "-------------------------" <<endl;
	}

    cout << "-------------------------- start calculating best sigma ---------------------------" << endl;
	vector<TGraph*> chi2_sigma(uimanager.output_nbins(), NULL);

	int best_par_idx = 0;
	TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
	c1->Divide(2, 2);
	TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
	c2->Divide(2, 2);
	c1->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf["));

	TFile *outroot = new TFile(TString("alt_model_correct" + uimanager.output_filename("fit") + ".root"), "RECREATE");
	vector<double> xangle(uimanager.output_nbins(), 0), best_chi2(uimanager.output_nbins(), 0), sigma_change(uimanager.output_nbins(), 0.), best_sigma(uimanager.output_nbins(), 0), best_sigma_err(uimanager.output_nbins(), 0), mean_change(uimanager.output_nbins(), 0.), best_mean(uimanager.output_nbins(), 0), best_mean_err(uimanager.output_nbins(), 0);
	
	for (int i = 0; i < uimanager.output_nbins(); i++) {
		xangle[i] = (uimanager.get_output_angles()[i] + uimanager.get_output_angles()[i + 1]) / 2;
		TF1 *fitfcn = new TF1("parabolla", "[0]*(x-[1])*(x-[1])+[2]", uimanager.get_sigma_start(), uimanager.get_sigma_start() + uimanager.get_nsigma()*uimanager.get_sigma_step());
		fitfcn->SetParameters(50, 1, 1);
        fitfcn->SetParLimits(0, 0, 200);

		chi2_sigma[i] = new TGraph(uimanager.get_nsigma(), &sigma[0], &chi2[i][0]);
		chi2_sigma[i]->SetName(Form("chi2_sigma_%d", i + 1));
		chi2_sigma[i]->SetTitle(Form("chi2_sigma_%d", i + 1));
		chi2_sigma[i]->Fit("parabolla", "R");
		chi2_sigma[i]->SetTitle(Form("%d: #chi^{2} vs. #sigma change #theta: [%.2f, %.2f]", i + 1, float(uimanager.get_output_angles()[i]), float(uimanager.get_output_angles()[i + 1])));
		chi2_sigma[i]->SetMarkerStyle(20);
		fitfcn->SetLineColor(kRed);
		fitfcn->SetLineWidth(2);

		chi2_sigma[i]->Write();


		if (i / 4 >= uimanager.output_nbins() / 4) c2->cd(i % 4 + 1);
		else c1->cd(i % 4 + 1);
		chi2_sigma[i]->Draw("ap");
		if (i == uimanager.output_nbins() - 1) c2->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf"));
		else if (i % 4 == 3) c1->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf"));

		sigma_change[i] = fitfcn->GetParameter(1);
		best_chi2[i] = fitfcn->GetParameter(2);

		best_sigma[i] = mc_sigma[i] * sigma_change[i];
        best_sigma_err[i] = sqrt(fabs(1/fitfcn->GetParameter(0)/45))*mc_sigma[i];

		best_par_idx = int((sigma_change[i] - uimanager.get_sigma_start()) / uimanager.get_sigma_step());
		if (best_par_idx >= uimanager.get_nsigma()) best_par_idx = uimanager.get_nsigma() - 1;
		else if (best_par_idx < 0) best_par_idx = 0;

        best_mean[i] = mc_mean[i] + parameters[i][best_par_idx][1];
        best_mean_err[i] = mc_sigma[i]/sqrt(parameters[i][best_par_idx][0]);
        mean_change[i] = parameters[i][best_par_idx][1];

		cout << "sigma: " << mc_sigma[i] << " " << sigma_change[i] << " " << best_sigma[i] << " mean: " << mc_mean[i] << " " << parameters[i][best_par_idx][1] << " " << best_mean[i] << endl;

	}

	TGraph *best_chi2_angle = new TGraph(uimanager.output_nbins(), &xangle[0], &best_chi2[0]);
	best_chi2_angle->SetName("best_chi2_angle");

	best_chi2_angle->SetTitle("#Chi^{2} vs. #theta");
	best_chi2_angle->SetMarkerStyle(20);

	TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1200);
	best_chi2_angle->Draw("ap");
	c3->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf"));
	best_chi2_angle->Write();

    vector<double> zeros(uimanager.output_nbins(), 0);
	TGraphErrors *best_sigma_angle = new TGraphErrors(uimanager.output_nbins(), &xangle[0], &best_sigma[0], &zeros[0], &best_sigma_err[0]);
	TGraph *mc_sigma_angle = new TGraph(uimanager.output_nbins(), &xangle[0], &mc_sigma[0]);
	TGraph *sigma_change_angle = new TGraph(uimanager.output_nbins(), &xangle[0], &sigma_change[0]);
    sigma_change_angle->SetName("sigma_change_angle");
    sigma_change_angle->SetTitle("change of #sigma vs. #theta");

	best_sigma_angle->SetMarkerStyle(21);
	best_sigma_angle->SetMarkerColor(kRed);

	mc_sigma_angle->SetMarkerStyle(20);
	mc_sigma_angle->SetMarkerColor(kBlue);
    
    float p[3] = {0.007764 - 0.000953*1.5 + 0.0006715*2.25, -0.000953 + 0.0006715*3, 0.0006715*2};
    TF1 *sigma_trend = new TF1("sigma_trend", Form("%.10f + %.10f*(x-1.5) + %.10f*(x-1.5)*(x-1.5) + [0]*(x-1.5)*(x-1.5)*(x-1.5)", p[0], p[1], p[2]), 1.5, max_angle);
    sigma_trend->SetLineColor(kGreen);
    sigma_trend->SetLineWidth(2);
    best_sigma_angle->Fit("sigma_trend", "R");

	sigma_change_angle->SetMarkerStyle(20);
	sigma_change_angle->SetMarkerColor(kBlue);

	c3->cd();
	TMultiGraph* mg = new TMultiGraph("mg", "M.C. #sigma vs. #theta before and after correction");
	mg->Add(best_sigma_angle);
	mg->Add(mc_sigma_angle);
	mg->Draw("ap");
	c3->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf"));

    c3->cd();
    sigma_change_angle->Draw("ap");
	c3->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf"));

	mg->Write();
    sigma_change_angle->Write();

	TGraphErrors *best_mean_angle = new TGraphErrors(uimanager.output_nbins(), &xangle[0], &best_mean[0], &zeros[0], &best_mean_err[0]);
	TGraphErrors *mc_mean_angle = new TGraphErrors(uimanager.output_nbins(), &xangle[0], &mc_mean[0], &zeros[0], &zeros[0]);
	TGraph *mean_change_angle = new TGraph(uimanager.output_nbins(), &xangle[0], &mean_change[0]);
    mean_change_angle->SetName("mean_change_angle");
    mean_change_angle->SetTitle("change of mean vs. #theta");

	best_mean_angle->SetMarkerStyle(21);
	best_mean_angle->SetMarkerColor(kRed);

    TF1 *mean_trend = new TF1("mean_trend", "pol3", 0, 2.5);
    best_mean_angle->Fit("mean_trend", "R");
    mean_trend->SetLineColor(kGreen);
    mean_trend->SetLineWidth(2);
    
	ofstream output1("alt_model_correct" + uimanager.output_filename("fit") + ".dat");
	ofstream output2("alt_model_correct_parameters" + uimanager.output_filename("fit") + ".dat");
    for (int i = 0 ; i < uimanager.output_nbins(); i++) {
        if (xangle[i] < 1.5) output1 << "1. " << best_chi2[i] << endl;
        else output1 << sigma_trend->Eval(xangle[i]) / mc_sigma[i] << " " << best_chi2[i] << endl;

		best_par_idx = int((sigma_change[i] - uimanager.get_sigma_start()) / uimanager.get_sigma_step());
		if (best_par_idx >= uimanager.get_nsigma()) best_par_idx = uimanager.get_nsigma() - 1;
		else if (best_par_idx < 0) best_par_idx = 0;

		for (int j = 0; j < fitting.npar() - 1; j++) 
            if (j != 1) output2 << parameters[i][best_par_idx][j] << " " << errors[i][best_par_idx][j] << " ";
            else output2 << mean_trend->Eval(xangle[i]) - mc_mean[i] << " " << errors[i][best_par_idx][j] << " ";
		output2 << parameters[i][best_par_idx][fitting.npar() - 1] << " " << errors[i][best_par_idx][fitting.npar() - 1] << " " << chi2[i][best_par_idx] << endl;
    }
	output1.close();
	output2.close();

	mc_mean_angle->SetMarkerStyle(20);
	mc_mean_angle->SetMarkerColor(kBlue);

	mean_change_angle->SetMarkerStyle(20);
	mean_change_angle->SetMarkerColor(kBlue);

	c3->cd();
	TMultiGraph* mg1 = new TMultiGraph("mg1", "M.C. mean vs. #theta before and after correction");
	mg1->Add(best_mean_angle);
	mg1->Add(mc_mean_angle);
	mg1->Draw("ap");
	c3->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf"));

    c3->cd();
    mean_change_angle->Draw("ap");
	c3->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf"));

    mg1->Write();
    mean_change_angle->Write();

	c1->SaveAs(TString("alt_model_correct" + uimanager.output_filename("fit") + ".pdf]"));
	outroot->Close();

	return 0;
}
