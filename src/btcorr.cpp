#include "primex.h"
#include "alt_model_fit.h"

using namespace std;

int main(int argc, char* argv[])
{
	UImanager uimanager(argc, argv, "btcorr");
    uimanager.set_btc(0);

	TFile inroot(TString("./fitroot/rotatedmass" + uimanager.input_filename("data") + ".root"));
    if (!inroot.IsOpen()) exit(open_err("./fitroot/rotatedmass" + uimanager.input_filename("data") + ".root"));
	TH1F* hrotdbest((TH1F*)inroot.Get("hrotdbest"));
	TH1F* hrotdbest2((TH1F*)inroot.Get("hrotdbest2"));
	TH1F* hrotdbest3((TH1F*)inroot.Get("hrotdbest3"));
	TH1F* hsidebest((TH1F*)inroot.Get("hsidebest"));
	TH1F* hsidebest2((TH1F*)inroot.Get("hsidebest2"));
    TH1F* htdiff_best((TH1F*)inroot.Get("tdiff_best"));
	hrotdbest->SetDirectory(0);
	hrotdbest2->SetDirectory(0);
	hrotdbest3->SetDirectory(0);
	hsidebest->SetDirectory(0);
	hsidebest2->SetDirectory(0);
    htdiff_best->SetDirectory(0);
    
    TF1 *tdiff_best((TF1*)htdiff_best->GetListOfFunctions()->At(1));
    float ratio = tdiff_best->Integral(-uimanager.get_tdiff_cut(), uimanager.get_tdiff_cut())/tdiff_best->GetParameter(0);

    vector<TH1F*> hrotdbest_bin(uimanager.btcorr_nbins(), NULL), hrotdbest2_bin(uimanager.btcorr_nbins(), NULL);
    for (int i = 0; i < uimanager.btcorr_nbins(); ++i) {
        hrotdbest_bin[i] = (TH1F*)inroot.Get(Form("hrotdbest_%d", i + 1));
        hrotdbest2_bin[i] = (TH1F*)inroot.Get(Form("hrotdbest2_%d", i + 1));
        hrotdbest_bin[i]->SetDirectory(0);
        hrotdbest2_bin[i]->SetDirectory(0);
    }

	inroot.Close();

    vector<TH1F*> omega_bin(uimanager.btcorr_nbins(), NULL);
	TFile inroot1(TString("./fitroot/rotatedmass" + uimanager.input_filename("omega") + ".root"));
    if (!inroot1.IsOpen()) exit(open_err("./fitroot/rotatedmass" + uimanager.input_filename("omega") + ".root"));
	TH1F* omega((TH1F*)inroot1.Get("homg_elas"));
	omega->SetDirectory(0);
    omega->Scale(ratio);
    for (int i = 0; i < uimanager.btcorr_nbins(); ++i) {
        omega_bin[i] = (TH1F*)inroot1.Get(Form("omega_bin_%d", i+1));
        omega_bin[i]->Scale(ratio);
    }
	inroot.Close();
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector<vector<double>> init_par;	
	fitrod fitting(uimanager, init_par);
	vector<double> parameters(fitting.npar(), 0.);
	fitting.sethist(hrotdbest, hsidebest, hrotdbest, omega, NULL);
	fitting.initialize(parameters, 0.);
	TFitResultPtr r = fitting.fitting(parameters, 0., "");
    
    uimanager.set_btc(2);
    uimanager.set_poly_bkg(0);
    uimanager.set_sub_omg(0);
    uimanager.set_fit_range(-0.2, 0.2);
    
    fitting.set_uimanager(uimanager);
    fitting.define_par(1);

	vector<double> parameters1(fitting.npar(), 0.);

    for (int i = 1; i < 12; ++i) parameters1[i] = parameters[i];
    parameters1[0] = hrotdbest2->GetEntries()/2;
    parameters1[12] = hrotdbest2->GetEntries()/2;
    parameters1[13] = 0;
    parameters1[14] = 0.08;
    parameters1[15] = 0.5;
    parameters1[16] = 0;
    parameters1[17] = 0.08;

	fitting.sethist(hrotdbest2, hsidebest, hrotdbest, omega);
	fitting.fitting(parameters1, 0., "");
    
    for (int i = 1; i < 12; ++i) parameters1[i] = parameters[i];
    parameters1[0] = hrotdbest3->GetEntries()/2;
    parameters1[12] = hrotdbest3->GetEntries()/2;
    parameters1[13] = 0;
    parameters1[14] = 0.08;
    parameters1[15] = 0.5;
    parameters1[16] = 0;
    parameters1[17] = 0.08;

	fitting.sethist(hrotdbest3, hsidebest, hrotdbest, omega);
	fitting.fitting(parameters1, 0., "");

    ofstream output("btcorr.dat");
    for (int i = 0; i < uimanager.btcorr_nbins(); ++i) {
        uimanager.set_btc(0);
        uimanager.set_poly_bkg(1);
        uimanager.set_sub_omg(1);
        uimanager.set_fit_range(-0.1, 0.1);

        fitting.set_uimanager(uimanager);
        fitting.define_par(0);

        fitting.sethist(hrotdbest_bin[i], hsidebest, hrotdbest, omega_bin[i]);
        fitting.initialize(parameters, 0.);
        r = fitting.fitting(parameters, 0., "");
        
        uimanager.set_btc(2);
        uimanager.set_poly_bkg(0);
        uimanager.set_sub_omg(0);
        uimanager.set_fit_range(-0.2, 0.2);
        
        fitting.set_uimanager(uimanager);
        fitting.define_par(1);

        for (int j = 1; j < 12; ++j) parameters1[j] = parameters[j];
        parameters1[0] = hrotdbest2->GetEntries()/2;
        parameters1[12] = hrotdbest2->GetEntries()/2;
        parameters1[13] = 0;
        parameters1[14] = 0.08;
        parameters1[15] = 0.5;
        parameters1[16] = 0;
        parameters1[17] = 0.08;

        fitting.sethist(hrotdbest2_bin[i], hsidebest, hrotdbest, omega_bin[i]);
        fitting.fitting(parameters1, 0., "");

        TF1 bkg("bkg", Form("[0]*%f*([3]/[2]*exp(-(x-[1])*(x-[1])/2/[2]/[2])+(1-[3])/[5]*exp(-(x-[1]-[4])*(x-[1]-[4])/2/[5]/[5]))", 1/sqrt(2*3.14159265359)), uimanager.get_hist_limits()[0], uimanager.get_hist_limits()[1]);
        bkg.SetParameters(parameters1[12], parameters1[13], parameters1[14], parameters1[15], parameters1[16], parameters1[17]);
        double total_bkg = bkg.Integral(uimanager.get_hist_limits()[0], uimanager.get_hist_limits()[1]);

        output << parameters1[0]/parameters[0]*parameters1[12]/total_bkg << " ";
        for (int j = 13; j < 17; j++) output << parameters1[j] << " ";
        output << parameters1[17] << endl;
    }
    output.close();

	TFile outroot(TString("btcorr" + uimanager.input_filename("fit") + ".root"), "RECREATE");
	hrotdbest->Write();
	hrotdbest2->Write();
	hrotdbest3->Write();
    for (int i = 0; i < uimanager.btcorr_nbins(); i++) {
        hrotdbest_bin[i]->Write();
        hrotdbest2_bin[i]->Write();
    }
	outroot.Close();

    return 0;
}
