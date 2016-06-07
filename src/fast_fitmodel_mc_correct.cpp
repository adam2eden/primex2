#include "primex.h"

int main(int argc, char* argv[]) {
    UImanager uimanager(argc, argv, "fast_fitmodel_mc_correct");

    TFile *inroot = new TFile("./fitroot/rotatedmass_angle.root");
    if (!inroot) exit(open_err(string("rotatedmass_angle.root")));
    TTree *rotd = (TTree*)inroot->Get("rotd");
    float rotdm, theta, elas;
    int status;

    rotd->SetBranchAddress("rotd", &rotdm);
    rotd->SetBranchAddress("theta", &theta);
    rotd->SetBranchAddress("elas", &elas);
    rotd->SetBranchAddress("status", &status);

	vector<TH1F*> haltinvm(uimanager.output_nbins(), NULL), haltinvm_inelas(uimanager.output_nbins(), NULL);

	vector<double> sigma_change(uimanager.output_nbins(), 1.);
	ifstream input_sigma_change("alt_model_correct" + uimanager.input_filename("fit") + ".dat");
	if (!input_sigma_change.is_open()) {
        open_err("alt_model_correct" + uimanager.input_filename("fit") + ".dat");
        cout << "Use original sigma instead." << endl;
    }
    else {
        cout << "Use alt_model_correct" + uimanager.input_filename("fit") + ".dat" << endl;
        for (int i = 0; i < uimanager.input_nbins(); i++) {
            double temp, chi2;
            input_sigma_change >> temp >> chi2;
            if (uimanager.input_nbins() == uimanager.output_nbins()) sigma_change[i] = temp; 
            else {
                int start = int(uimanager.get_input_angles()[i] / (max_angle / nangle) + 0.5);
                int end = int(uimanager.get_input_angles()[i + 1] / (max_angle / nangle) + 0.5);
                for (int j = start; j < end; j++) sigma_change[j] = temp;
            }
        }
    }

    int nbins = 4000;
    for (int i = 0; i < uimanager.output_nbins(); i++) {
        haltinvm[i] = new TH1F(Form("haltinvm_%d", i), Form("rot-d m_{#gamma#gamma} w/ rot-d elas. cut #theta [%3.2f,%3.2f]", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1]), nbins, -0.2, 0.2);
        haltinvm[i]->SetDirectory(0);
        haltinvm_inelas[i] = new TH1F(Form("haltinvm_inelas_%d", i), Form("rot-d m_{#gamma#gamma} w/ rot-d elas. cut #theta [%3.2f,%3.2f] elas < %.3f or elas > %.3f", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1], uimanager.elas_low_cut(), uimanager.elas_high_cut()), nbins, -0.2, 0.2);
        haltinvm_inelas[i]->SetDirectory(0);
    }

    int ncandidates = rotd->GetEntries();
    cout << "number of pi0 candidates: " << ncandidates << endl;
	int step = 1;
    for (int i = 0; i < ncandidates; i++) {
		if (int(1.*i / ncandidates * 10) >= step) {
			cout << step * 10 << " %" << endl;
			step++;
		}
        rotd->GetEntry(i);
        
		int cc = 0;
        for(int m=0; m<nangle; m++)
		if (uimanager.get_output_angles()[m] < theta && theta <= uimanager.get_output_angles()[m + 1]) {
                cc = m;
                break;
            }

        if (cc >= uimanager.output_nbins()) continue;
        
        if (uimanager.elas_low_cut() < elas && elas < uimanager.elas_high_cut()) {
            if (uimanager.isouter()) haltinvm[cc]->Fill(rotdm*sigma_change[cc]);
            else if (!status) haltinvm[cc]->Fill(rotdm*sigma_change[cc]);
        }
        else {
            if (uimanager.isouter()) haltinvm_inelas[cc]->Fill(rotdm*sigma_change[cc]);
            else if (!status) haltinvm_inelas[cc]->Fill(rotdm*sigma_change[cc]);
        }
    }
    inroot->Close();

    TFile *outroot = new TFile(TString("fitmodel_mc" + uimanager.output_filename("fit") + ".root"), "RECREATE");
    outroot->cd();
    for (int i = 0; i < uimanager.output_nbins(); i++) {
        haltinvm[i]->Write();
        haltinvm_inelas[i]->Write();
    }
    outroot->Close();

    return 0;
}
