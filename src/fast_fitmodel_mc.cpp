#include "primex.h"

int main(int argc, char* argv[]) {
    UImanager uimanager(argc, argv, "fast_fitmodel_mc");

    TFile *inroot = new TFile("./fitroot/rotatedmass_angle.root");
    if (!inroot) exit(open_err(string("rotatedmass_angle.root")));
    TTree *rotd = (TTree*)inroot->Get("rotd");
    float rotdm, theta;
    int status;

    rotd->SetBranchAddress("rotd", &rotdm);
    rotd->SetBranchAddress("theta", &theta);
    rotd->SetBranchAddress("status", &status);

	vector<vector<TH1F*> > haltinvm(uimanager.output_nbins(), vector<TH1F*>(uimanager.get_nsigma(), NULL));

    vector<double> sigma_change(nangle, 0.);
    int nbins = 4000;
    for (int i = 0; i < uimanager.output_nbins(); i++) {
        for (int j = 0; j < uimanager.get_nsigma(); j++) {
			haltinvm[i][j] = new TH1F(Form("haltinvm_%d_%.2f", i, uimanager.get_sigma_start() + j*uimanager.get_sigma_step()), Form("rot-d m_{#gamma#gamma} w/ rot-d elas. cut #theta [%3.2f,%3.2f] sigma_adj: %.2f", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1], uimanager.get_sigma_start() + j*uimanager.get_sigma_step()), nbins, -0.2, 0.2);
            haltinvm[i][j]->SetDirectory(0);
        }
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

        for (int j = 0; j < uimanager.get_nsigma(); j++) {
			if (uimanager.isouter()) {
				haltinvm[cc][j]->Fill(rotdm*(uimanager.get_sigma_start() + j*uimanager.get_sigma_step()));
			}
			else if (!status) haltinvm[cc][j]->Fill(rotdm*(uimanager.get_sigma_start() + j*uimanager.get_sigma_step()));
        }
    }
    inroot->Close();

    TFile *outroot = new TFile(TString("fitmodel_mc" + uimanager.output_filename("mc") + ".root"), "RECREATE");
    outroot->cd();
    for (int i = 0; i < uimanager.output_nbins(); i++) {
        for (int j = 0; j < uimanager.get_nsigma(); j++) {
            haltinvm[i][j]->Write();
        }
    }
    outroot->Close();

    return 0;
}
