#include "primex.h"

using namespace std;

int main(int argc, char* argv[]) {
	UImanager uimanager(argc, argv, "efftable");
//work directory
    TString workdir(getenv("WORKDIR"));
	TString fname;
	if (uimanager.target() == 0) fname = "dfpmc_si.dat";
	else fname = "dfpmc_c12.dat";

    ifstream incs(workdir + "tables/" + fname);
	if (!incs.is_open()) {
		cerr << "cross section file doesn't exist" << endl;
		exit(1);
	}

	int nech, angle, proc;
	vector<vector<vector<double> > > pi0_cs(181, vector<vector<double> >(nmc1 + 1, vector<double>(6, 0)));

	for (int i = 1; i <= 180; i++) {
		for (int j = 1; j <= nmc1; j++) {
			for (int k = 1; k <= 5; k++) {
				incs >> nech >> angle >> proc;
                if (i != nech || j != angle || k != proc) {
					cerr << "corrupted file" << endl;
					exit(1);
				}
                else incs >> pi0_cs[nech][angle][proc];
                cout << nech << " " << angle << " " << proc << pi0_cs[nech][angle][proc] << endl;
			}
		}	
	}

    TString fname0(Form("tmt0_lelas%.2f.dat", uimanger.get_lelas_cut()));
	ifstream tmt0(workdir + "tables/"+fname0);
    TString fname1(Form("tmt1_lelas%.2f.dat", uimanger.get_lelas_cut()));
	ifstream tmt1(workdir + "tables/"+fname1);
	if (!tmt1.is_open() || !tmt0.is_open()) {
		cerr << "efficieny file doesn't exist" << endl;
		exit(1);
	} else {
        cout << "read file " + workdir + "tables/" << fname0 << "and " << fname1 << endl;
    }
	
    int rec1, rec2;
	double eff;
	vector<vector<vector<double> > > efftable0(181, vector<vector<double> >(nrec1 + 1, vector<double>(6, 0)));	
	vector<vector<vector<double> > > efftable1(181, vector<vector<double> >(nrec1 + 1, vector<double>(6, 0)));	
    while (1) {
        tmt0 >> nech >> angle >> rec1 >> rec2; 
        if (tmt0.eof()) break;
        for (int k = rec1; k <= rec2; k++) {
            tmt0 >> eff;
            if (angle <= nmc1 && k <= nrec1)
                for (int l = 1; l <= 5; l++) efftable0[nech][k][l] += eff * pi0_cs[nech][angle][l];
        }
    }
    while (1) {
        tmt1 >> nech >> angle >> rec1 >> rec2; 
        if (tmt1.eof()) break;
        for (int k = rec1; k <= rec2; k++) {
            tmt1 >> eff;
            if (angle <= nmc1 && k <= nrec1)
                for (int l = 1; l <= 5; l++) efftable1[nech][k][l] += eff * pi0_cs[nech][angle][l];
        }
    }

	TString outputname0 = workdir + "/tables/dfp0", outputname1 = workdir + "/tables/dfp1", trail = "";
    if (silicon) trail = "_si.dat";
    else trail = "_c12.dat";
	
	uimanager.setwo(0);
    ofstream output0(outputname0 + uimanager.input_filename("mc") + trail);
	uimanager.setwo(1);
	ofstream output1(outputname1 + uimanager.input_filename("mc") + trail);
	output0 << setprecision(5) << fixed << scientific;
	output1 << setprecision(5) << fixed << scientific;
	for (int i = 1; i <= 180; i++) {
		for (int j = 1; j <= nrec1; j++) {
			for (int k = 1; k <= 5; k++) {
				output0 << setw(3) << i << setw(5) << j << setw(2) << k << setw(13) << efftable0[i][j][k] << endl;
				output1 << setw(3) << i << setw(5) << j << setw(2) << k << setw(13) << efftable1[i][j][k] << endl;
			}
		}
	}

	return 0;
}
