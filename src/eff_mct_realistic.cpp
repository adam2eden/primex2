#include "primex.h"

int main (int argc, char* argv[]) {
	UImanager uimanager(argc, argv, "eff_mct");
	TString workdir(getenv("WORKDIR"));

	const int kMax = 100;

	int runnumber, eventid;
	float Etot_cl, Etot;

	int nph, tid[kMax], eid[kMax];
	float phE[kMax];
	int npi0, id1[kMax], id2[kMax], type1[kMax], type2[kMax];
	float invm[kMax], angle[kMax], hx1[kMax], hy1[kMax], hx2[kMax], hy2[kMax], E1[kMax], E2[kMax];

	float MC_angle, MC_E, MC_z, MC_x1, MC_y1, MC_x2, MC_y2, MC_E1, MC_E2;

	TChain *event = new TChain("event");

	event->SetBranchAddress("runnumber", &runnumber);
	event->SetBranchAddress("eventid", &eventid);
	event->SetBranchAddress("Etot_cl",&Etot_cl);
	event->SetBranchAddress("Etot",&Etot);

	event->SetBranchAddress("nph", &nph);
	event->SetBranchAddress("phE", phE);
	event->SetBranchAddress("tid", tid);
	event->SetBranchAddress("eid", eid);
	event->SetBranchAddress("npi0", &npi0);
	event->SetBranchAddress("id1",id1);
	event->SetBranchAddress("id2",id2);
	event->SetBranchAddress("type1",type1);
	event->SetBranchAddress("type2",type2);
	event->SetBranchAddress("hx1",hx1);
	event->SetBranchAddress("hy1",hy1);
	event->SetBranchAddress("hx2",hx2);
	event->SetBranchAddress("hy2",hy2);
	event->SetBranchAddress("E1",E1);
	event->SetBranchAddress("E2",E2);
	event->SetBranchAddress("invm",invm);
	event->SetBranchAddress("angle",angle);
	event->SetBranchAddress("MC_angle",&MC_angle);
	event->SetBranchAddress("MC_E",&MC_E);
	event->SetBranchAddress("MC_z",&MC_z);
	event->SetBranchAddress("MC_x1",&MC_x1);
	event->SetBranchAddress("MC_y1",&MC_y1);
	event->SetBranchAddress("MC_x2",&MC_x2);
	event->SetBranchAddress("MC_y2",&MC_y2);
	event->SetBranchAddress("MC_E1",&MC_E1);
	event->SetBranchAddress("MC_E2",&MC_E2);

	for (int i = 0; i < int(uimanager.get_runs().size()); i++) event->Add(Form("./mc/ge%d.cout.root", uimanager.get_runs()[i]));
    vector<vector<int>> eff_ech_mct(180, vector<int>(540, 0));
    vector<vector<int>> ech_mct(180, vector<int>(540, 0));
    vector<int> eff_mct(540, 0);
    vector<int> mct_total(540, 0);

    TFile outroot("eff_mct_realistic.root", "RECREATE");
    int pass = 0;
    float theta_rec[100]; 
    TTree *ev = new TTree("ev", "one M.C. event");
    ev->Branch("theta_mc", &MC_angle, "theta_mc/F");
    ev->Branch("eid", &eid[0], "eid/I");
    ev->Branch("phE", &phE[0], "phE/F");
    ev->Branch("nrec", &pass, "nrec/I");
    ev->Branch("theta_rec", theta_rec, "theta_rec[nrec]/F");
    

    int nevents = event->GetEntries();
    cout<<"number of events: "<<nevents<<endl;
	for(int i=1;i<=nevents;i++){
		if(i%1000000==0)cout<<i<<endl;
		event->GetEntry(i);
        pass = 0;
        if(MC_angle >= 2.7) continue;
		for(int j=0;j<npi0;j++){
			for(int k=0;k<nph;k++){
                int match1 = 0, match2 = 0;
                if(fabs(hx1[j]-MC_x1)>=4||fabs(hx2[j]-MC_x2)>=4||fabs(hy1[j]-MC_y1)>=4||fabs(hy2[j]-MC_y2)>=4||fabs(1-E1[j]/MC_E1)>=0.3||fabs(1-E2[j]/MC_E2)>=0.3) match1 = 1;
                if(fabs(hx2[j]-MC_x1)>=4||fabs(hx1[j]-MC_x2)>=4||fabs(hy2[j]-MC_y1)>=4||fabs(hy1[j]-MC_y2)>=4||fabs(1-E2[j]/MC_E1)>=0.3||fabs(1-E1[j]/MC_E2)>=0.3) match2 = 1;
                if (match1 && match2) continue;

				if(E1[j]+E2[j]<3.5||E1[j]+E2[j]>8)continue;
                if(E1[j]<0.5 || E2[j]<0.5) continue;

				if(id1[j]<1000||id2[j]<1000)continue;
				int irow1, irow2, icol1, icol2;
				irow1 = (id1[j]-1001)/34+1;
				icol1 = (id1[j]-1001)%34+1;
				if(icol1>=16 && icol1<=19 && irow1>=16 && irow1<=19) continue;
				irow2 = (id2[j]-1001)/34+1;
				icol2 = (id2[j]-1001)%34+1;
				if(icol2>=16 && icol2<=19 && irow2>=16 && irow2<=19) continue;
				int isouterlayer = (icol1 == 1) || (icol1 == 34) || (irow1 == 1) || (irow1 == 34) || (icol2 == 1) || (icol2 == 34) || (irow2 == 1) || (irow2 == 34);
                if(angle[j] < 0 || angle[j] > 3.2) continue;
                theta_rec[pass] = angle[j];
                ++pass;
			}
		}
        int cc = int(MC_angle/0.005);
        mct_total[cc]++;
        ech_mct[eid[0] - 1][cc]++;
/*        if(MC_angle < 0.2) {
            cout << setw(2) << cc << setw(4) << eid[0] << setprecision(3) << fixed << setw(6) << MC_angle << setw(3) << pass << endl; 
        }
*/        if(pass > 0)
        {
            eff_mct[cc]++;
            eff_ech_mct[eid[0] - 1][cc]++;
        }
        ev->Fill();
	}

    ofstream output("eff_ech_mc.dat");
    outroot.cd();
    TH1D heff_mct("eff_mct", "eff_mct", 540, 0, 2.7);
    for(int i = 1; i <= 540; ++i)
    {
        if(mct_total[i - 1] == 0) cout << "zero denominator for M.C. angle: " << i << endl;
        else {
            cout << i << " " << eff_mct[i - 1] << " " << mct_total[i - 1] << " " << eff_mct[i - 1]*1./mct_total[i - 1] << " " << sqrt(mct_total[i - 1])/mct_total[i - 1]*eff_mct[i - 1]*1./mct_total[i - 1] << endl;
            heff_mct.SetBinContent(i, eff_mct[i - 1]*1./mct_total[i - 1]);
            heff_mct.SetBinError(i, sqrt(mct_total[i - 1])/mct_total[i - 1]*eff_mct[i - 1]*1./mct_total[i - 1]);
        }
    }
    for(int i = 1; i <= 180; ++i)
    {
        for(int j = 1; j <= 540; ++j)
        {
            output << setw(3) << i << setw(4) << j << setw(10) << setprecision(7) << fixed << eff_ech_mct[i - 1][j - 1]*1./ech_mct[i - 1][j - 1] << endl;
        }
    }
    heff_mct.Write(); 
    ev->Write();
	outroot.Close();
    output.close();

	return 0;
}
