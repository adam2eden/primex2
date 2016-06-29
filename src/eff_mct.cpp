#include <primex.h>

int main(int argc, char* argv[]){
    vector<double> eid_dist = {0.0024208768, 0.0013084378, 0.0080143990, 0.0028857995, 0.0063417418, 0.0026269207, 0.0059270122, 0.0030106680, 0.0063838612, 0.0031747853, 0.0068374140, 0.0032393836, 0.0057597513, 0.0028512922, 0.0055873314, 0.0036027598, 0.0070705900, 0.0041495445, 0.0061263916, 0.0035052080, 0.0062953467, 0.0037479553, 0.0071386490, 0.0034913690, 0.0053078084, 0.0035455046, 0.0084149965, 0.0034658228, 0.0063346998, 0.0037505208, 0.0087057448, 0.0023106382, 0.0048839772, 0.0022146080, 0.0086812465, 0.0022498460, 0.0049359289, 0.0026215449, 0.0097910838, 0.0030694448, 0.0057011592, 0.0023492927, 0.0066916243, 0.0033049374, 0.0090496372, 0.0041227858, 0.0062285564, 0.0035706614, 0.0086107986, 0.0039306612, 0.0071202812, 0.0039015458, 0.0076017892, 0.0040982112, 0.0086700492, 0.0028680902, 0.0054940915, 0.0033979283, 0.0100612000, 0.0050620340, 0.0068545612, 0.0048460092, 0.0088548628, 0.0038326959, 0.0060979868, 0.0040158675, 0.0082774170, 0.0066623804, 0.0082180099, 0.0057133803, 0.0067992773, 0.0043726112, 0.0086418612, 0.0037498824, 0.0063548702, 0.0044739530, 0.0098256272, 0.0051580481, 0.0066690851, 0.0036691487, 0.0080661980, 0.0040915626, 0.0106631361, 0.0040991265, 0.0061937359, 0.0026620984, 0.0073910720, 0.0036334209, 0.0100982487, 0.0030275222, 0.0060826101, 0.0029543403, 0.0093005265, 0.0025709382, 0.0061023870, 0.0024706965, 0.0084049675, 0.0034990533, 0.0091763726, 0.0055113271, 0.0095589677, 0.0041280773, 0.0055169117, 0.0039466803, 0.0107112255, 0.0050447502, 0.0074413133, 0.0044523413, 0.0070669365, 0.0055834491, 0.0104544826, 0.0037378741, 0.0078883942, 0.0060580595, 0.0090152022, 0.0046068552, 0.0058981738, 0.0035763625, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0108579145, 0.0049906186, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0096771154, 0.0020747040, 0.0065395393, 0.0034850938, 0.0116043957, 0.0024547296, 0.0075223201, 0.0036390055, 0.0089209425, 0.0050220103, 0.0091729961, 0.0038423475, 0.0070865167, 0.0040228251, 0.0096624614, 0.0049473269, 0.0108678150, 0.0042415759, 0.0066914476, 0.0034631168, 0.0093252136, 0.0072740644, 0.0100549450, 0.0055336092, 0.0078818702, 0.0054673729, 0.0087701664, 0.0036696746, 0.0109633875, 0.0025023653, 0.0050732393, 0.0021195253, 0.0102841584, 0.0057396692, 0.0105363124, 0.0059488728, 0.0146054764, 0.0000000000, 0.0000000000, 0.0000000000, 0.0167939520, 0.0061681174, 0.0090204535, 0.0060559558, 0.0105734373, 0.0058989005, 0.0069586290, 0.0051160933, 0.0081325989, 0.0010565207, 0.0002536957, 0.0000457327};
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
	
	int nmc = 640, nrec = 160;
	for (int i = 0; i < int(uimanager.get_runs().size()); i++) event->Add(Form("./mc/ge%d.cout.root", uimanager.get_runs()[i]));
    vector<vector<int>> eff_ech_mct(180, vector<int>(nmc, 0));
    vector<vector<int>> ech_mct(180, vector<int>(nmc, 0));
    vector<double> eff_mct(nmc, 0);
    vector<double> mct_total(nmc, 0);

	TFile outroot("eff_mct.root", "RECREATE");
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
        if(MC_angle >= 3.2) continue;
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
        //mct_total[cc]++;
        ech_mct[eid[0] - 1][cc]++;
        if(pass > 0)
        {
            //eff_mct[cc]++;
            eff_ech_mct[eid[0] - 1][cc]++;
        }
        ev->Fill();
	}
	ev->Write();
	outroot.Close();
	
    TFile outroot1("eff_mct_hist.root", "RECREATE");
    TH1D heff_mct("eff_mct", "eff_mct", nmc, 0, 3.2);
	TH1D heff_mct1("eff_mct1", "eff_mct normalized by echn dist.", nmc, 0, 3.2);
    for(int i = 1; i <= nmc; ++i)
    {
		if(mct_total[i - 1] == 0) cout << "zero denominator for M.C. angle: " << i << endl;
		else
		{
			heff_mct.SetBinContent(i, eff_mct[i - 1]/mct_total[i - 1]);
			heff_mct.SetBinError(i, sqrt(mct_total[i - 1])/mct_total[i - 1]*eff_mct[i - 1]*1./mct_total[i - 1]);
			for(int j = 0; j < 180; ++j)
			{
				double scale = eid_dist[j];
				mct_total[i - 1] += ech_mct[j][i - 1] * scale;
				eff_mct[i - 1] += eff_ech_mct[j][i - 1] * scale;
			}
			if(mct_total[i - 1] == 0) cout << "zero denominator for M.C. angle: " << i << endl;
			heff_mct1.SetBinContent(i, eff_mct[i - 1]/mct_total[i - 1]);
			heff_mct1.SetBinError(i, sqrt(mct_total[i - 1])/mct_total[i - 1]*eff_mct[i - 1]*1./mct_total[i - 1]);
		}
    }

	heff_mct.Write();
	heff_mct1.Write(); 
    outroot1.Close();
	
	return 0;
}
