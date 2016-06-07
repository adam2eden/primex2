#include "primex.h"

using namespace std;

int main(int argc, char* argv[]){
	UImanager uimanager(argc, argv, "rotatedmass_omg");
	gStyle->SetOptStat(0);

	const int kMax = 100;

	int runnumber, eventid;
	float Etot_cl, Etot;
	int nph, tid[kMax];
	float phE[kMax];
	int npi0, id1[kMax], id2[kMax], type1[kMax], type2[kMax];
	float invm[kMax], angle[kMax], hx1[kMax], hy1[kMax], hx2[kMax], hy2[kMax], E1[kMax], E2[kMax];

    char* pPath = getenv("WORKDIR");
    if (!pPath) {
        cerr << "Please set env WORKDIR" << endl;
        exit(1);
    }

    TString workdir(pPath);
    TString fname;
    if(!uimanager.target()) fname = "combine_si.root";
    else fname = "combine_c12.root";
    cout << "use " + fname << endl;
    TFile *inroot = new TFile(workdir + "/omega/" + fname);
	TTree *event = (TTree*)inroot->Get("event");
    
	event->SetBranchAddress("runnumber", &runnumber);
	event->SetBranchAddress("eventid", &eventid);
	event->SetBranchAddress("Etot_cl",&Etot_cl);
	event->SetBranchAddress("Etot",&Etot);

	event->SetBranchAddress("nph", &nph);
	event->SetBranchAddress("phE", phE);
	event->SetBranchAddress("tid", tid);
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

    vector<TH1F*> homg(uimanager.output_nbins(), NULL), homg_inelas(uimanager.output_nbins(), NULL), homg_bin(uimanager.output_nbins(), NULL);

	for(int i=0;i<uimanager.output_nbins();i++){
		homg[i] = new TH1F(Form("homg_%d", i), Form("rotated m_{#gamma#gamma} omg bkg. #theta [%3.2f,%3.2f]", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1]), mdiv, -0.2, 0.2);
		homg_inelas[i] = new TH1F(Form("homg_inelas_%d", i), Form("rotated m_{#gamma#gamma} omg bkg. #theta [%3.2f,%3.2f] elas < %.3f or elas > %.3f", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1], uimanager.elas_low_cut(), uimanager.elas_high_cut()), mdiv, -0.2, 0.2);
	}
    for (int i = 0; i < uimanager.btcorr_nbins(); i++) homg_bin[i] = new TH1F(Form("omega_bin_%d", i+1), Form("omega_bin_%d", i+1), mdiv, -0.2, 0.2);
    TH1F homg_elas("homg_elas", Form("omega distribution, #theta [0.00, %3.2f]", max_angle), mdiv, -0.2, 0.2);
    TH1F homg_inelas_total("homg_inelas", Form("inelastic omega distribution, #theta [0.00, %3.2f] elas. < %.3f or elas > %.3f", max_angle, uimanager.elas_low_cut(), uimanager.elas_high_cut()), mdiv, -0.2, 0.2);

	float mc_mass = 0.134875;
	float mc_elas = 1.;
	float costheta = 0.70710678118;
	float sintheta = 0.70710678118;

	int nevents = event->GetEntries();
	cout<<"number of events: "<<nevents<<endl;
	for(int i=1;i<=nevents;i++){
		if(i%1000000==0)cout<<i<<endl;
		event->GetEntry(i);
		for(int j=0;j<npi0;j++){
			for(int k=0;k<nph;k++){
                if(E1[j]<0.5 || E2[j]<0.5) continue;
				if(E1[j]+E2[j]<3.5||E1[j]+E2[j]>8) continue;

				if(id1[j]<1000||id2[j]<1000)continue;
				int irow1, irow2, icol1, icol2;
				irow1 = (id1[j]-1001)/34+1;
				icol1 = (id1[j]-1001)%34+1;
				if(icol1>=16 && icol1<=19 && irow1>=16 && irow1<=19) continue;
				irow2 = (id2[j]-1001)/34+1;
				icol2 = (id2[j]-1001)%34+1;
				if(icol2>=16 && icol2<=19 && irow2>=16 && irow2<=19) continue;
				int isouterlayer = (icol1 == 1) || (icol1 == 34) || (irow1 == 1) || (irow1 == 34) || (icol2 == 1) || (icol2 == 34) || (irow2 == 1) || (irow2 == 34);

				if(angle[j]>max_angle)continue;

				float elas = (E1[j]+E2[j])/phE[k]/mc_elas;
				float cinvm = invm[j]/mc_mass;

				int cc = 0, cc1 = 0;
                for(int m=0; m<uimanager.output_nbins(); m++) {
					if (uimanager.get_output_angles()[m]<angle[j] && angle[j] <= uimanager.get_output_angles()[m + 1]) {
                            cc = m;
                            break;
					}
				}
                for(int m = 0; m < uimanager.btcorr_nbins(); m++) {
                    if (uimanager.get_btcorr_angles()[m] < angle[j] && angle[j] <= uimanager.get_btcorr_angles()[m + 1]) {
                        cc1 = m;
                        break;
                    }
                }

				if(cc<0||cc>=uimanager.output_nbins())continue;

				if(uimanager.get_lelas_cut()<cinvm*sintheta+elas*costheta){
					if(-0.2<cinvm*costheta-elas*sintheta&&cinvm*costheta-elas*sintheta<0.2){
                        if (uimanager.elas_low_cut() < elas && elas < uimanager.elas_high_cut()) {
                            if (!isouterlayer) {
                                homg[cc]->Fill(cinvm*costheta-elas*sintheta);
                                if (0 <= cc1 && cc1 < uimanager.btcorr_nbins()) homg_bin[cc1]->Fill(cinvm*costheta-elas*sintheta);
                            }
                            else if (uimanager.isouter()) {
                                homg[cc]->Fill(cinvm*costheta-elas*sintheta);
                                if (0 <= cc1 && cc1 < uimanager.btcorr_nbins()) homg_bin[cc1]->Fill(cinvm*costheta-elas*sintheta);
                            }
                            homg_elas.Fill(cinvm*costheta-elas*sintheta);
                        }
                        else {
                            if (!isouterlayer) homg_inelas[cc]->Fill(cinvm*costheta-elas*sintheta);
                            else if (uimanager.isouter()) homg_inelas[cc]->Fill(cinvm*costheta-elas*sintheta);
                            homg_inelas_total.Fill(cinvm*costheta-elas*sintheta);
                        }
					}
				}
			}
		}
	}

	TFile* outroot = new TFile(TString("rotatedmass" + uimanager.output_filename("omega") + ".root"), "RECREATE");
	for(int i=0;i<uimanager.output_nbins();i++){
		homg[i]->Write();
		homg_inelas[i]->Write();
	}
    for (int i = 0; i < uimanager.btcorr_nbins(); i++) homg_bin[i]->Write();
    homg_elas.Write();
    homg_inelas_total.Write();

	outroot->Close();

	return 0;
}
