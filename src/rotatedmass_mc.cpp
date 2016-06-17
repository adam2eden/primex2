#include "primex.h"

int main (int argc, char* argv[]) {
	UImanager uimanager(argc, argv, "rotatedmass_mc");

	gStyle->SetOptStat(0);

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

    TString outrootname("fitmodel_mc" + uimanager.output_filename("mc") + ".root"), inrootname;
    if(uimanager.ismc()) {
        outrootname = "rotatedmass" + uimanager.output_filename("data") + ".root";
        inrootname = "fitroot/altfit" + uimanager.output_filename("fit") + ".root";
    }
    TFile* outroot = new TFile(outrootname, "RECREATE");

	vector<vector<TH1F*> > haltinvm(uimanager.output_nbins(), vector<TH1F*>(uimanager.get_nsigma(), NULL));
    vector<TH1F*> haltinvm_mc(uimanager.output_nbins(), NULL), haltinvm_mc_org(uimanager.output_nbins(), NULL), hbkg(uimanager.output_nbins(), NULL), hside(uimanager.output_nbins(), NULL);
	TH1F *hinvm, *helas;

    int nbins = mdiv;
	if (uimanager.get_runs().size() > 10)nbins = 40*nbins;
	for (int i = 0; i < uimanager.output_nbins(); i++) {
        haltinvm_mc[i] = new TH1F(Form("hrotd_%d", i), Form("rot-d m_{#gamma#gamma} w/ rot-d elas. cut #theta [%3.2f,%3.2f]", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1]), nbins, -0.2, 0.2);
        haltinvm_mc_org[i] = new TH1F(Form("hrotd_org_%d", i), Form("rot-d m_{#gamma#gamma} w/ rot-d elas. cut #theta [%3.2f,%3.2f]", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1]), nbins, 0, 0);
		for (int j = 0; j < uimanager.get_nsigma(); j++) {
			haltinvm[i][j] = new TH1F(Form("haltinvm_%d_%.2f", i, uimanager.get_sigma_start() + j*uimanager.get_sigma_step()), Form("rot-d m_{#gamma#gamma} w/ rot-d elas. cut #theta [%3.2f,%3.2f] sigma_adj: %.2f", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1], uimanager.get_sigma_start() + j*uimanager.get_sigma_step()), nbins, -0.2, 0.2);
        }
	}
	hinvm = new TH1F("hinvm", "hinvm", 4*nbins, 0, 0);
	helas = new TH1F("helas", "helas", 4*nbins, 0, 0);

    TFile* inroot;
    vector<float> npi(uimanager.output_nbins(), 0);
    if(uimanager.ismc()) {
        inroot = new TFile(inrootname);
        if(!inroot->IsOpen()) {
            cout << "couldn't find " << inrootname << endl;
            exit(1);
        }
        for(int i = 0; i < uimanager.output_nbins(); i++) hbkg[i] = (TH1F*)inroot->Get(Form("hbkg_%d", i));

        TTree *fitresult = (TTree*)inroot->Get("fitresult");

        const int kMax = 22;
        int npar, NDF, Npfits, index_angle;
        double angles, parameters[kMax], errors[kMax], chi2, Npi0, Npi0_err;

        fitresult->SetBranchAddress("index_angle",&index_angle);
        fitresult->SetBranchAddress("npar",&npar);
        fitresult->SetBranchAddress("angles",&angles);
        fitresult->SetBranchAddress("parameters",parameters);
        fitresult->SetBranchAddress("Npi0",&Npi0);
        fitresult->SetBranchAddress("Npi0_err",&Npi0_err);
        fitresult->SetBranchAddress("errors",errors);
        fitresult->SetBranchAddress("chi2",&chi2);

        for(int i = 0; i < fitresult->GetEntries(); ++i) {
            fitresult->GetEntry(i);
            npi[i] = parameters[0];
        }
    }

	float costheta = 0.70710678118;
	float sintheta = 0.70710678118;
    float mc_mass, mc_elas;
    
	ifstream input("fitmodel.dat");
    int nevents = event->GetEntries();
    cout << nevents << endl;
	if (!input.is_open()) {
		cout << "------------------- can't find fitmodel.dat with peak center --------------------" << endl;
        cout << "--------------------------- start seeking peak center ---------------------------" << endl;
		if (uimanager.get_runs().size() >= 1) {
            for(int i=1;i<=nevents;i++){
                event->GetEntry(i);
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

                        if(angle[j] >= max_angle)continue;

                        float elas = (E1[j]+E2[j])/phE[k];

                        if (!isouterlayer) {
                            hinvm->Fill(invm[j]);
                            helas->Fill(elas);
                        }
                        else if (uimanager.isouter()) {
                            hinvm->Fill(invm[j]);
                            helas->Fill(elas);
                        }
                    }
                }
            }

            TF1 *f1 = new TF1("f1", "gaus", 0.133, 0.139); 
            TF1 *f3 = new TF1("f3", "gaus", 0.99, 1.02); 

            hinvm->Fit("f1", "NRQ");
            helas->Fit("f3", "NRQ");
            
            TF1 *f11 = new TF1("f11", "gaus", f1->GetParameter(1)-f1->GetParameter(2), f1->GetParameter(1)+f1->GetParameter(2));
            TF1 *f33 = new TF1("f33", "gaus", f3->GetParameter(1)-f3->GetParameter(2), f3->GetParameter(1)+f3->GetParameter(2));

            hinvm->Fit("f11", "QR+");
            helas->Fit("f33", "QR+");

            mc_mass = f11->GetParameter(1);//0.134976;//
            mc_elas = f33->GetParameter(1);//1.;

            ofstream output("fitmodel.dat");
            output << mc_mass << " " << mc_elas << endl;
            output.close();
        }
    }
    else {
        input >> mc_mass >> mc_elas;
        input.close();
    }

    cout << "inv. mass peak (no tran.): " << mc_mass << endl;
    cout << "elas. peak (no tran.):     " << mc_elas << endl;
    cout << "--------------------------- end of seeking peak center ---------------------------" << endl;

    TFile *dataroot = new TFile("rotatedmass_angle.root", "RECREATE");
    TTree *rotd = new TTree("rotd", "rotated mass, angle, status of transitional region after cut");
    float rotdm, theta, elas;
    int status;
    TVectorD v(1);
    v[0] = 3.5;
    v.Write("LE1plusE2");
    v[0] = 8;
    v.Write("HE1plusE2");
    v[0] = 0.5;
    v.Write("least_cluster_energy");
    v[0] = 0.2;
    v.Write("rotdm_range");
    v[0] = 1.25;
    v.Write("lelas");

    rotd->Branch("rotd", &rotdm, "rotd/F");
    rotd->Branch("theta", &theta, "theta/F");
    rotd->Branch("elas", &elas, "elas/F");
    rotd->Branch("status", &status, "status/I");

    vector<int> dndt(uimanager.output_nbins(), 0);
	cout<<"number of events: "<<nevents<<endl;
	for(int i=1;i<=nevents;i++){
		if(i%1000000==0)cout<<i<<endl;
		event->GetEntry(i);
        int cc1 = 0;
        for(int m = 0; m < uimanager.output_nbins(); m++) {
            if (uimanager.get_output_angles()[m] < MC_angle && MC_angle <= uimanager.get_output_angles()[m + 1]) {
                cc1 = m;
                break;
            }
        }
        if (cc1 >= 0 && cc1 < uimanager.output_nbins()) dndt[cc1]++;

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
                status = isouterlayer;

				elas = (E1[j]+E2[j])/MC_E/mc_elas;
				float cinvm = invm[j]/mc_mass;

                if(elas < uimanager.elas_low_cut() || elas > uimanager.elas_high_cut()) continue;

				int cc = 0;
                for(int m = 0; m < uimanager.output_nbins(); m++) {
					if (uimanager.get_output_angles()[m] < angle[j] && angle[j] <= uimanager.get_output_angles()[m + 1]) {
                        cc = m;
                        break;
                    }
                }

				if (cc<0 || cc >= uimanager.output_nbins()) continue;

				if(angle[j] > max_angle) continue;

				if (cinvm*costheta - elas*sintheta < uimanager.get_hist_limits()[0] || cinvm*costheta - elas*sintheta > uimanager.get_hist_limits()[1]) continue;

				if (cinvm*sintheta+elas*costheta > uimanager.get_lelas_cut()) {
					for (int m = 0; m < uimanager.get_nsigma(); m++) {
					    if (!isouterlayer) {
                            haltinvm[cc][m]->Fill((cinvm*costheta-elas*sintheta)*(uimanager.get_sigma_start()+m*uimanager.get_sigma_step()));
                            if(m == 0) haltinvm_mc[cc]->Fill(cinvm*costheta-elas*sintheta);
                            if(m == 0) haltinvm_mc_org[cc]->Fill(cinvm*costheta-elas*sintheta);
                        }
                        else if (uimanager.isouter()) {
                            haltinvm[cc][m]->Fill((cinvm*costheta-elas*sintheta)*(uimanager.get_sigma_start()+m*uimanager.get_sigma_step()));
                            if(m == 0) haltinvm_mc[cc]->Fill(cinvm*costheta-elas*sintheta);
                            if(m == 0) haltinvm_mc_org[cc]->Fill(cinvm*costheta-elas*sintheta);
                        }
					}
                    rotdm = cinvm*costheta-elas*sintheta;
                    theta = angle[j];
                    rotd->Fill();
				}
			}
		}
	}

	cout<<"Finished processing input files. Write results."<<endl;
    rotd->Write();
    dataroot->Close();

    TFile f("fitroot/" + outrootname);
    if(!f.IsOpen()) exit(open_err(outrootname));
    TVectorD* ratios = (TVectorD*)f.Get("ratios");
    vector<double> vec = {(*ratios)[0], (*ratios)[1]};
    for(int i = 0; i < uimanager.output_nbins(); ++i) {
        hside[i] = (TH1F*)f.Get(Form("hside_%d", i));
        hside[i]->SetDirectory(0);
    }
    f.Close();

    outroot->cd();
    TVectorD ratios1(2);
    ratios1[0] = vec[0];
    ratios1[1] = vec[1];
    ratios1.Write("ratios");
	
	ofstream output_mc("Npi0_mc.dat");
	ofstream output_MC("Npi0_MC.dat");
	for (int i = 0; i< uimanager.output_nbins(); i++) {
        if(!uimanager.ismc()) {
            for (int j = 0; j < uimanager.get_nsigma(); j++) {
                haltinvm[i][j]->Write();
            }
        }
        else {
            hbkg[i]->Scale(haltinvm_mc[i]->GetEntries()/npi[i]);
			output_mc << haltinvm_mc[i]->GetEntries() << endl;
			output_MC << dndt[i] << endl;
            haltinvm_mc[i]->Add(hbkg[i], 1);
            haltinvm_mc[i]->Write();
            haltinvm_mc_org[i]->Write();
            hside[i]->Write();
        }
    }
	output_mc.close();
    hinvm->Write();
    helas->Write();

	outroot->Close();

	return 0;
}
