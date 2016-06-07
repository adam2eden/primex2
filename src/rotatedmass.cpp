#include "primex.h"
#include "stdlib.h"

struct tid_eid_pair {
    int tid;
    int eid;
};

bool operator==(const tid_eid_pair& lhs, const tid_eid_pair& rhs) {
    return lhs.tid == rhs.tid && lhs.eid == rhs.eid;
}

namespace std
{
    template<> struct hash<tid_eid_pair>
    {
        int operator()(tid_eid_pair const& s) const
        {
            return s.tid * 208 + s.eid; 
        }
    };
}

double bw = 0.5;

Double_t DoubleGaussianFit(Double_t *x, Double_t *para) {
    double t1 = (x[0] - para[1])/para[2];
    double t2 = (x[0] - para[1] - para[4])/para[5];
    return para[0]/sqrt(2*3.14159265359)*bw*((1 - para[3])/fabs(para[2])*exp(-0.5*t1*t1) + para[3]/fabs(para[5])*exp(-0.5*t2*t2));
}

Double_t PolyBkg(Double_t *x, Double_t *para){
    double bkg = para[0];
    return bkg;
}

Double_t fitfcn(Double_t *x, Double_t *para){
    return DoubleGaussianFit(x, para) + PolyBkg(x, &para[6]);
}

double btdiff = 0;
vector<double> tdiff_evs;
bool abscmp(const int i, const int j){ return fabs(tdiff_evs[i]-btdiff) < fabs(tdiff_evs[j]-btdiff); };

int main(int argc, char* argv[]){
    TH1::SetDefaultSumw2(kTRUE);
	UImanager uimanager(argc, argv, "rotatedmass");
	
    const int kMax = 100;
	int runnumber, eventid;
	float Etot_cl, Etot;

	int nph, tid[kMax], eid[kMax];
	float phE[kMax], tagmt[kMax], tdiff_tid_offset[kMax], tdiff_hycal_offset[kMax];

    vector<float> tdiff(kMax, 100);
	int npi0, id1[kMax], id2[kMax], type1[kMax], type2[kMax], dtstat[kMax];// ntrigphoton, trigphoton_id[kMax], trigphoton_tid[kMax], trigphoton_eid[kMax];
	float invm[kMax], angle[kMax], hx1[kMax], hy1[kMax], hx2[kMax], hy2[kMax], E1[kMax], E2[kMax], dt[kMax];
	TFile *inroot = new TFile("./fitroot/combine.root");
	TTree *event = (TTree*)inroot->Get("event");

	event->SetBranchStatus("ntagm",0);
	event->SetBranchStatus("tagm_*",0);
	event->SetBranchStatus("nhycal",0);
	event->SetBranchStatus("total_sum",0);
	event->SetBranchStatus("ntrigphoton",0);
	event->SetBranchStatus("trigphoton_*",0);
	event->SetBranchStatus("ncluster",0);
	event->SetBranchStatus("hycalcluster_*",0);

	event->SetBranchAddress("runnumber", &runnumber);
	event->SetBranchAddress("eventid", &eventid);
	event->SetBranchAddress("Etot_cl",&Etot_cl);
	event->SetBranchAddress("Etot",&Etot);

    /*event->SetBranchAddress("ntrigphoton", &ntrigphoton);
    event->SetBranchAddress("trigphoton_tid", trigphoton_tid);
    event->SetBranchAddress("trigphoton_id", trigphoton_id);
    event->SetBranchAddress("trigphoton_eid", trigphoton_eid);
*/
	event->SetBranchAddress("nph", &nph);
	event->SetBranchAddress("phE", phE);
	event->SetBranchAddress("tdiff", &tdiff[0]);
	event->SetBranchAddress("tdiff_tid_offset", tdiff_tid_offset);
	event->SetBranchAddress("tagmt", tagmt);
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
	event->SetBranchAddress("tdiff_hycal_offset", tdiff_hycal_offset);
    event->SetBranchAddress("dtstat", dtstat);
    event->SetBranchAddress("dt", dt);

    TString outrootname("rotatedmass" + uimanager.output_filename("data") + ".root");
    TFile* outroot = new TFile(outrootname, "RECREATE");
	vector<TH1F*> hrotd(uimanager.output_nbins(), NULL), hside(uimanager.output_nbins(), NULL), hrotdbest_bin(uimanager.btcorr_nbins(), NULL), hrotdbest2_bin(uimanager.btcorr_nbins(), NULL), hrotdbest3_bin(uimanager.btcorr_nbins(), NULL), hside_bin(uimanager.btcorr_nbins(), NULL), hrotdbest_out_dt(5, NULL), hrotdbest_out_dt20(uimanager.btcorr_nbins(), NULL);
	TH1F *hrotdbest, *hrotdbest2, *hrotdbest3, *hinvm, *hinvmbest, *hinvmbest2, *helas, *helasbest, *helasbest2, *hsidebest, *hsidebest2, *hinelas, *htdiff, *htdiff_all, *htdiff_best, *hdt, *hrotd_elas, *hrotd_inelas;
    TH1I *multiplicity;
	TH2F *h2d;
    TH1I *nhits_tcs[18];

    for (int i = 0; i < 18; i++)
        nhits_tcs[i] = new TH1I(Form("nhits_tcs_%d", i), Form("T counter %d multiplicity (1 means only best tdiff, >1 means best tdiff and accidentals hits on same counter in one event; hybrid mass and hybrid elasticity cuts applied)", i+1), 10, 1, 11);
	for (int i = 0; i<uimanager.output_nbins(); i++){
		hrotd[i] = new TH1F(Form("hrotd_%d", i), Form("rotated m_{#gamma#gamma} #theta [%3.2f,%3.2f]", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1]), mdiv, -0.2, 0.2);
		hside[i] = new TH1F(Form("hside_%d", i), Form("rotated m_{#gamma#gamma} coincident accidentals #theta [%3.2f,%3.2f]", uimanager.get_output_angles()[i], uimanager.get_output_angles()[i + 1]), mdiv, -0.2, 0.2);
	}
    h2d = new TH2F("h2d", "h2d", 5*mdiv, 0, 0, 5*mdiv, 0, 0);
    hrotdbest = new TH1F("hrotdbest", "rotated m_{#gamma#gamma} w/ best tdiff", mdiv, -0.2, 0.2);
    hrotdbest2 = new TH1F("hrotdbest2", "rotated m_{#gamma#gamma} w/ 2nd best tdiff", mdiv, -0.2, 0.2);
    hrotdbest3 = new TH1F("hrotdbest3", "rotated m_{#gamma#gamma} w/ 3rd best tdiff", mdiv, -0.2, 0.2);
    hsidebest = new TH1F("hsidebest", "rotated m_{#gamma#gamma} sideband w/ best tdiff when only one beam candidate", mdiv, -0.2, 0.2);
    hsidebest2 = new TH1F("hsidebest2", "rotated m_{#gamma#gamma} sideband w/ 2nd best tdiff when number of beam candidates >= 2", mdiv, -0.2, 0.2);
    hinvmbest = new TH1F("hinvmbest", "m_{#gamma#gamma} w/ best tdiff", 240, 0.08, 0.2);
    hinvmbest2 = new TH1F("hinvmbest2", "m_{#gamma#gamma} w/ 2nd best tdiff", 240, 0.08, 0.2);
    helasbest = new TH1F("helasbest", "elasticity w/ best tdiff", 400, 0.8, 1.2);
    helasbest2 = new TH1F("helasbest2", "elasticity w/ 2nd best tdiff", 400, 0.8, 1.2);

    for(int i = 1; i <= 5; ++i) hrotdbest_out_dt[i - 1] = new TH1F(Form("hrotdbest_out_dt_%d", i), Form("rotated m_{#gamma#gamma} w/ best tdiff outside dt window: [%d, %d] ns", -10*i, 10*i), mdiv, -0.2, 0.2);

    htdiff = new TH1F("htdiff","Timing distribution under cuts used in rotatedmass analysis", int(50/bw+0.5),-25,25);
    htdiff_best = new TH1F("tdiff_best", "best tdiff distribution when single beam candidate", int(50/bw+0.5), -25, 25);
    hdt = new TH1F("dt", "timming difference between 2 clusters", 1000, -50, 50);

    for(int i = 1; i <= uimanager.btcorr_nbins(); ++i) {
        hrotdbest_bin[i-1] = new TH1F(Form("hrotdbest_%d", i), Form("rotated m_{#gamma#gamma} w/ best tdiff #theta: [%.2f, %.2f]", uimanager.get_btcorr_angles()[i - 1], uimanager.get_btcorr_angles()[i]), mdiv, -0.2, 0.2);
        hrotdbest2_bin[i-1] = new TH1F(Form("hrotdbest2_%d", i), Form("rotated m_{#gamma#gamma} w/ 2nd best tdiff #theta: [%.2f, %.2f]", uimanager.get_btcorr_angles()[i - 1], uimanager.get_btcorr_angles()[i]), mdiv, -0.2, 0.2);
        hrotdbest3_bin[i-1] = new TH1F(Form("hrotdbest3_%d", i), Form("rotated m_{#gamma#gamma} w/ 3rd best tdiff #theta: [%.2f, %.2f]", uimanager.get_btcorr_angles()[i - 1], uimanager.get_btcorr_angles()[i]), mdiv, -0.2, 0.2);
        hside_bin[i-1] = new TH1F(Form("hside_bin_%d", i), Form("rotated m_{#gamma#gamma} w/ side band best tdiff when only one beam candidate #theta: [%.2f, %.2f]", uimanager.get_btcorr_angles()[i - 1], uimanager.get_btcorr_angles()[i]), mdiv, -0.2, 0.2);
        hrotdbest_out_dt20[i-1] = new TH1F(Form("hrotdbest_out_dt20_%d", i), Form("rotated m_{#gamma#gamma} out of dt window [-20, 20]ns, #theta: [%.2f, %.2f]", uimanager.get_btcorr_angles()[i - 1], uimanager.get_btcorr_angles()[i]), mdiv, -0.2, 0.2);
    }

    htdiff_all = new TH1F("htdiff_all","Timing distribution under cuts used in rotatedmass analysis all tdiff", int(50/bw+0.5),-25,25);
    hrotd_elas = new TH1F("hrotd_elas", Form("rotated m_{#gamma#gamma} elasticity: [%.3f, %.3f]", uimanager.elas_low_cut(), uimanager.elas_high_cut()), mdiv, -0.2, 0.2);
    hrotd_inelas = new TH1F("hrotd_inelas", Form("rotated m_{#gamma#gamma} elas. < %.3f < or elas. > %.3f", uimanager.elas_low_cut(), uimanager.elas_high_cut()), mdiv, -0.2, 0.2);

    hinvm = new TH1F("hinvm", "hinvm", 10*mdiv, 0, 0);
    helas = new TH1F("helas", "helas", 10*mdiv, 0, 0);
    hinelas = new TH1F("hinelas", "hinelas", 10*mdiv, 0, 0);
    multiplicity = new TH1I("multiplicity", "multiplicity", 100, 0, 10);

    hinvm->SetBit(TH1::kCanRebin); 
    helas->SetBit(TH1::kCanRebin); 
    hinelas->SetBit(TH1::kCanRebin); 

	double costheta = 0.70710678118;
	double sintheta = 0.70710678118;
    cout << "--------------------------- start seeking peak center ---------------------------" << endl;
	int nevents = event->GetEntries();
    cout << "events: " << nevents << endl;
    for (int i = 1; i <= nevents; i++) {
        if (i%1000000 == 0) cout << i << endl;
		event->GetEntry(i);
        if(!uimanager.target()) {
            if(64802<=runnumber && runnumber<64854) continue;//cut unstable runs (from single arm)
            if(runnumber<64716 || runnumber>64988) continue;
        } else {
            if(runnumber<65006 || runnumber>65112) continue;
        }
        for(int j = 0;j < npi0; j++){
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

            if(angle[j]>=max_angle)continue;

			for(int k=0;k<nph;k++){
                float elas = (E1[j] + E2[j])/phE[k];
                float tdiff_ev = tdiff[k] - tdiff_tid_offset[k] - tdiff_hycal_offset[j];
                if(int((tid[k]-1)/2)>=18)continue;
				if (-uimanager.get_tdiff_cut() < tdiff_ev && tdiff_ev < uimanager.get_tdiff_cut()){
                    if (!isouterlayer) {
                        hinvm->Fill(invm[j]);
                        helas->Fill(elas);
                    }
                    else if (uimanager.isouter()) {
                        hinvm->Fill(invm[j]);
                        helas->Fill(elas);
                    }
                }
                htdiff->Fill(tdiff_ev);
            }
        }
    }

    TF1 *f1 = new TF1("f1", "gaus", 0.133, 0.139); 
    TF1 *f3 = new TF1("f3", "gaus", 0.99, 1.02); 
    TF1 *f5 = new TF1("f5", fitfcn, -10, 20, 7);

    f5->SetParNames("N_{pk}", "mean1", "#sigma_{1}", "fraction", "distance", "#sigma_{2}");
    f5->SetParameter(0, htdiff->GetEntries());
    f5->SetParameter(1, 0);
    f5->SetParameter(2, 1);
    f5->SetParameter(3, 0.2);
    f5->SetParameter(4, 0);
    f5->SetParameter(5, 2);
    f5->SetParameter(6, 2000);

    hinvm->Fit("f1", "R");
    helas->Fit("f3", "R");
    htdiff->Fit("f5", "R");

    TF1 *f11 = new TF1("f11", "gaus", f1->GetParameter(1)-f1->GetParameter(2), f1->GetParameter(1)+f1->GetParameter(2));
    TF1 *f33 = new TF1("f33", "gaus", f3->GetParameter(1)-f3->GetParameter(2), f3->GetParameter(1)+f3->GetParameter(2));

    hinvm->Fit("f11", "RQ+");
    helas->Fit("f33", "RQ+");

	float da_mass = f11->GetParameter(1);//0.134976;//
	float da_elas = f33->GetParameter(1);//1.;
    btdiff = f5->GetParameter(1) + f5->GetParameter(3)*f5->GetParameter(4);
   
    TF1 *f6 = new TF1("f6", DoubleGaussianFit, min_tdiff, max_tdiff, 6);
    f6->SetParameters(f5->GetParameters());

	double signal_in_cut = f6->Integral(-uimanager.get_tdiff_cut(), uimanager.get_tdiff_cut()) / bw;
	double bkg_in_cut = htdiff->Integral(htdiff->FindBin(-uimanager.get_tdiff_cut()), htdiff->FindBin(uimanager.get_tdiff_cut())) - signal_in_cut;
	double bkg_in_cut1 = f5->GetParameter(6) / bw * 2 * uimanager.get_tdiff_cut();
    double signal_in_sideband = f6->Integral(min_tdiff, max_tdiff)/bw - signal_in_cut;
	double total_in_sideband = htdiff->Integral(htdiff->FindBin(min_tdiff), htdiff->FindBin(-uimanager.get_tdiff_cut())) + htdiff->Integral(htdiff->FindBin(uimanager.get_tdiff_cut()), htdiff->FindBin(max_tdiff));
    double ratio1 = bkg_in_cut/total_in_sideband;
    double ratio2 = signal_in_sideband/signal_in_cut;

    TVectorD v(2);
    if(!uimanager.best_tdiff()){
        v[0] = ratio1;
        v[1] = ratio2;
    }

	if (!uimanager.best_tdiff()) {
        cout << "inv. mass peak (no tran.):   " << da_mass << endl;
        cout << "elas. peak (no tran.):       " << da_elas << endl;
        cout << "tdiff. peak :                " << btdiff << endl;
        cout << "bkg/SB in cut :              " << ratio1 <<endl;
        cout << "signal in SB/signal in cut : " << ratio2 << endl;
        cout << "signal :                     " << signal_in_cut << endl;
        cout << "total in cut :               " << htdiff->Integral(htdiff->FindBin(-uimanager.get_tdiff_cut()), htdiff->FindBin(uimanager.get_tdiff_cut())) << endl;
        cout << "bkg in cut :                 " << bkg_in_cut << endl;
        cout << "bkg in cut (method2) :       " << bkg_in_cut1 << endl;
        cout << "total out of time :          " << total_in_sideband << endl;
    }

    cout << "--------------------------- end of seeking peak center ---------------------------" << endl;
	nevents = event->GetEntries();
    
    double rotd, Tdiff, pangle, _invm, _elas, e0, e1, e2, eph, x1, x2, y1, y2;
    int Tid, multi, done, status, outer, _tc, _ec, _tc_count, _tc_ec_count;
    unordered_map<tid_eid_pair, int> tc_ec_count;
    unordered_map<int, int> tc_count;
    tid_eid_pair _tc_ec;

    TTree *kinematics = new TTree("kinematics", "kinematics");
    kinematics->Branch("rotd", &rotd, "rotd/D");
    kinematics->Branch("invm", &_invm, "invm/D");
    kinematics->Branch("elas", &_elas, "elas/D");
    kinematics->Branch("e0", &e0, "e0/D");
    kinematics->Branch("e1", &e1, "e1/D");
    kinematics->Branch("e2", &e2, "e2/D");
    kinematics->Branch("eph", &eph, "eph/D");
    kinematics->Branch("x1", &x1, "x1/D");
    kinematics->Branch("x2", &x2, "x2/D");
    kinematics->Branch("y1", &y1, "y1/D");
    kinematics->Branch("y2", &y2, "y2/D");
    kinematics->Branch("tdiff", &Tdiff, "tdiff/D");
    kinematics->Branch("Tid", &Tid, "Tid/I");
    kinematics->Branch("outer", &outer, "outer/I");
    kinematics->Branch("angle", &pangle, "angle/D");
    if (uimanager.best_tdiff()) kinematics->Branch("status", &status, "status/I");

    TTree *tc_ec = new TTree("tc_ec", "tc_ec");
	tc_ec->Branch("runnumber", &runnumber, "runnumber/I");
    tc_ec->Branch("tid", &_tc, "tid/I");
    tc_ec->Branch("eid", &_ec, "eid/I");
    tc_ec->Branch("tc_ec_count", &_tc_ec_count, "tc_ec_count/I");

    TTree *tc = new TTree("tc", "tc");
	tc->Branch("runnumber", &runnumber, "runnumber/I");
    tc->Branch("tid", &_tc, "tid/I");
	tc->Branch("tc_count", &_tc_count, "tc_count/I");

    TTree *beam = new TTree("beam", "beam candidate");
    beam->Branch("nph", &nph, "nph/I");
    beam->Branch("phE", phE, "phE[nph]/F");
    beam->Branch("tdiff", &tdiff[0], "tdiff[nph]/F");
    beam->Branch("multi", &multi, "multi/I");

    cout<<"number of events: "<<nevents<<endl;
	for(int i=1;i<=nevents;i++){
		if(i%1000000 == 0)cout<<i<<endl;
		event->GetEntry(i);
        multi = 0; Tid = 99; done = 0;
        e1 = 99; e2 = 99; eph = 99; x1 = 99; x2 = 99; y1 = 99; y2 = 99; Tdiff = 99; pangle = 99; rotd = 99; _elas = 99; _invm = 99; 
        int set = 0;

        tc_ec_count.clear();
        for (int j = 0; j < nph; ++j) {
            _tc_ec.tid = tid[j];
            _tc_ec.eid = eid[j];
            ++tc_ec_count[_tc_ec];
        }

        for (unordered_map<tid_eid_pair, int>::const_iterator iter = tc_ec_count.begin(); iter != tc_ec_count.end(); ++iter) {
            _tc = iter->first.tid;
            _ec = iter->first.eid;
            _tc_ec_count = iter->second;
            tc_ec->Fill();
        }

        tc_count.clear();
        for (int j = 0; j < nph; ++j) {
            ++tc_count[tid[j]];
        }

        for (unordered_map<int, int>::const_iterator iter = tc_count.begin(); iter != tc_count.end(); ++iter) {
            _tc = iter->first;
            _tc_count = iter->second;
            tc->Fill();
        }

		for(int j=0;j<npi0;j++){

            int cc = 0, cc1 = 0;

            for(int m = 0; m < uimanager.output_nbins(); m++) {
				if (uimanager.get_output_angles()[m] < angle[j] && angle[j] <= uimanager.get_output_angles()[m + 1]) {
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
            int irow1, irow2, icol1, icol2;
            irow1 = (id1[j]-1001)/34+1;
            icol1 = (id1[j]-1001)%34+1;
            if(icol1>=16 && icol1<=19 && irow1>=16 && irow1<=19) continue;
            irow2 = (id2[j]-1001)/34+1;
            icol2 = (id2[j]-1001)%34+1;
            if(icol2>=16 && icol2<=19 && irow2>=16 && irow2<=19) continue;
            int isouterlayer = (icol1 == 1) || (icol1 == 34) || (irow1 == 1) || (irow1 == 34) || (icol2 == 1) || (icol2 == 34) || (irow2 == 1) || (irow2 == 34);
            if(!uimanager.target()) {//silicon
                if(64802<=runnumber && runnumber<64854) continue;//cut unstable runs (from single arm)
                if(runnumber<64716 || runnumber>64988) continue;
            } else {
                if(runnumber<65006 || runnumber>65112) continue;
            }
            if(E1[j]+E2[j]<3.5||E1[j]+E2[j]>8)continue;
            if(E1[j]<0.5 || E2[j]<0.5) continue;

            if(id1[j]<1000||id2[j]<1000)continue;

            if(angle[j]>=max_angle)continue;

            if(cc < 0 || cc >= uimanager.output_nbins())continue;

            x1 = hx1[j];
            x2 = hx2[j];
            y1 = hy1[j];
            y2 = hy2[j];
            e1 = E1[j];
            e2 = E2[j];

            vector<int> index, tcs;
            tdiff_evs.clear();
            vector<double> rotds, elasticity, invmass, photon, pi0_e0, pi0_tdiff, pi0_angle;

            for(int k=0;k<nph;k++){
                double tdiff_ev = tdiff[k] - tdiff_tid_offset[k] - tdiff_hycal_offset[j];
                tdiff_evs.push_back(tdiff_ev);
                double elas = (E1[j]+E2[j])/phE[k]/da_elas;
                double cinvm = invm[j]/da_mass;
                rotd = cinvm*costheta-elas*sintheta;

                tcs.push_back((tid[k]-1)/2);
                pi0_tdiff.push_back(tdiff_ev);
                rotds.push_back(rotd); 
                elasticity.push_back(elas);
                invmass.push_back(invm[j]);
                photon.push_back(phE[k]);
                pi0_angle.push_back(angle[j]);
                pi0_e0.push_back(E1[j] + E2[j] - phE[k]);

                if(int((tid[k]-1)/2)>=18)continue;
                if(cinvm*sintheta+elas*costheta<uimanager.get_lelas_cut()) continue;

                h2d->Fill(cinvm, elas);
                if(-0.2 < rotd && rotd < 0.2){
                    set++;
                    Tdiff = tdiff_ev;
                    Tid = tid[k];
                    pangle = angle[j];
                    _invm = invm[j];
                    _elas = (E1[j]+E2[j])/phE[k];
                    e0 = E1[j] + E2[j] - phE[k];
                    eph = phE[k];

                    if(!uimanager.best_tdiff()) {
                        kinematics->Fill();
						if (-uimanager.get_tdiff_cut() < tdiff_ev && tdiff_ev < uimanager.get_tdiff_cut()){
                            if (uimanager.elas_low_cut() < elas && elas < uimanager.elas_high_cut()) {
                                if (!isouterlayer) {
                                    hrotd[cc]->Fill(rotd);
                                    hrotd_elas->Fill(rotd);
                                }
                                else if (uimanager.isouter()) {
                                    hrotd[cc]->Fill(rotd);
                                    hrotd_elas->Fill(rotd);
                                }
                            }
                            else {
                                if (!isouterlayer) hrotd_inelas->Fill(rotd);
                                else if (uimanager.isouter()) hrotd_inelas->Fill(rotd);
                            }
                        }
                        if((min_tdiff<tdiff_ev && tdiff_ev<-6) || (6<tdiff_ev && tdiff_ev<max_tdiff)){
                            if (uimanager.elas_low_cut() < elas && elas < uimanager.elas_high_cut()) {
                                if (!isouterlayer) {
                                    hside[cc]->Fill(rotd);
                                    if(0 <= cc1 && cc1 < uimanager.btcorr_nbins())
                                        hside_bin[cc1]->Fill(rotd);
                                }
                                else if (uimanager.isouter()) {
                                    hside[cc]->Fill(rotd);
                                    if(0 <= cc1 && cc1 < uimanager.btcorr_nbins())
                                        hside_bin[cc1]->Fill(rotd);
                                }
                            }
                        }
                    } else {
                        index.push_back(k);
						if (fabs(tdiff_ev) < uimanager.get_tdiff_cut()) multi++;
                    }

                    if (uimanager.elas_low_cut() < elas && elas < uimanager.elas_high_cut()) {
                        if (!isouterlayer)
                            htdiff_all->Fill(tdiff_ev);
                        else if (uimanager.isouter())
                            htdiff_all->Fill(tdiff_ev);
                    }
                }
            }

            if(uimanager.best_tdiff() && index.size()) {
                sort(index.begin(), index.end(), abscmp);

                for(int k = 0; k < index.size(); k++) {
                    Tdiff = pi0_tdiff[index[k]];
                    Tid = tcs[index[k]];
                    pangle = pi0_angle[index[k]];                    
                    _invm = invmass[index[k]];
                    _elas = elasticity[index[k]];
                    e0 = pi0_e0[index[k]];
                    eph = photon[index[k]];
                    if (k == 0) status = 1;
                    else status = 0;
                    outer = isouterlayer; 
                    kinematics->Fill();
                }
                //2nd best tdiff elas. cut applied to hrotdbest2 hrotdbest2_bin and the 3rd
                //hinvmbest2 helasbest2 are not subject to elas. cuts
                if(index.size()>=2 && fabs(tdiff_evs[index[1]])<uimanager.get_tdiff_cut()) {
                    if (uimanager.elas_low_cut() < elasticity[index[0]] && elasticity[index[0]] < uimanager.elas_high_cut()) {
                        if (!isouterlayer) {
                            hrotdbest2->Fill(rotds[index[1]]);
                            if(0 <= cc1 && cc1 < uimanager.btcorr_nbins())
                                hrotdbest2_bin[cc1]->Fill(rotds[index[1]]);
                        }
                        else if (uimanager.isouter()) {
                            hrotdbest2->Fill(rotds[index[1]]);
                            if(0 <= cc1 && cc1 < uimanager.btcorr_nbins())
                                hrotdbest2_bin[cc1]->Fill(rotds[index[1]]);
                        }
                    }
                    if (!isouterlayer) {
                        hinvmbest2->Fill(invmass[index[1]]);
                        helasbest2->Fill(elasticity[index[1]]);
                    }
                    else if (uimanager.isouter()) {
                        hinvmbest2->Fill(invmass[index[1]]);
                        helasbest2->Fill(elasticity[index[1]]);
                    }
                }

                if(index.size() >= 3 && fabs(tdiff_evs[index[2]]) < uimanager.get_tdiff_cut()) {
                    if (uimanager.elas_low_cut() < elasticity[index[0]] && elasticity[index[0]] < uimanager.elas_high_cut()) {
                        if (!isouterlayer) {
                            hrotdbest3->Fill(rotds[index[1]]);
                            if(0 <= cc1 && cc1 < uimanager.btcorr_nbins())
                                hrotdbest3_bin[cc1]->Fill(rotds[index[1]]);
                        }
                        else if (uimanager.isouter()) {
                            hrotdbest3->Fill(rotds[index[1]]);
                            if(0 <= cc1 && cc1 < uimanager.btcorr_nbins())
                            hrotdbest3_bin[cc1]->Fill(rotds[index[1]]);
                        }
                    }
                }
                //best tdiff elas. cut applied to hrotdbest hrotd[] hrotd_elas hrotd_inelas htdiff_best
                //hinvmbest helasbest are not subjected to elas. cut
				if (fabs(tdiff_evs[index[0]])<uimanager.get_tdiff_cut()){
                    if (!done) {
                        multiplicity->Fill(int(index.size()));
                        done = 1;
                    }
                    if (!isouterlayer) {
                        hinvmbest->Fill(invmass[index[0]]);
                        helasbest->Fill(elasticity[index[0]]);
                    }
                    else if (uimanager.isouter()) {
                        hinvmbest->Fill(invmass[index[0]]);
                        helasbest->Fill(elasticity[index[0]]);
                    }

                    if (uimanager.elas_low_cut() < elasticity[index[0]] && elasticity[index[0]] < uimanager.elas_high_cut()) {
                        if (!isouterlayer) {
                            hrotdbest->Fill(rotds[index[0]]);
                            if(!dtstat[j]) {
                                hdt->Fill(dt[j]);
                                for(int k = 1; k <= 5; ++k)
                                    if(fabs(dt[j]) > k*10) hrotdbest_out_dt[k - 1]->Fill(rotds[index[0]]);
                            }
                            if(0 <= cc1 && cc1 < uimanager.btcorr_nbins()) {
                                hrotdbest_bin[cc1]->Fill(rotds[index[0]]);
                                if(fabs(dt[j]) > 20 && !dtstat[j]) hrotdbest_out_dt20[cc1]->Fill(rotds[index[0]]);
                            }
                            hrotd[cc]->Fill(rotds[index[0]]);
                            hrotd_elas->Fill(rotds[index[0]]);
                        }
                        else if (uimanager.isouter()) {
                            hrotdbest->Fill(rotds[index[0]]);
                            if(!dtstat[j]) {
                                hdt->Fill(dt[j]);
                                for(int k = 1; k <= 5; ++k)
                                    if(fabs(dt[j]) > k*10) hrotdbest_out_dt[k - 1]->Fill(rotds[index[0]]);
                            }
                            if(0 <= cc1 && cc1 < uimanager.btcorr_nbins()) {
                                hrotdbest_bin[cc1]->Fill(rotds[index[0]]);
                                if(fabs(dt[j]) > 20 && !dtstat[j]) hrotdbest_out_dt20[cc1]->Fill(rotds[index[0]]);
                            }
                            hrotd[cc]->Fill(rotds[index[0]]);
                            hrotd_elas->Fill(rotds[index[0]]);
                        }
                    }
                    else {
                        if (!isouterlayer) hrotd_inelas->Fill(rotds[index[0]]);
                        else if (uimanager.isouter()) hrotd_inelas->Fill(rotds[index[0]]);
                    }
                }
                //best tdiff distibution subject to elas. cut
                if (index.size() == 1 && uimanager.elas_low_cut() < elasticity[index[0]] && elasticity[index[0]] < uimanager.elas_high_cut()) {
                    if (!isouterlayer)
                        htdiff_best->Fill(tdiff_evs[index[0]]);
                    else if (uimanager.isouter())
                        htdiff_best->Fill(tdiff_evs[index[0]]);
                }
                //side band with only one beam candidate all subject to elas. cut 
                if((min_tdiff<tdiff_evs[index[0]]&& tdiff_evs[index[0]]<-6) || (6<tdiff_evs[index[0]] && tdiff_evs[index[0]]<max_tdiff)){
                    if (uimanager.elas_low_cut() < elasticity[index[0]] && elasticity[index[0]] < uimanager.elas_high_cut()) {
                        if (!isouterlayer) {
                            if (index.size() == 1) {
                                hside[cc]->Fill(rotds[index[0]]);
                                hsidebest->Fill(rotds[index[0]]);
                                if (0 <= cc1 && cc1 < uimanager.btcorr_nbins()) hside_bin[cc1]->Fill(rotds[index[0]]);
                            }
                            if (index.size() >= 2) hsidebest2->Fill(rotds[index[1]]);
                        }
                        else if (uimanager.isouter()) {
                            if (index.size() == 1) {
                                hside[cc]->Fill(rotds[index[0]]);
                                hsidebest->Fill(rotds[index[0]]);
                                if (0 <= cc1 && cc1 < uimanager.btcorr_nbins()) hside_bin[cc1]->Fill(rotds[index[0]]);
                            }
                            if (index.size() >= 2) hsidebest2->Fill(rotds[index[1]]);
                        }
                    }
                }
            }
            if(set) beam->Fill();
        }
    }

    if(uimanager.best_tdiff()){
        TF1 *f7 = new TF1("f7", fitfcn, -10, 10, 8);

        f7->SetParNames("N_{pk}", "mean1", "#sigma_{1}", "fraction", "distance", "#sigma_{2}");
        f7->SetParameter(0, htdiff_best->GetEntries()/2);
        f7->SetParameter(1, 0);
        f7->SetParameter(2, 1);
        f7->SetParameter(3, 0.2);
        f7->SetParameter(4, 0);
        f7->SetParameter(5, 1);
        f7->SetParameter(6, 2000);
        f7->SetParameter(7, 10);

        f7->SetParLimits(3,0,1);
        f7->SetParLimits(5,0,3);

        htdiff_best->Fit("f7", "R");
        
        TF1 *f8 = new TF1("f8", DoubleGaussianFit, -10, 10, 6);
        f8->SetParameters(f7->GetParameters());
        TF1 *f9 = new TF1("f9", PolyBkg, -10, 10, 2);
        f9->SetParameters(&(f7->GetParameters()[6]));
        htdiff_best->GetListOfFunctions()->Add(f8);
        htdiff_best->GetListOfFunctions()->Add(f9);

		double signal_in_cut_best = f8->Integral(-uimanager.get_tdiff_cut(), uimanager.get_tdiff_cut()) / bw;
		double bkg_in_cut_best = htdiff_best->Integral(htdiff_best->FindBin(-uimanager.get_tdiff_cut()), htdiff_best->FindBin(uimanager.get_tdiff_cut())) - signal_in_cut_best;

		double bkg_in_cut1_best = f7->GetParameter(6) / bw * 2 * uimanager.get_tdiff_cut();
        double signal_in_sideband_best = f8->Integral(min_tdiff, max_tdiff)/bw - signal_in_cut_best;

		double total_in_sideband_best = htdiff_best->Integral(htdiff_best->FindBin(min_tdiff), htdiff_best->FindBin(-uimanager.get_tdiff_cut())) + htdiff_best->Integral(htdiff_best->FindBin(uimanager.get_tdiff_cut()), htdiff_best->FindBin(max_tdiff));
        double ratio1_best = bkg_in_cut_best/total_in_sideband_best;
        double ratio2_best = signal_in_sideband_best/signal_in_cut_best;

        cout << "best tdiff selected" << endl;
        cout << "inv. mass peak (no tran.):   " << da_mass << endl;
        cout << "elas. peak (no tran.):       " << da_elas << endl;
        cout << "tdiff. peak :                " << btdiff << endl;
        cout << "bkg/SB in cut :              " << ratio1_best <<endl;
        cout << "signal in SB/signal in cut : " << ratio2_best << endl;
        cout << "signal :                     " << signal_in_cut_best << endl;
		cout << "total in cut :               " << htdiff_best->Integral(htdiff_best->FindBin(-uimanager.get_tdiff_cut()), htdiff_best->FindBin(uimanager.get_tdiff_cut())) << endl;
        cout << "bkg in cut :                 " << bkg_in_cut_best << endl;
        cout << "bkg in cut (method2) :       " << bkg_in_cut1_best << endl;
        cout << "total out of time :          " << total_in_sideband_best << endl;
        
        v[0] = ratio1_best;
        v[1] = ratio2_best;
        v[0] *= multiplicity->GetMaximum()*1./multiplicity->GetEntries();
        v[1] *= multiplicity->GetMaximum()*1./multiplicity->GetEntries();

        htdiff_best->Write();
    }

    v.Write("ratios");
    h2d->Write();
    htdiff_all->Write();
    htdiff->Write();
    hinvm->Write();
    helas->Write();
    hrotd_elas->Write();
    hrotd_inelas->Write();
    if (uimanager.best_tdiff()) {
        multiplicity->Write();
        hrotdbest->Write();
        hdt->Write();
        for(int i = 0; i < 5; ++i) hrotdbest_out_dt[i]->Write();
        hrotdbest2->Write();
        hrotdbest3->Write();
        for (int i = 0; i < uimanager.btcorr_nbins(); i++) {
            hrotdbest_bin[i]->Write();
            hrotdbest2_bin[i]->Write();
            hrotdbest3_bin[i]->Write();
            hside_bin[i]->Write();
            hrotdbest_out_dt20[i]->Write();
        }
        hsidebest->Write();
        hsidebest2->Write();
        hinvmbest->Write();
        hinvmbest->Write();
        hinvmbest2->Write();
        helasbest->Write();
        helasbest2->Write();
        kinematics->Write();
    }
    beam->Write();

    for(int i = 0; i < uimanager.output_nbins(); i++) {
        hrotd[i]->Write();
        hside[i]->Write();
    }

    tc_ec->Write();
	tc->Write();

	outroot->Close();

	return 0;
}
