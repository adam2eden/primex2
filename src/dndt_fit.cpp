#include "primex.h"

using namespace std;

TH1F *hcoulm = new TH1F("hcoulm", "primakoff", nangle, 0, max_angle);
TH1F *hncohe = new TH1F("hncohe", "nuclear coherent", nangle, 0, max_angle);
TH1F *hcosfiint = new TH1F("hcosfiint", "cosine component of interference", nangle, 0, max_angle);
TH1F *hsinfiint = new TH1F("hsinfiint", "sine component of interference", nangle, 0, max_angle);
TH1F *hninco = new TH1F("hninco", "nuclear incoherent", nangle, 0, max_angle);

Double_t dndt(Double_t *x, Double_t *par){
	Double_t xx = x[0];
	Int_t bin = hcoulm->GetXaxis()->FindBin(xx);
	Double_t coulm = par[0]*hcoulm->GetBinContent(bin);
	Double_t ncohe = par[1]*hncohe->GetBinContent(bin);
	Double_t fiint = TMath::Sqrt(par[0]*par[1])*(TMath::Cos(par[2])*hcosfiint->GetBinContent(bin)+TMath::Sin(par[2])*hsinfiint->GetBinContent(bin));
	Double_t ninco = par[3]*hninco->GetBinContent(bin);
	return coulm+ncohe+fiint+ninco;
}

int main(int argc, char* argv[]){
    UImanager uimanager(argc, argv, "dndt_fit");
//work directory    
    TString workdir(getenv("WORKDIR"));
//***************************set constant ****************************************
    Double_t acdntl_corr = 1., respf_corr = 1.;
	Double_t eff_corr = uimanager.get_adcerr() * uimanager.get_br_corr() * uimanager.get_ta_corr() * acdntl_corr * respf_corr;//	adcerr x tdiff_cut_eff x veto_eff x pi0 decay products abs. (set to 1, since included in current rmat) x bkg_in_fit (not in current rmat)
    if(uimanager.best_tdiff()) eff_corr *= uimanager.get_bit_corr();
    if(uimanager.ismc()) eff_corr = uimanager.get_ta_corr(); 

	Double_t f_l = uimanager.flux() * uimanager.get_tgt_lumi() * eff_corr;

//************ start reading files *******************************
//read tables/eflux.dat [tables/effcor.dat] tables/dfp1_xx.dat or current_foler/dfp1_xx.dat
//get echn flux
	Double_t flw[180];
	ifstream eflux(workdir+"tables/eflux.dat");
	if (!eflux.is_open()) {
		cout << "eflux doesn't exist" << endl;
		exit(1);
 	}
	for(int i=0;i<180;i++)eflux>>flw[i];
//get fitting correction factor
    float effcor[nangle], err[nangle];
    ifstream feffcor;
	bool mc_corr = false;
	if (!uimanager.target()) feffcor.open(workdir+"tables/effcor.dat");
    else feffcor.open(workdir+"tables/effcor_c12.dat");	
	if (feffcor.is_open()) {
		cout << "MC correction exist" << endl;
		mc_corr = true;
		for(int i=0;i<nangle;i++)
			feffcor>>effcor[i]>>err[i];
	} else {
        cout << "fitting w/o M.C. correction" << endl;
    }
//get efficiency table
    TString dfp;
    if (!uimanager.target()) dfp = Form("dfp%d", int(uimanager.isouter())) + uimanager.input_filename("mc") + "_si.dat";
    else dfp = Form("dfp%d", int(uimanager.isouter())) + uimanager.input_filename("mc") + "_c12.dat";

	ifstream profile;
    profile.open(dfp);
    if (profile.is_open()) {
        cout << "use local efficiency table: " << dfp << endl;
    } else {
        //profile.open(workdir+"efficiency2/tables/"+dfp);
        profile.open(workdir+"/tables/"+dfp);
        if(!profile.is_open()) cout << " can't open table " << dfp << endl;
        cout << " use efficiency table in tables: " << dfp << endl;
    }
	
    const int phy = 5;
	const int ech = 180;
	Double_t dfprob[phy][nangle][ech];
	Double_t coulm[nangle], ncohe[nangle], cosfiint[nangle], sinfiint[nangle], ninco[nangle];
	for(int i=0;i<nangle;i++){
		coulm[i]=0;
        ncohe[i]=0;
		cosfiint[i]=0;
		sinfiint[i]=0;
		ninco[i]=0;
	}

	for(int i = 0;i < 180;i++)
		for(int j = 0;j < nangle;j++)
			for(int k = 0;k < 5;k++){
				Double_t f_l_j=f_l;
				if (mc_corr && j < nangle) f_l_j *= effcor[j];
				int ic = 0;
				int jc = 0;
				int kc = 0;
				profile>>ic>>jc>>kc>>dfprob[k][j][i];
                if(ic!=i+1||jc!=j+1||kc!=k+1){
					cout << i << " " << j << " " << k << endl;
					cout << ic << " " << jc << " " << kc << " " << dfprob[k][j][i] << endl;
					cout<<"Bad data given for dfprob"<<endl;
					exit(1);
				}
				if(kc==1)coulm[j]+=flw[i]*dfprob[k][j][i]*f_l_j;
				else if(kc==2)cosfiint[j]+=flw[i]*dfprob[k][j][i]*f_l_j;
				else if(kc==3)ncohe[j]+=flw[i]*dfprob[k][j][i]*f_l_j;
				else if(kc==4)ninco[j]+=flw[i]*dfprob[k][j][i]*f_l_j;
				else sinfiint[j]+=flw[i]*dfprob[k][j][i]*f_l_j;
            }

	for(int i=1;i<=nangle;i++){
		hcoulm->SetBinContent(i,coulm[i-1]);
		hcosfiint->SetBinContent(i,cosfiint[i-1]);
		hncohe->SetBinContent(i,ncohe[i-1]);
		hninco->SetBinContent(i,ninco[i-1]);
		hsinfiint->SetBinContent(i,sinfiint[i-1]);
    }
//read yield
    TString inrootname("pi0alt" + uimanager.input_filename("fit") + ".root");

	TFile *yield = new TFile(inrootname);
    if (!yield->IsOpen()) exit(open_err(inrootname));
	TH1F *yieldhist = (TH1F*)yield->Get("hyield");
	//TH1F *yieldhist = (TH1F*)yield->Get("hyield_mc_nocut");
	yieldhist->SetDirectory(0);
	if (!uimanager.target() && !uimanager.isouter()) yieldhist->SetTitle("#pi^{0} yield (Si, crystal w/o tran.)");
    else if (!uimanager.target()) yieldhist->SetTitle("#pi^{0} yield (Si, crystal w/ tran.)");
    else if (!uimanager.isouter()) yieldhist->SetTitle("#pi^{0} yield (C12, crystal w/o tran.)");
    else yieldhist->SetTitle("#pi^{0} yield (C12, crystal w/ tran.)");
    
    double nyield = 0, nyield_err = 0;
    for(int i = 1; i<=125; i++) {
        nyield += yieldhist->GetBinContent(i);
        nyield_err += yieldhist->GetBinError(i)*yieldhist->GetBinError(i);
    }
    nyield_err = sqrt(nyield_err);
//***************** end reading files **************************************
    
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	gStyle->SetPaperSize(12,24);
	Double_t parameter[4] = {7.9, 1.0, 1.0, 0.8};

	TF1 *ftot = new TF1("ftot", dndt, 0., dndt_fit_range, 4);
	ftot->SetParNames("#Gamma","C1","#phi","C2");

	ftot->SetParameter(0,parameter[0]);
	ftot->SetParameter(1,parameter[1]);
	ftot->SetParameter(2,parameter[2]);
	ftot->SetParameter(3,parameter[3]);

    ftot->SetParLimits(0, 3, 10);
    ftot->SetParLimits(1, 0.3, 1.5);
    ftot->SetParLimits(2, 0, 3.1415936/2);
    ftot->SetParLimits(3, 0., 10.0);

	yieldhist->Fit("ftot", "RMBE0");

	Double_t chi2 = ftot->GetChisquare()/ftot->GetNDF();
	parameter[0]=ftot->GetParameter(0);
	parameter[1]=ftot->GetParameter(1);
	parameter[2]=ftot->GetParameter(2);
	parameter[3]=ftot->GetParameter(3);
	
    Double_t e1 = ftot->GetParError(0);
    Double_t e2 = ftot->GetParError(1);
    Double_t e3 = ftot->GetParError(2);
    Double_t e4 = ftot->GetParError(3);

	//cout<<chi2<<" "<<parameter[0]<<" "<<parameter[1]<<" "<<parameter[2]<<" "<<parameter[3]<<endl;
    TString fname = "fitresult" + uimanager.input_filename("fit") + ".dat";
    ofstream output(fname);
    output<<nyield<<" "<<nyield_err<<" "<<chi2<<" "<<parameter[0]<<" "<<e1<<" "<<parameter[1]<<" "<<e2<<" "<<parameter[2]<<" "<<e3<<" "<<parameter[3]<<" "<<e4<<endl;

	TH1F *fithist = new TH1F("fithist", "fithist", fit_nangle, 0, dndt_fit_range);
	TH1F *fitcoulm = new TH1F("fitcolum", "fitcolum", fit_nangle, 0, dndt_fit_range);
	TH1F *fitncohe = new TH1F("fitncohe", "fitncohe", fit_nangle, 0, dndt_fit_range);
	TH1F *fitit = new TH1F("fitit", "fitit", fit_nangle, 0, dndt_fit_range);
	TH1F *fitninco = new TH1F("fitninco", "fitninco", fit_nangle, 0, dndt_fit_range);

	yieldhist->SetLineWidth(2);
	fithist->SetLineWidth(2);
	fitcoulm->SetLineWidth(2);
	fitncohe->SetLineWidth(2);
	fitit->SetLineWidth(2);
	fitninco->SetLineWidth(2);

	float accfit[fit_nangle], acchist[fit_nangle], acchisterr[fit_nangle], diff[fit_nangle];

	for(int i=1;i<=fit_nangle;i++){
		fithist->SetBinContent(i,ftot->Eval(max_angle/nangle*(i-0.5)));
		fitcoulm->SetBinContent(i,parameter[0]*hcoulm->GetBinContent(i));
		fitncohe->SetBinContent(i,parameter[1]*hncohe->GetBinContent(i));
		fitit->SetBinContent(i,TMath::Sqrt(parameter[0]*parameter[1])*(TMath::Cos(parameter[2])*hcosfiint->GetBinContent(i)+TMath::Sin(parameter[2])*hsinfiint->GetBinContent(i)));
		fitninco->SetBinContent(i,parameter[3]*hninco->GetBinContent(i));
		if (i == 1) {
			accfit[0] = fithist->GetBinContent(1);
			acchist[0] = yieldhist->GetBinContent(1);
			acchisterr[0] = yieldhist->GetBinError(1) * yieldhist->GetBinError(1);
		}
		if (i != fit_nangle) {
			acchist[i] = yieldhist->GetBinContent(i) + acchist[i - 1];
			accfit[i] = fithist->GetBinContent(i) + accfit[i - 1];
			acchisterr[i] = acchisterr[i - 1] + yieldhist->GetBinError(i) * yieldhist->GetBinError(i);
		}
	}

	for(int i = 0; i < fit_nangle; ++i) {
		
		diff[i] = acchist[i] - accfit[i];
		acchisterr[i] = sqrt(acchisterr[i]);
	}

	TString outfilename("fyield" + uimanager.input_filename("fit"));

	float angles[fit_nangle], ex[fit_nangle];
	for(int i=0;i<fit_nangle;i++) {
		angles[i] = 0.02*(i+0.5);
        ex[i] = 0;
    }

	TGraphErrors ge(fit_nangle, angles, diff, ex, acchisterr);
	TGraphErrors ge1(fit_nangle, angles, acchist, ex, acchisterr);
	TGraph ge2(fit_nangle, angles, accfit);

	TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
	c1->SaveAs(outfilename+".pdf[");
	if (mc_corr) {
		TGraphErrors *g1 = new TGraphErrors(fit_nangle,angles,effcor,ex,err);
		g1->SetMarkerStyle(4);
		g1->SetTitle("fitting procedure systematics correction by M.C. (w/ rot-d elas. cut eff.)");
		c1->cd();
		g1->Draw("ap");
		g1->GetXaxis()->SetRangeUser(0,max_angle);	
		g1->SetMarkerStyle(21);
		g1->SetMarkerColor(kBlue);
		g1->SetMarkerSize(1);
		c1->SaveAs(outfilename+".pdf");
	}	
	yieldhist->SetMinimum(0);
	yieldhist->SetMaximum(yieldhist->GetMaximum()*1.05);
	yieldhist->Draw("e1");
	fitcoulm->SetLineColor(kBlue);
	fitcoulm->Draw("sameC");
	fitncohe->SetLineColor(32);
	fitncohe->Draw("sameC");
	fitit->SetLineColor(kBlack);
	fitit->Draw("sameC");
	fitninco->SetLineColor(kYellow);
	fitninco->Draw("sameC");
	fithist->SetLineColor(kRed);
	fithist->Draw("sameC");

	TLegend *leg = new TLegend(0.1,0.7,0.45,0.9);
	leg->SetFillColor(0);
	leg->SetTextSize(0.03);
	leg->AddEntry(fitcoulm,"Primakoff","L");
	leg->AddEntry(fitncohe,"Nuclear Coherent","L");
	leg->AddEntry(fitit,"Interference","L");
	leg->AddEntry(fitninco,"Nuclear Incoherent","L");
	//leg->Draw();

	c1->SaveAs(outfilename+".pdf");

	ge.SetMarkerStyle(20);
	ge.SetMarkerColor(kBlue);
	ge.SetTitle("accumulated yield - fit vs. #theta");
	ge.Draw("ap");
	c1->SaveAs(outfilename + ".pdf");

	ge1.SetMarkerStyle(20);
	ge1.SetMarkerColor(kBlue);
	ge1.SetTitle("accumulated yield and fitting vs. #theta");
	ge1.Draw("ap");
	ge2.SetLineColor(kRed);
	ge2.Draw("sameC");
	c1->SaveAs(outfilename + ".pdf");

	c1->SaveAs(outfilename+".pdf]");

	TFile *ftyd = new TFile(outfilename+".root", "RECREATE");
	yieldhist->Write();
	fithist->Write();
	fitcoulm->Write();
	fitncohe->Write();
	fitit->Write();
	fitninco->Write();
    hcoulm->Write();
    hncohe->Write();
    hninco->Write();
	hcosfiint->Write();
    hsinfiint->Write();
    ftyd->Close();
    
	return 0;
}
