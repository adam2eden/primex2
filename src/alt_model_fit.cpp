#include "primex.h"
#include "alt_model_fit.h"

fitrod::fitrod(const UImanager& uimanager, vector<vector<double> >& input_par)
:_npar(22), user_defined(false), uimanager(uimanager), h(NULL), accidentals(NULL), mcpk(NULL), omega(NULL)
{
	init_par = input_par;
}

void fitrod::define_par(bool input) {
    user_defined = input;
}

void fitrod::setinitpar(vector<vector<double> >& input_par) {
    init_par = input_par;
}

void fitrod::sethist(TH1F* _h, TH1F* _accidentals, TH1F* _mcpk, TH1F *_omega, TH1F *_hbkg) {
	h = _h;
	accidentals = _accidentals;
	mcpk = _mcpk;
	omega = _omega;
	hbkg = _hbkg;
}

void fitrod::initialize(vector<double>& par, double angle) {
	if (int(par.size()) != _npar) exit(input_err());

	if (uimanager.btdiff_correction() == 0) {
        par[14] = 0;
        par[15] = 0;
        par[16] = 0.006;
        par[17] = 0.2;
        par[18] = 0;
        par[19] = 0.01;
    }
    else if (uimanager.btdiff_correction() == 1) {// use the correction from program btcorr
        int cc = 0;
        for (int i = 0; i < uimanager.btcorr_nbins(); ++i)
            if ( uimanager.get_btcorr_angles()[i] <= angle && angle < uimanager.get_btcorr_angles()[i + 1])
                cc = i;
        for (int i = 0; i < 6; ++i) par[i + 14] = uimanager.get_btcorr_pars()[cc][i];
	}
	else if (uimanager.btdiff_correction() == 2) {// calculate the correction
		par[14] = h->GetEntries() / 2;
		par[15] = 0;
		par[16] = 0.006;
		par[17] = 0.2;
		par[18] = 0;
		par[19] = 0.01;
	}

	if (uimanager.get_method() == 1 && init_par.size() == 0) {
		par[0] = h->GetEntries()/2.;
		if (!uimanager.target()) {//for silicon
			par[7] = 1;
			par[8] = 5;
			par[9] = 1;
			par[10] = 0;
			par[11] = 0;
			par[12] = 0;
		}
		else { // for carbon
			par[7] = 3.08;
			par[8] = 0.05;
			par[9] = 0.0024;
			par[10] = 0;
			par[11] = 0;
			par[12] = 0;
		}
		if (par[1] < -99.) {
            par[1] = find_mean(par, angle);
            if(fabs(par[1]) > 0.01) par[1] = 0;
        }
        else par[1] = 0.;
	}
	else if (uimanager.get_method() == 2 && init_par.size() == 0) {
		par[0] = h->GetEntries() / 3.;
		par[1] = 0;
		par[2] = 0.006;
		par[3] = 0.2;
		par[4] = 0;
		par[5] = 0.01;
		if (uimanager.subtract_acc()) par[6] = accidentals->GetEntries();
		else par[6] = 0;
		par[7] = 1;
		par[8] = 5;
		par[9] = 1;
		par[10] = 0;
        par[11] = 0;
        par[12] = 0;
	}
	else if (init_par.size()) {
		int iangle = 0;
		for (int i = 0; i < uimanager.input_nbins(); i++) {
			if (uimanager.get_input_angles()[i] < angle && angle < uimanager.get_input_angles()[i + 1]) {
                iangle = i;
                break;
            }
		}
		par = init_par[iangle];
	}

	if (uimanager.subtract_acc()) par[6] = accidentals->GetEntries();
	else par[6] = 0;

	par[13] = uimanager.flux() / uimanager.flux_omega()*omega->GetEntries();
}

void fitrod::calc_Npi0(vector<double>& par, double& Npi0, double& Npi0_err) {
	if (int(par.size()) != _npar) exit(input_err());
	Npi0 = 0; Npi0_err = 0;

    float bw = h->GetBinWidth(1);
    float x1 = -0.05, x2 = 0.05, ratio = 1;

	int lbin = h->FindBin(x1);
	int hbin = h->FindBin(x2);
	double lratio = (h->GetBinLowEdge(lbin) + h->GetBinWidth(lbin) - x1) / h->GetBinWidth(lbin);
	double rratio = (x2 - h->GetBinLowEdge(hbin)) / h->GetBinWidth(hbin);

	for (int k = lbin; k <= hbin; k++) {
		if (k == lbin) {
			Npi0 += h->GetBinContent(k)*lratio;
			Npi0_err += pow(h->GetBinError(k)*lratio, 2);
		}
		else if (k == hbin) {
			Npi0 += h->GetBinContent(k)*rratio;
			Npi0_err += pow(h->GetBinError(k)*rratio, 2);
		}
		else {
			Npi0 += h->GetBinContent(k);
			Npi0_err += pow(h->GetBinError(k), 2);
		}
	}
    float Nbest2 = 0, Nside = 0, Nomg = 0, Nbkg = 0;
    if (h->GetListOfFunctions()->FindObject("best2")) {
        TF1* best2((TF1*)h->GetListOfFunctions()->FindObject("best2"));
        Nbest2 = best2->Integral(x1, x2) / bw;
    }
    if (h->GetListOfFunctions()->FindObject("accfcn")) {
        TF1* side((TF1*)h->GetListOfFunctions()->FindObject("accfcn"));
        Nside = side->Integral(x1, x2) / bw;
    }
    if (h->GetListOfFunctions()->FindObject("omgfcn")) {
        TF1* omg((TF1*)h->GetListOfFunctions()->FindObject("omgfcn"));
        Nomg = omg->Integral(x1, x2) / bw;
    }
    if (h->GetListOfFunctions()->FindObject("bkg")) {
        TF1* bkg((TF1*)h->GetListOfFunctions()->FindObject("bkg"));
        Nbkg = bkg->Integral(x1, x2) / bw;
    }
    Npi0 -= Nbest2;
    Npi0 -= Nside;
    Npi0 -= Nomg;
    Npi0 -= Nbkg;
	Npi0_err = sqrt(Npi0 + 2 * (Nbest2 + Nside + Nomg + Nbkg));

    TF1* pkfcn = (TF1*)h->GetListOfFunctions()->FindObject("pkfcn");
    //ratio = pkfcn->Integral(x1, x2) / pkfcn->Integral(uimanager.fit_low_limit(), uimanager.fit_high_limit());
    Npi0 /= ratio;
    Npi0_err /= ratio;
}

void fitrod::setpar(TF1* fitfcn, fitmethod *fit, vector<double>& par, bool fixed, bool par_initiated, double bw, double angle) {
	if (uimanager.get_method() == 1) {
		fitfcn->FixParameter(1, par[1]);
		fitfcn->FixParameter(2, 0);
		fitfcn->FixParameter(3, 0);
		fitfcn->FixParameter(4, 0);
		fitfcn->FixParameter(5, 0);
	}
	else if (!fixed && !par_initiated) {
		fitfcn->SetParameter(1, par[1]);
		fitfcn->SetParameter(2, par[2]);
		fitfcn->SetParameter(3, par[3]);
		fitfcn->SetParameter(4, par[4]);
		fitfcn->SetParameter(5, par[5]);
	
		fitfcn->SetParLimits(2, bw / sqrt(3), par[2] * 1.5);
		fitfcn->SetParLimits(3, 0., 1.);
		fitfcn->SetParLimits(4, -0.01, 0.01);
		fitfcn->SetParLimits(5, bw / sqrt(3), par[5] * 1.5);
	}
	else if (!fixed && par_initiated) {
		fitfcn->SetParameter(1, par[1]);
		fitfcn->SetParameter(2, fabs(par[2]));
		fitfcn->SetParameter(3, par[3]);
		fitfcn->SetParameter(4, par[4]);
		fitfcn->SetParameter(5, fabs(par[5]));
	
		fitfcn->SetParLimits(1, par[1] - fabs(par[2]), par[1] + fabs(par[2]));
		fitfcn->SetParLimits(2, fabs(par[2])*0.5, fabs(par[2])*1.5);
		fitfcn->SetParLimits(3, 0.5*par[3], 1.);
		fitfcn->SetParLimits(4, TMath::Min(par[4] * 0.2, par[4] * 2), TMath::Max(par[4] * 0.2, par[4] * 2));
		fitfcn->SetParLimits(5, fabs(par[5])*0.5, fabs(par[5])*1.5);
	}
	else if (fixed) {
		fitfcn->SetParameter(1, par[1]);
		fitfcn->SetParLimits(1, par[1] - fabs(par[2]), par[1] + fabs(par[2]));
		fitfcn->FixParameter(2, par[2]);
		fitfcn->FixParameter(3, par[3]);
		fitfcn->FixParameter(4, par[4]);
		fitfcn->FixParameter(5, par[5]);
	}
    //for acc. omega and polynomial background
    if (uimanager.btdiff_correction() == 1) {
        fit->init_bkg(&par[7]);
        fitfcn->FixParameter(6, par[6]);//acc.
        fitfcn->SetParameter(7, par[7]);
        fitfcn->SetParameter(8, par[8]);
        fitfcn->SetParameter(9, par[9]);
        fitfcn->SetParameter(10, par[10]);
        //fitfcn->SetParameter(11, par[11]);
        //fitfcn->SetParameter(12, par[12]);
        fitfcn->FixParameter(11, 0);
        fitfcn->FixParameter(12, 0);
        fitfcn->FixParameter(13, par[13]);//omega.

        int lbin1 = h->FindBin(uimanager.get_hist_limits()[0]);
        int lbin2 = h->FindBin(-0.03);
        int hbin1 = h->FindBin(0.03);
        int hbin2 = h->FindBin(uimanager.get_hist_limits()[1]);
        int count_bins = 0;
        double height = 0.;
        for (int i = lbin1; i <= hbin2; i++) {
            if (lbin2 <= i || i <= hbin1) continue;
            count_bins++;
            height += h->GetBinContent(i);
        }
        height /= count_bins;
        height *= 5;
        fitfcn->SetParLimits(7, 0., height);
        fitfcn->SetParLimits(8, -height / 10, height / 10);
        fitfcn->SetParLimits(9, -height / 100, height / 100);
        fitfcn->SetParLimits(10, -height / 1000, height / 1000);
        //fitfcn->SetParLimits(11, -height / 1000, height / 1000);
        //fitfcn->SetParLimits(12, -height / 1000, height / 1000);
    }

	//include best tdiff correction (subtract 2nd best tdiff)
	if (uimanager.btdiff_correction() == 0) {//don't use the correction
		fitfcn->FixParameter(14, 0.);
		for (int i = 15; i < 20; i++)
			fitfcn->FixParameter(i, 0.5);
	}
	else if (uimanager.btdiff_correction() == 1) {// use the correction
		for (int i = 14; i < 20; i++)
			fitfcn->FixParameter(i, par[i]);
	}
	else {// calculate the correction
        for (int i = 6; i < 14; i++) fitfcn->FixParameter(i, 0.);
		for (int i = 14; i < 20; i++)
			fitfcn->SetParameter(i, par[i]);
        fitfcn->SetParLimits(15, -0.1, 0.1);
        fitfcn->SetParLimits(16, 0.05, 0.1);
        fitfcn->SetParLimits(17, 0, 1);
        fitfcn->SetParLimits(18, -0.15, 0.15);
        fitfcn->SetParLimits(19, 0.05, 0.2);
	}

    if (uimanager.btdiff_correction()) {
        fitfcn->SetParName(14, "N_{2nd best}");
        fitfcn->SetParName(15, "mean1");
        fitfcn->SetParName(16, "#sigma_{1}");
        fitfcn->SetParName(17, "fraction");
        fitfcn->SetParName(18, "mean2");
        fitfcn->SetParName(19, "#sigma_{2}");
    }

    fitfcn -> SetParameter(0, par[0]);
    fitfcn -> SetParLimits(0,100,10*par[0]);
};

double fitrod::find_mean(vector<double> par, double angle)
{
	const int count = 50;
	vector<vector<double> > archi2(2, vector<double>(count+1, 0.));
	vector<vector<double> > findpar(count+1, vector<double>(_npar, 0.));
	double scanrange = 0.008, mid = 0, x2 = 100., y2;

	for (int i = 0; i <= count; i++)
	{
		par[1] = mid + (i - count / 2.) / count*scanrange;
		TFitResultPtr r = fitting(par, angle, "N0Q");
		double chi2 = r->Chi2();
		int NDF = r->Ndf();
		for (int j = 0; j<_npar; j++) findpar[i][j] = r->Parameter(j);
		archi2[0][i] = findpar[i][1];
		archi2[1][i] = chi2 / NDF;
		if (fabs(par[1]) < fabs(x2)) {
			x2 = par[1];
			y2 = archi2[1][i];
		}
	}

	TGraph mean_chi2(count+1, &archi2[0][0], &archi2[1][0]);
	mean_chi2.SetMarkerStyle(20);
	mean_chi2.SetTitle(Form("#theta: %.2f", angle));
	mean_chi2.GetXaxis()->SetTitle("peak position");
	mean_chi2.GetYaxis()->SetTitle("#chi^{2}/NDF");

	double x1 = archi2[0][0];
	double x3 = archi2[0][count];
	double y1 = archi2[1][0];
	double y3 = archi2[1][count];
	double b = -(x1*x1*(y3 - y2) - x3*x3*(y1 - y2)) / 2 / (x3*(y1 - y2) - x1*(y3 - y2));
	double a = (y1 - y2) / (x1*x1 - 2 * b*x1);
	double c = y1 - a*(x1 - b)*(x1 - b);

	double fmin = max(b - 0.002, mid - scanrange / 2);
	double fmax = min(b + 0.002, mid + scanrange / 2);

	TF1 fitfcn("fitfcn", "[0]*(x-[1])**2+[2]", fmin, fmax);

	fitfcn.SetParameters(a, b, c);
	mean_chi2.Fit("fitfcn", "QR");
    mean_chi2s.push_back(mean_chi2);
    
    double ret = fitfcn.GetParameter(1);
    if (isnan(ret)) ret = 0;
    return ret; 
}

TFitResultPtr fitrod::fitting(vector<double>& par, double angle, TString option){
	TString op("SR");
    bool add_fit_func = false;
    if (string(option).find("0") == string::npos && string(option).find("N") == string::npos) add_fit_func = true;
    op = op + option;

	Double_t bw = h->GetXaxis()->GetBinWidth(1);
	double ratio_mc = 1, Nmcpk = 0;
	if (uimanager.get_method() != 2) {
		ratio_mc = h->GetBinWidth(1) / mcpk->GetBinWidth(1);
		int nbins = mcpk->GetNbinsX();
		for (int i = 1; i <= nbins; i++) Nmcpk += mcpk->GetBinContent(i);
	}
	else Nmcpk = h->GetEntries();
	double ratio_acc = h->GetBinWidth(1)/accidentals->GetBinWidth(1);
	double ratio_omg = h->GetBinWidth(1)/omega->GetBinWidth(1);

	fitmethod *fit = new fitmethod(accidentals, mcpk, omega, h, uimanager, accidentals->GetEntries() / ratio_acc, Nmcpk / ratio_mc,
		omega->GetEntries() / ratio_omg, angle);
    
    TF1 *fitfcn = new TF1("fitmc", fit, &fitmethod::fitfcn, uimanager.fit_low_limit(), uimanager.fit_high_limit(), _npar - 2);
	fitfcn -> SetLineColor(kRed);
	fitfcn -> SetLineWidth(1);

    if (uimanager.get_method() != 2) {
        if (uimanager.btdiff_correction() != 2) fitfcn -> SetParNames("N_{#pi^{0}}","offset","N_{acc.}");
        else fitfcn -> SetParNames("N_{substituted}","offset","N_{acc.}");
    } else {
        if (uimanager.btdiff_correction() != 2) fitfcn -> SetParNames("N_{#pi^{0}}","mean1","#sigma_{1}","fraction","distance","#sigma_{2}");
        else fitfcn -> SetParNames("N_{substituted}","mean1","#sigma_{1}","fraction","distance","#sigma_{2}");
    }
	fitfcn -> SetParName(6,"N_{acc.}");
	fitfcn -> SetParName(7,"p0");
	fitfcn -> SetParName(8,"p1");
	fitfcn -> SetParName(9,"p2");
	fitfcn -> SetParName(10,"p3");
	fitfcn -> SetParName(11,"p4");
	fitfcn -> SetParName(12,"p5");
	fitfcn -> SetParName(13,"N_{#omega}");
    
    bool fixed = false, par_initiated = false;
	TFitResultPtr r;
    if(init_par.size()) {
        double chi2 = 0;
        fixed = true; 
        setpar(fitfcn, fit, par, fixed, par_initiated, bw, angle);
        r = h -> Fit(fitfcn, op+"LRMBEQ");
        chi2 = r->Chi2()/r->Ndf(); 
        for(int k = 0; k < _npar - 2; k++){
            par[k] = r->Parameter(k);
        }
        if(uimanager.get_method() == 2) {
            fixed = false;
            par_initiated = true;
			setpar(fitfcn, fit, par, fixed, par_initiated, bw, angle);
            r = h -> Fit(fitfcn, op+"LRMBEQ");//2nd fit
            if (r->Chi2()/r->Ndf() > chi2) {
                fixed = true;
				setpar(fitfcn, fit, par, fixed, par_initiated, bw, angle);
                r = h -> Fit(fitfcn, op+"LRMBEQ");
            }
        }
        for(int k = 0;k < _npar - 2; k++){
            par[k] = r->Parameter(k);
        }
    } else {
        if (user_defined) fixed = true;
		setpar(fitfcn, fit, par, fixed, par_initiated, bw, angle);
        r = h -> Fit(fitfcn, op+"LRMBE");//1st fit
        for(int k = 0; k < _npar - 2; k++){
            par[k] = r->Parameter(k);
        }
    }
///////////////////////////////////////////////////////////////////////////////////////
    if (add_fit_func) {
        TF1 *pkfcn = NULL, *accfcn = NULL;
        if (uimanager.get_method() != 2) {
            pkfcn = new TF1("pkfcn", fit, &fitmethod::mcfit, uimanager.fit_low_limit(), uimanager.fit_high_limit(), 2);
            pkfcn->SetParameters(par[0],par[1]);
        } else {
            pkfcn = new TF1("pkfcn", fit, &fitmethod::DoubleGaussianFit, uimanager.fit_low_limit(), uimanager.fit_high_limit(), 6);
            pkfcn->SetParameters(par[0],par[1],par[2],par[3],par[4],par[5]);
        }
        if (fabs(par[6]) > 1.e-3) {
            accfcn = new TF1("accfcn", fit, &fitmethod::accfit, uimanager.fit_low_limit(), uimanager.fit_high_limit(), 1);
            accfcn->SetParameter(0,fitfcn->GetParameter(6));
            accfcn->SetLineColor(38);
            h -> GetListOfFunctions() -> Add(accfcn);
        }
        pkfcn->SetLineColor(kGreen);

        TF1 *bkg = new TF1("bkg", fit, &fitmethod::CubBkg, uimanager.fit_low_limit(), uimanager.fit_high_limit(), 6);
        bkg->SetParameters(par[7], par[8], par[9], par[10], par[11], par[12]);
        bkg->SetLineColor(kBlue);
        
        TF1 *omgfcn = new TF1("omgfcn", fit, &fitmethod::omegafit, uimanager.fit_low_limit(), uimanager.fit_high_limit(), 1);
        omgfcn->SetLineColor(kBlack);
        omgfcn->SetParameter(0,par[13]);
        
        h -> SetLineColor(kBlack);
        h -> GetListOfFunctions() -> Add(pkfcn);

        if (uimanager.use_poly_bkg()) h -> GetListOfFunctions() -> Add(bkg);
        if (uimanager.sub_omg()) h -> GetListOfFunctions() -> Add(omgfcn);

        if (uimanager.btdiff_correction() == 1) {
            TF1 *best2 = new TF1("best2", fit, &fitmethod::besttdiff, uimanager.fit_low_limit(), uimanager.fit_high_limit(), 6);
            best2->SetParameters(par[0]*par[14], par[15], par[16], par[17], par[18], par[19]);
            best2->SetLineColor(kPink);

            h -> GetListOfFunctions() -> Add(best2);
        }
        
        par[20] = bkg->Integral(uimanager.fit_low_limit(), uimanager.fit_high_limit()) / bw;
        par[21] = omgfcn->Integral(uimanager.fit_low_limit(), uimanager.fit_high_limit()) / bw;

		if (hbkg) {
			for (int i = 0; i < mdiv; ++i) {
				float x = hbkg->GetBinCenter(i);
				if (x < uimanager.fit_low_limit() + hbkg->GetBinWidth(i) || x > uimanager.fit_high_limit() - hbkg->GetBinWidth(i))
					hbkg->SetBinContent(i, h->GetBinContent(i));
				else {
					float y = h->GetBinContent(i) - pkfcn->Eval(x) - omgfcn->Eval(x);
					if (uimanager.btdiff_correction() == 1) y -= best2->Eval(x);
					hbkg->SetBinContent(i, y);
					hbkg->SetBinError(i, sqrt(fabs(y));
				}
			}
		}

        if (uimanager.get_method() == 2) {
            TF1 *g1 = new TF1("g1", Form("[0]*%f*(1-[3])/[2]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])", bw / sqrt(2 * 3.14159265359)), uimanager.get_hist_limits()[0], uimanager.get_hist_limits()[1]);
            g1 -> SetParameters(par[0],par[1],par[2],par[3]);
            g1 -> SetLineColor(kBlue);
            g1 -> SetLineStyle(3);
            TF1 *g2 = new TF1("g2", Form("[0]*%f*[3]/[2]*exp(-0.5*(x-[1]-[4])*(x-[1]-[4])/[2]/[2])", bw / sqrt(2 * 3.14159265359)), uimanager.get_hist_limits()[0], uimanager.get_hist_limits()[1]);
            g2 -> SetParameters(par[0],par[1],par[5],par[3],par[4]);
            g2 -> SetLineColor(kBlue);
            g2 -> SetLineStyle(3);
            
            h -> GetListOfFunctions() -> Add(g1);
            h -> GetListOfFunctions() -> Add(g2);

            if (uimanager.btdiff_correction() == 2) {
                TF1 *g3 = new TF1("g3", Form("[0]*%f*(1-[3])/[2]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])", bw / sqrt(2 * 3.14159265359)), uimanager.get_hist_limits()[0], uimanager.get_hist_limits()[1]);
                g3 -> SetParameters(par[14],par[15],par[16],par[17]);
                g3 -> SetLineColor(kPink);
                g3 -> SetLineStyle(3);
                TF1 *g4 = new TF1("g4", Form("[0]*%f*[3]/[2]*exp(-0.5*(x-[1]-[4])*(x-[1]-[4])/[2]/[2])", bw / sqrt(2 * 3.14159265359)), uimanager.get_hist_limits()[0], uimanager.get_hist_limits()[1]);
                g4 -> SetParameters(par[14],par[15],par[19],par[17],par[18]);
                g4 -> SetLineColor(kPink);
                g4 -> SetLineStyle(3);

                h -> GetListOfFunctions() -> Add(g3);
                h -> GetListOfFunctions() -> Add(g4);
            }
        }
    }

    if (!add_fit_func) {
        delete fit;
        delete fitfcn;
    }

	return r;
}

void fitrod::make_plot() {
	TString pdfname("altfit_mean" + uimanager.output_filename("fit") + ".pdf");

	TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
	c1->Divide(2, 2);
	TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1200);
	c2->Divide(2, 2);
	c1->SaveAs(pdfname + "[");

	for (int i = 0; i < mean_chi2s.size(); ++i){
		if (i / 4 >= mean_chi2s.size() / 4)c2->cd(i % 4 + 1);
		else c1->cd(i % 4 + 1);
		mean_chi2s[i].SetMinimum(0);
		mean_chi2s[i].Draw("ap");

		if (i == mean_chi2s.size() - 1)c2->SaveAs(pdfname);
		else if (i % 4 == 3)c1->SaveAs(pdfname);
	}
	c1->SaveAs(pdfname + "]");
}
