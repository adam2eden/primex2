#ifndef alt_model_fit
#define alt_model_fit

#include "TF1.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "primex.h"

using namespace std;
class fitmethod {
private:
	TH1F *accidentals, *mcpk, *omega, *hdata;
	UImanager uimanager;
	double acc_norm, mc_norm, omg_norm, bw, anchor;
	int large_angle;
	
public:

	fitmethod(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, UImanager& uimanager, double acc_norm, double mc_norm, double omg_norm, double angles) :
		accidentals{ h1 }, mcpk{ h2 }, omega{ h3 }, hdata{ h4 }, uimanager(uimanager), acc_norm{ acc_norm }, mc_norm{ mc_norm }, omg_norm{ omg_norm }, bw{ h4->GetBinWidth(1) }
	//method: 0 mc, acc. 1 mc no acc subtraction 2. double gaussian no acc subtraction
	{
		if (uimanager.get_method() == 1) {
			if (angles >= 0.6) {
                anchor = -0.03;
				large_angle = 1;
                //anchor = 0;
				//large_angle = 0;
			}
			else {
				anchor = 0;
				large_angle = 0;
			}
		}
		else {
			if (angles >= 3) {
				anchor = -0.04;
				large_angle = 1;
			}
			else {
				anchor = 0;
				large_angle = 0;
			}
		}
	};

	Double_t accfit(Double_t *x, Double_t *para){
		Double_t t = x[0];
		int bin = accidentals->FindBin(t);
		Double_t acc = para[0] * (accidentals->GetBinContent(bin)) / acc_norm;
		return acc;
	}

	Double_t mcfit(Double_t *x, Double_t *para){
		Double_t t = x[0] - para[1];
		int bin = mcpk->FindBin(t);
		Double_t mc = para[0] * (mcpk->GetBinContent(bin)) / mc_norm;
		return mc;
	}

	Double_t DoubleGaussianFit(Double_t *x, Double_t *para) {
		double t1 = (x[0] - para[1]) / para[2];
		double t2 = (x[0] - para[1] - para[4]) / para[5];
		return para[0] / sqrt(2 * 3.14159265359)*bw*((1 - para[3]) / fabs(para[2])*exp(-0.5*t1*t1) + para[3] / fabs(para[5])*exp(-0.5*t2*t2));
	}

	Double_t besttdiff(Double_t *x, Double_t *para) {
		double t1 = (x[0] - para[1]) / para[2];
		double t2 = (x[0] - para[1] - para[4]) / para[5];
		return para[0] / sqrt(2 * 3.14159265359)*bw*((1 - para[3]) / fabs(para[2])*exp(-0.5*t1*t1) + para[3] / fabs(para[5])*exp(-0.5*t2*t2));
	}

	Double_t omegafit(Double_t *x, Double_t *para){
		Double_t t = x[0];
		int bin = omega->FindBin(t);
		Double_t omg = para[0] * (omega->GetBinContent(bin)) / omg_norm;
		return omg;
	}

	Double_t CubBkg(Double_t *x, Double_t *para){
		Double_t t = (x[0] - anchor) / bw;
		Double_t bkg = 0;
        if (!large_angle)bkg = para[0] + para[1] * t + para[2] * t*t + para[3] * t*t*t;
        //if (!large_angle)bkg = para[0] + para[1] * t + para[2] * t*t;
        //if (!large_angle)bkg = para[0] + para[1] * t + para[2] * t*t + para[3] * t*t*t + para[4] * t * t * t * t;
        //if (!large_angle)bkg = para[0] + para[1] * t + para[2] * t*t + para[3] * t*t*t + para[4] * t * t * t * t + para[5]*t*t*t*t*t;
		else if (t>0)bkg = para[0] + para[1] * t + para[2] * t*t + para[3] * t*t*t;
		else if (uimanager.target() == 0 && !uimanager.subtract_acc()) bkg = para[0] + para[1] * t + para[2] * t*t;
		else if (uimanager.target() == 1 || uimanager.subtract_acc()) bkg = para[0] + para[1] * t;
		return bkg;
	}

	void init_bkg(double *para) {
		double count_xmin = 0, count_xmax = 0, count_anchor = 0, count = 0;
		double t_xmin = (uimanager.fit_low_limit() - anchor) / bw, t_xmax = (uimanager.fit_high_limit() - anchor) / bw, t = (0.075 - anchor) / bw;
		int bin1 = hdata->FindBin(uimanager.fit_low_limit());
		int bin2 = hdata->FindBin(uimanager.fit_high_limit());
		int bin3 = hdata->FindBin(anchor);
		int bin4 = hdata->FindBin(0.075);
		for (int i = -5; i <= 5; i++) {
			count_xmin += hdata->GetBinContent(bin1 + i);
			count_xmax += hdata->GetBinContent(bin2 + i);
			if (anchor <= -0.03 || anchor >= 0.03)count_anchor += hdata->GetBinContent(bin3 + i);
			else count_anchor += (hdata->GetBinContent(hdata->FindBin(-0.03) + i) + hdata->GetBinContent(hdata->FindBin(0.03) + i)) / 2;
			count += hdata->GetBinContent(bin4 + i);
		}
		count_xmin /= 11;
		count_xmax /= 11;
		count_anchor /= 11;
		count /= 11;

		para[0] = count_anchor;

		double term1, term2;
		term1 = (t_xmax*count_xmin - t_xmin*count_xmax + para[0] * (t_xmin - t_xmax)) / t_xmin / t_xmax / (t_xmin - t_xmax);
		term2 = (t_xmax*count - t*count_xmax + para[0] * (t - t_xmax)) / t / t_xmax / (t - t_xmax);
		para[3] = (term1 - term2) / (t_xmin - t);
		para[2] = term1 - para[3] * (t_xmin + t_xmax);
		para[1] = (count_xmin - para[0] - para[2] * t_xmin*t_xmin - para[3] * t_xmin*t_xmin*t_xmin) / t_xmin;
	}

	Double_t fitfcn(Double_t *x, Double_t *para){
		//total 20 parameters
        Double_t ret = 0;
		if (uimanager.get_method() == 1)
            ret = mcfit(x, para) + accfit(x, &para[6]);
        else
            ret = DoubleGaussianFit(x, para) + accfit(x, &para[6]);
        if (uimanager.use_poly_bkg()) ret += CubBkg(x, &para[7]);

        if (uimanager.sub_omg()) ret += omegafit(x, &para[13]);

        if (uimanager.btdiff_correction() == 1) ret += para[0]*besttdiff(x, &para[14]);
        else if (uimanager.btdiff_correction() == 2) ret += besttdiff(x, &para[14]);

        return ret;
	}
};

class fitrod {
	int _npar;
    bool user_defined;
	
	UImanager uimanager;
	vector<vector<double> > init_par;
    vector<TGraph> mean_chi2s;

	TH1F *h, *accidentals, *mcpk, *omega, *hbkg;
private:
	double find_mean(vector<double> par, double angle);
public:
	fitrod(const UImanager& uimanager, vector<vector<double> >& init_par);

    void set_uimanager(const UImanager& _uimanager) { uimanager = _uimanager; }
    void define_par(bool input);
    void sethist(TH1F* _h, TH1F* _accidentals, TH1F* _mcpk, TH1F *_omega, TH1F *_hbkg);
	void setpar(TF1* fitfcn, fitmethod* fit, vector<double>& par, bool fixed, bool par_initiated, double bw, double angle);
    void setinitpar(vector<vector<double> >& input_par);
	TFitResultPtr fitting(vector<double>& par, double angle, TString option);
	void initialize(vector<double>& par, double angle);
	void calc_Npi0(vector<double>& par, double& Npi0, double& Npi0_err);

	int npar() { return _npar; }

    void make_plot();
};

#endif
