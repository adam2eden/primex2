#ifndef primex
#define primex

#include <vector>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <map>
#include <unordered_map>

#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
#include "TRandom2.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTree.h"
#include "TChain.h"
#include "TPaveText.h"
#include "TCutG.h"
#include "TAxis.h"
#include "TVectorD.h"

using namespace std;

const double max_angle = 2.5;
const int nangle = 125;

const int fit_nangle = 125, nmc = 640, nrec = 160;
const int nmc1 = 540, nrec1 = 135;
const double min_tdiff = -14, max_tdiff = 25, dndt_fit_range = 2.5;
const int nregion = 2, mdiv = 100;

class UImanager {
	int _nsigma, _method, _mc, _num_runs, _outer, _btc, _use_poly_bkg, _sub_omg;

	double _sigma_start, _sigma_step, _lelas, _tsi, _ta_corr, _adcerr, _bit_corr,  _br_corr, _tgt_lumi, _flux, _flux_omega, _xmin, _xmax, _tdiff_cut;
	bool _subacc, _best_tdiff, _use_Npi0;
    string prog;

    map<pair<string, string>, bool> _input_option_list;

	vector<int> _run;
	vector<double> _input_angles, _output_angles, _btcorr_angles, _limits, _elas, _invm_cut, _z;
    vector<vector<double> > _btcorr_par;
	vector<map<int, double> > _flux_map;
    
    bool input_option(char* arg, string opt);
    void help();
public:
	UImanager(int argc, char* argv[], string prog);
	
	int target();
	int input_nbins() { return int(_input_angles.size() - 1); }
	int output_nbins() { return int(_output_angles.size() - 1); }
    int btcorr_nbins() { return int(_btcorr_angles.size() - 1); }

	const vector<double>& get_input_angles() { return _input_angles; }
	const vector<double>& get_output_angles() { return _output_angles; }
    const vector<double>& get_btcorr_angles() { return _btcorr_angles; }
	const vector<int>& get_runs() { return _run; }
	const vector<double>& get_hist_limits() { return _limits; }
    const vector<vector<double> >& get_btcorr_pars() { return _btcorr_par; }

	int get_num_runs() { return _num_runs; }
	int get_method() { return _method; }
	int get_nsigma() { return _nsigma; }
	int btdiff_correction() { return _btc; }
    int use_poly_bkg() { return _use_poly_bkg; }
    int sub_omg() { return _sub_omg; };

    void set_btc(int btc) { _btc = btc; }
    void set_poly_bkg(int bkg) { _use_poly_bkg = bkg; }
    void set_sub_omg(int omg) { _sub_omg = omg; }
    void set_fit_range(double xmin, double xmax) { _xmin = xmin; _xmax = xmax; }
	void set_wo(int wo) { _outer = wo; }

	double get_adcerr() { return _adcerr; }
	double get_bit_corr() { return _bit_corr; }
	double get_br_corr() { return _br_corr; }
	double get_lelas_cut() { return _lelas; }
	double elas_low_cut() { return _elas[0]; }
	double elas_high_cut() { return _elas[1]; }
	double invm_low_cut() { return _invm_cut[0]; }
	double invm_high_cut() { return _invm_cut[1]; }
	double get_tgt_lumi() { return _tgt_lumi; }
	double get_sigma_start() { return _sigma_start; }
	double get_sigma_step() { return _sigma_step; }
    double get_tdiff_cut() { return _tdiff_cut; }
    double get_tsi() { return _tsi; }
    double get_ta_corr() { return _ta_corr; }
    double getz() { return _z[target()]; }

	double fit_low_limit() { return _xmin; }
	double fit_high_limit() { return _xmax; }

	int fitting_method() { return _method; }
	bool subtract_acc() { return _subacc; }
	bool ismc() { return _mc; }
	bool isouter() { return _outer; }
    bool best_tdiff() { return _best_tdiff; }
    bool use_Npi0() { return _use_Npi0; }

	double flux() { return _flux; }
	double flux(int _run_number) { return _flux_map[target()][_run_number]; }
	double flux_omega() { return _flux_omega; }
	
	string input_filename(const string filetype);
	string output_filename(const string filetype);
};


bool isequal(float x, float y);
bool isposint(const string&);
int open_err (const string fname);
int open_err(const TString fname);
int input_err ();
int wrong_folder ();

#endif
