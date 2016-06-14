#include "primex.h"

vector<vector<double> > angles{ { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.9, 2.1, 2.3, 2.5 },
{ 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.2, 1.5, 1.8, 2.1, 2.5 } };

vector<vector<double> > btcorr_angles{ {0., 0.5, 1.0, 1.5, 2.5}, {0., 0.5, 1.0, 1.5, 2.5} };

const vector<vector<double>> flux_target = {{ 5.3034887516e12, 2.4482987753e12 }, { 5.3e12, 2.4e12 }};
const vector<double> flux_omega_target = { 410.1e12, 263.0e12 };
const vector<double> tsi = { 0.04973456e-6, 0.17700914e-6 };
const vector<vector<double>> ta_corr = {{ 0.960488*1.00187, 0.970*1.0365 }, {1, 1}};
const vector<double> bit_corr = { 0.9885, 0.98569 };
const vector<double> tgt_lumi = { 0.04973456e-3, 0.17700914e-3 };

UImanager::UImanager(int argc, char* argv[], string prog)
:_nsigma(31), _method(1), _mc(0), _num_runs(0), _outer(1), _btc(1), _use_poly_bkg(1), _sub_omg(1),
_sigma_start(0.85), _sigma_step(0.01), _lelas(1.15), _tsi(tsi[target()]), _adcerr(0.995), _bit_corr(bit_corr[target()]), _br_corr(1. - 1.2e-2), _tgt_lumi(tgt_lumi[target()]), _flux(0), _flux_omega(0), _xmin(-0.1), _xmax(0.1), _tdiff_cut(6.),
_subacc(true), _best_tdiff(true), _use_Npi0(false), prog(prog)
{
	_limits = { -0.2, 0.2 };
	_elas = { -100, 100 };
	_invm_cut = { 0, 1 };
    _z = {702, 702.18};

    _btcorr_angles = btcorr_angles[target()];
    
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-b"    ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-out"  ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-elas" ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-invm" ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass",              "-lelas"), true) );

    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-b"    ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-m"    ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-in"   ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-elas" ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-invm" ), true) );
    _input_option_list.insert( make_pair( make_pair("btcorr",                   "-lelas"), true) );

    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-r"    ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-mc"   ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-out"  ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-elas" ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-invm" ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_mc",           "-lelas"), true) );

    _input_option_list.insert( make_pair( make_pair("rotatedmass_omg",          "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_omg",          "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_omg",          "-out"  ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_omg",          "-elas" ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_omg",          "-invm" ), true) );
    _input_option_list.insert( make_pair( make_pair("rotatedmass_omg",          "-lelas"), true) );

    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc",         "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc",         "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc",         "-out"  ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc",         "-elas" ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc",         "-lelas"), true) );

    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-b"    ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-m"    ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-l"    ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-in"   ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-acc"  ), true) );
	_input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-btc"  ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-out"  ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-elas" ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-invm" ), true) );
    _input_option_list.insert( make_pair( make_pair("fast_fitmodel_mc_correct", "-lelas"), true) );

    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-b"    ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-m"    ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-l"    ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-acc"  ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-btc"  ), true) );
	_input_option_list.insert( make_pair( make_pair("alt_model_corrret",        "-btc"  ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-out"  ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-elas" ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-invm" ), true) );
    _input_option_list.insert( make_pair( make_pair("alt_model_correct",        "-lelas"), true) );

    _input_option_list.insert( make_pair( make_pair("altfit",                   "-b"    ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-m"    ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-l"    ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-r"    ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-in"   ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-acc"  ), true) );
	_input_option_list.insert( make_pair( make_pair("altfit",                   "-btc"  ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-out"  ), true) );
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-elas" ), true) );       
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-invm" ), true) );       
    _input_option_list.insert( make_pair( make_pair("altfit",                   "-lelas"), true) );       

    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-b"    ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-m"    ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-l"    ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-r"    ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-in"   ), true) );
	_input_option_list.insert( make_pair( make_pair("plotalt",                  "-mc"   ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-acc"  ), true) );
	_input_option_list.insert( make_pair( make_pair("plotalt",                  "-btc"  ), true) );
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-elas" ), true) );       
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-invm" ), true) );       
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-npi0" ), true) );       
    _input_option_list.insert( make_pair( make_pair("plotalt",                  "-lelas"), true) );       

    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-b"    ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-f"    ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-m"    ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-l"    ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-r"    ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-t"    ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-in"   ), true) );
	_input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-mc"   ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-wo"   ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-acc"  ), true) );
	_input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-btc"  ), true) );
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-elas" ), true) ); 
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-invm" ), true) ); 
    _input_option_list.insert( make_pair( make_pair("dndt_fit",                 "-lelas"), true) ); 

	map<int, double> flux_map_si = { { 64716, 26.98007683 }, { 64717, 35.27877433 }, { 64718, 30.67485698 }, { 64719, 46.34924759 }, { 64720, 40.79595001 },
	{ 64721, 25.06975045 }, { 64722, 25.99939587 }, { 64723, 41.77658505 }, { 64724, 46.10743272 }, { 64725, 25.10158572 }, { 64726, 51.31902668 }, { 64727, 35.79349243 },
	{ 64728, 52.04475919 }, { 64729, 47.99889222 }, { 64730, 34.75644096 }, { 64731, 37.70315338 }, { 64732, 32.02591986 }, { 64733, 35.14052011 }, { 64734, 30.41480142 },
	{ 64735, 41.04003273 }, { 64736, 38.94370050 }, { 64737, 33.75222240 }, { 64738, 49.47490513 }, { 64739, 43.79289900 }, { 64740, 48.03947775 }, { 64741, 44.01384653 },
	{ 64742, 52.23303475 }, { 64743, 70.45428551 }, { 64751, 23.29301464 }, { 64756, 16.62722158 }, { 64757, 44.98019248 }, { 64758, 20.56092133 }, { 64759, 40.08394689 },
	{ 64760, 36.94683140 }, { 64761, 37.90874270 }, { 64762, 37.04898496 }, { 64763, 36.21053457 }, { 64764, 20.10613809 }, { 64767, 35.08118866 }, { 64768, 25.08062774 },
	{ 64769, 35.57638287 }, { 64770, 36.25859915 }, { 64771, 31.77146302 }, { 64772, 37.63284622 }, { 64773, 35.57282689 }, { 64779, 37.06284621 }, { 64780, 34.09496591 },
	{ 64782, 14.15583265 }, { 64784, 35.97715582 }, { 64785, 32.14214913 }, { 64786, 21.75931883 }, { 64802, 39.11343255 }, { 64803, 43.84187137 }, { 64804, 41.74657143 },
	{ 64808, 16.69669220 }, { 64809, 13.25534224 }, { 64839, 44.00606058 }, { 64841, 41.83109022 }, { 64850, 2.16885883 }, { 64852, 17.36536171 }, { 64854, 32.52440714 },
	{ 64856, 40.19723941 }, { 64857, 41.25243045 }, { 64858, 40.63380840 }, { 64859, 41.44339554 }, { 64860, 40.93464042 }, { 64861, 36.71787561 }, { 64863, 25.42014391 },
	{ 64864, 43.33557556 }, { 64865, 52.07718034 }, { 64866, 43.51278096 }, { 64867, 47.85746071 }, { 64868, 44.89913571 }, { 64869, 44.65176967 }, { 64870, 43.48247663 },
	{ 64871, 39.94455166 }, { 64872, 47.50804065 }, { 64873, 2.37141206 }, { 64874, 41.79692736 }, { 64875, 36.85035360 }, { 64887, 42.92478740 }, { 64888, 1.46275860 },
	{ 64889, 29.30952191 }, { 64890, 44.99314081 }, { 64891, 37.12997275 }, { 64892, 41.89844311 }, { 64893, 42.14793882 }, { 64894, 45.09471900 }, { 64896, 44.68275330 },
	{ 64897, 43.72364997 }, { 64907, 44.99857641 }, { 64908, 46.23832514 }, { 64909, 40.43249687 }, { 64911, 45.09889611 }, { 64912, 45.09158136 }, { 64913, 17.65935774 },
	{ 64914, 11.90173005 }, { 64915, 45.10728224 }, { 64917, 45.62759055 }, { 64918, 45.21398451 }, { 64919, 46.35152515 }, { 64920, 44.43399858 }, { 64921, 46.16176888 },
	{ 64923, 39.37130658 }, { 64924, 8.53994316 }, { 64925, 41.76490799 }, { 64926, 43.30522181 }, { 64927, 41.40531425 }, { 64928, 42.96792858 }, { 64929, 20.37017855 },
	{ 64930, 31.48676870 }, { 64931, 42.22550925 }, { 64932, 19.33828104 }, { 64933, 33.79106328 }, { 64935, 40.70319064 }, { 64936, 45.23365759 }, { 64937, 46.17474445 },
	{ 64938, 47.44896760 }, { 64939, 44.25182551 }, { 64943, 46.63920049 }, { 64945, 45.89737482 }, { 64946, 45.06022067 }, { 64947, 46.74060751 }, { 64950, 46.68856266 },
	{ 64951, 42.64461330 }, { 64952, 44.93455751 }, { 64959, 32.22686107 }, { 64960, 47.41555840 }, { 64961, 47.08019570 }, { 64962, 45.87728098 }, { 64963, 49.34552149 },
	{ 64964, 43.12358705 }, { 64965, 33.32509155 }, { 64969, 7.91273979 }, { 64970, 43.15591763 }, { 64971, 43.39338151 }, { 64972, 40.03853530 }, { 64973, 20.79153207 },
	{ 64974, 43.83737247 }, { 64975, 40.18494183 }, { 64976, 45.70468721 }, { 64977, 26.19547567 }, { 64981, 41.82513092 }, { 64982, 45.82538111 }, { 64983, 51.28283603 },
	{ 64984, 47.89583327 }, { 64985, 15.49300153 }, { 64986, 34.60435647 }, { 64988, 15.86538372 } };

	map<int, double> flux_map_c12 = { { 65006, 4.95266672 }, { 65007, 51.31765765 }, { 65008, 33.05463002 }, { 65009, 38.60310008 }, { 65010, 47.82890344 }, { 65011, 37.89612365 },
	{ 65012, 50.32880412 }, { 65013, 49.44093586 }, { 65014, 47.72058751 }, { 65015, 46.05756360 }, { 65016, 46.85946322 }, { 65017, 51.22289778 }, { 65018, 47.51206486 }, 
	{ 65019, 18.95614002 }, { 65029, 17.18278852 }, { 65030, 54.10555294 }, { 65031, 1.60429310 }, { 65032, 47.63329556 }, { 65033, 39.76412853 }, { 65034, 21.07356745 }, 
	{ 65035, 35.90009687 }, { 65036, 53.29654619 }, { 65037, 42.94676237 }, { 65038, 50.63151030 }, { 65039, 59.08703394 }, { 65040, 57.26591545 }, { 65041, 68.61143417 }, 
	{ 65042, 71.24192471 }, { 65043, 75.57033827 }, { 65044, 70.56807161 }, { 65046, 71.60824322 }, { 65047, 68.88778089 }, { 65048, 67.14873766 }, { 65049, 72.85338573 }, 
	{ 65050, 66.28122582 }, { 65051, 70.17714127 }, { 65052, 75.42806943 }, { 65053, 64.48869458 }, { 65055, 23.02382368 }, { 65059, 71.19792288 }, { 65096, 68.50947434 },
	{ 65097, 49.78323686 }, { 65107, 52.47671567 }, { 65108, 54.60669002 }, { 65109, 68.47537250 }, { 65110, 56.22178299 }, { 65111, 68.65845379 }, { 65112, 40.23722546 } };

	_flux_map.push_back(flux_map_si);
	_flux_map.push_back(flux_map_c12);

	int iarg = 1;
	while (argc > iarg) {
		if (input_option(argv[iarg], "-b")) {
			iarg++;
			_best_tdiff = 1;
		}
		else if (input_option(argv[iarg], "-f")) {
			iarg++;
			if (argc - iarg < 2) exit(input_err());
			_limits[0] = atof(argv[iarg++]);
			_limits[1] = atof(argv[iarg++]);
		}
		else if (input_option(argv[iarg], "-m")) {
			iarg++;
			if (argc - iarg < 1) exit(input_err());
			_method = atoi(argv[iarg++]);
		}
		else if (input_option(argv[iarg], "-l")) {
			iarg++;
			if (argc - iarg < 2) exit(input_err());
			_xmin = atof(argv[iarg++]);
			_xmax = atof(argv[iarg++]);
		}
        else if (input_option(argv[iarg], "-r")) {
			iarg++;
			if (argc - iarg < 1) exit(input_err());
			if (!isposint(string(argv[iarg]))) exit(input_err());
			_num_runs = atoi(argv[iarg++]);
			if (argc - iarg < _num_runs) exit(input_err());
			for (int i = 0; i < _num_runs; i++) {
				_run.push_back(atoi(argv[iarg]));
				_flux += _flux_map[target()][atoi(argv[iarg++])]*1.e9;
			}
		}
		else if (input_option(argv[iarg], "-t")) {
			iarg++;
			if (argc - iarg < 1) exit(input_err());
			_tdiff_cut = atof(argv[iarg++]);
		}
		else if (input_option(argv[iarg], "-in")) {
			iarg++;
			if (argc - iarg < 1) exit(input_err());
			if (string(argv[iarg]) == "fine") {
                iarg++;
				for (int i = 0; i < nangle + 1; i++) _input_angles.push_back(max_angle/nangle*i);
            }
			else if (string(argv[iarg++]) == "coarse") _input_angles = angles[target()];
			else exit(input_err());
		}
		else if (input_option(argv[iarg], "-mc")) {
			iarg++;
			_mc = 1;
		}
		else if (input_option(argv[iarg], "-wo")) {
			iarg++;
			_outer = 1;
		}
		else if (input_option(argv[iarg], "-btc")) {
			iarg++;
			if (argc - iarg < 1) exit(input_err());
			_btc = atoi(argv[iarg++]);
			if (_btc != 0 && _btc != 1 && _btc != 2) exit(input_err());
            if (_btc == 1) {
                vector<double> temp(6, 0.);
                ifstream in_btcorr_par("./btcorr.dat");
                if (!in_btcorr_par.is_open()) exit(open_err(string("./btcorr.dat")));
                for (int i = 0; i < btcorr_nbins(); i++) {
                    for (int j = 0; j < 6; j++) in_btcorr_par >> temp[j]; 
                    _btcorr_par.push_back(temp);
                }
            }
		}
		else if (input_option(argv[iarg], "-acc")) {
			iarg++;
			_subacc = true;
		}
		else if (input_option(argv[iarg], "-out")) {
			iarg++;
			if (argc - iarg < 1) exit(input_err());
			if (string(argv[iarg]) == "fine") {
                iarg++;
                for (int i = 0; i < nangle + 1; i++) _output_angles.push_back(max_angle / nangle*i);
            }
			else if (string(argv[iarg++]) == "coarse") _output_angles = angles[target()];
			else exit(input_err());
		}
		else if (input_option(argv[iarg], "-elas")) {
			++iarg;
			if (argc - iarg < 2) exit(input_err());
			_elas[0] = atof(argv[iarg++]);
			_elas[1] = atof(argv[iarg++]);
		}
		else if (input_option(argv[iarg], "-invm")) {
			++iarg;
			if (argc - iarg < 2) exit(input_err());
			_invm_cut[0] = atof(argv[iarg++]);
			_invm_cut[1] = atof(argv[iarg++]);
		}
		else if (input_option(argv[iarg], "-npi0")) {
			++iarg;
			_use_Npi0 = true;
		}
		else if (input_option(argv[iarg], "-lelas")) {
			iarg++;
			_lelas = atof(argv[iarg++]);
		}
        else if (string(argv[iarg]) == "-h" || string(argv[iarg]) == "--help") {
            iarg++;
            help();
            exit(1);
        }
		else {
            help();
            exit(input_err());
        }

		if ((_limits[0] >= _limits[1]) || _xmin >= _xmax) exit(input_err());
		if (_method != 1 && _method != 2) {
			cerr << "wrong fitting method, exit" << endl;
			exit(1);
		}
	}
    
    _ta_corr = ta_corr[_mc][target()]; 

	if (_output_angles.size() < _input_angles.size()) _output_angles = _input_angles;
	
    if (!_input_angles.size()) _input_angles = angles[target()];
	if (!_output_angles.size()) _output_angles = angles[target()];

	if (!_flux) _flux = flux_target[_mc][target()];
	if (!_flux_omega)_flux_omega = flux_omega_target[target()];
}

bool UImanager::input_option(char* arg, string opt) {
    return _input_option_list[pair<string, string>(prog, opt)] && string(arg) == opt;
}

void UImanager::help() {
    cout << endl << prog << " [options]" << endl;

    if (prog == "rotatedmass")
        cout << " Construct rotatedmass distributions w/ and w/o the transitional region." << endl;
    if (prog == "rotatedmass_mc")
        cout << " Construct M.C. peak shape for rotated mass w/ and w/o the transitional region." << endl;
    if (prog == "rotatedmass_omg")
        cout << " Construct omega rotated mass w/ and w/o the transitional region." << endl;
    if (prog == "fast_fitmodel_mc")
        cout << " Fast version of rotatedmass_mc, constructing rotatedmass distributions w/ and w/o the transitional region using prepared smaller root files." << endl;
    if (prog == "fast_fitmodel_mc_correct")
        cout << " Construct M.C. peak shape with optimized sigma" << endl;
    if (prog == "alt_model_correct")
        cout << " Find the correction to the M.C. peak sigma" << endl;
    if (prog == "altfit")
        cout << " Fit rotated mass distribution w/ M.C. peak shape or double Gaussian shape" << endl;

    cout << "Options:" << endl <<
                "-h | --help              Print this help" << endl << endl;
    if (input_option("-b", "-b"))
        cout << "-b                       Use best tdiff"  << endl << endl;
    if (input_option("-f", "-f"))
        cout << "-f number number         Set rotated mass distribution limits"  << endl << endl;
    if (input_option("-m", "-m"))
        cout << "-m (1|2)                 Set fitting method. 1: Use M.C. 2: Use double Gaussian shape" << endl << endl;
    if (input_option("-l", "-l"))
        cout << "-l number number         Set rotated mass fitting range"  << endl << endl;
    if (input_option("-r", "-r"))
        cout << "-r integer runnumber_list  Use specific run(s). integer: number of runs. runnumber list: list of run numbers delimited by space"  << endl << endl;
    if (input_option("-t", "-t"))
        cout << "-t number                Set tdiff cut window: [ -number, number ]"  << endl << endl;
    if (input_option("-in", "-in"))
        cout << "-in (coarse|fine)        Set number of angular bins in the input. "  << endl <<
                "                                |Silicon  |  Carbon" << endl <<
                "                         _______|_________|________" << endl <<
                "                         coarse |  21     |    14  " << endl <<
                "                         _______|_________|________" << endl <<
                "                         fine   | 125     |   125  " << endl;
    if (input_option("-mc", "-mc"))
        cout << "-mc                      Prepare realistic mc rotated mass"  << endl << endl;
    if (input_option("-wo", "-wo"))
        cout << "-wo                      Use transitional region"  << endl << endl;
	if (input_option("-btc", "-btc"))
		cout << "-btc (0|1|2)             Subtract 2nd best tdiff"  << endl << endl;
    if (input_option("-acc", "-acc"))
        cout << "-acc                     Subtract accidentals using side bands"  << endl << endl;
    if (input_option("-out", "-out"))
        cout << "-out (coarse|fine)       Set number of angular bins in the output."  << endl <<
                "                                |Silicon  |  Carbon" << endl <<
                "                         _______|_________|________" << endl <<
                "                         coarse |  21     |    14  " << endl <<
                "                         _______|_________|________" << endl <<
                "                         fine   | 125     |   125  " << endl;
    if (input_option("-elas", "-elas"))
        cout << "-elas number number      Set limits for regular elasticity cut as number (default [-100, 100])" << endl << endl;
    if (input_option("-invm", "-invm"))
        cout << "-invm number number      Set limits for regular invm. cut as number (default [-100, 100])" << endl << endl;
    if (input_option("-npi0", "-npi0"))
        cout << "-npi0                    Use hist counts - background fit as the number of pi0s" << endl << endl;
    if (input_option("-lelas", "-lelas"))
        cout << "-lelas number            Set lower rotated elasticity cut as number (default 1.25)" << endl << endl;
}

int UImanager::target(){
	char cwd[1024];
	getcwd(cwd, sizeof(cwd));
	string folder(cwd);
	if (folder.find("carbon") != string::npos) return 1;
	else if (folder.find("silicon") != string::npos) return 0;
	else exit(wrong_folder());
}

string UImanager::input_filename (const string filetype) {
	if (filetype == "data") {
		string fname = string(Form("_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f_tdiff%.2f", input_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas, _tdiff_cut));
		if (_best_tdiff) fname += "_btdiff";
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
		return fname;
	}
	else if (filetype == "omega") {
        string fname = string(Form("_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f_omg", input_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas));
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
        return fname;
    }
	else if (filetype == "mc") {
        string fname = string(Form("_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f", input_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas));
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
        return fname;
    }
	else if (filetype == "fit") {
		string fname = string(Form("_method%d_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f_tdiff%.2f_btc%d", _method, input_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas, _tdiff_cut, _btc));
		if (_best_tdiff) fname += "_btdiff";
		if (_subacc) fname += "_subacc";
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
		return fname;
	}
	else exit(input_err());
}

string UImanager::output_filename (const string filetype) {
	if (filetype == "data") {
		string fname = string(Form("_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f_tdiff%.2f", output_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas, _tdiff_cut));
		if (_best_tdiff) fname += "_btdiff";
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
		return fname;
	}
	else if (filetype == "omega") {
        string fname = string(Form("_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f_omg", output_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas));
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
        return fname;
    }
	else if (filetype == "mc") {
        string fname = string(Form("_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f", output_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas));
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
        return fname;
    }
	else if (filetype == "fit") {
		string fname;
        if (_btc) fname = string(Form("_method%d_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f_tdiff%.2f_btc%d", _method, output_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas, _tdiff_cut, _btc));
        else fname = string(Form("_method%d_nt%d_mt%.2f_o%d_lm%.2f_hm%.2f_le%.2f_tdiff%.2f", _method, output_nbins(), max_angle, _outer, _limits[0], _limits[1], _lelas, _tdiff_cut));
		if (_best_tdiff) fname += "_btdiff";
		if (_subacc) fname += "_subacc";
		if (_elas[0] > -50) fname += string(Form("_e%.3f_%.3f", _elas[0], _elas[1]));
		return fname;
	}
	else exit(input_err());
}

bool isposint(const string &s) {
	if (s.empty() || !isdigit(s[0])) return false;
	char *p;
	strtol(s.c_str(), &p, 10);
	return (*p == 0);
}

bool isequal(float x, float y) {
	return fabs(x - y) <= 1.e-3*fabs(x);
}

int open_err(const string fname) {
	cerr << "Could not open file " + fname << endl;
	return 1;
}

int open_err(const TString fname) {
	cerr << "Could not open file " + fname << endl;
	return 1;
}

int input_err() {
	cerr << "wrong input" << endl;
	return 1;
}

int wrong_folder() {
	cerr << "wrong folder" << endl;
	return 1;
}
