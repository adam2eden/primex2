#include "primex.h"

using namespace std;

int main(int argc, char* argv[]){
    UImanager uimanager(argc, argv, "plotalt");

	const double bw = 0.00034906585; //0.02 degree in radian

	TFile *inroot = new TFile (TString("altfit" + uimanager.input_filename("fit") + ".root"));

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
	fitresult->SetBranchAddress("NDF",&NDF);
	fitresult->SetBranchAddress("Npfits",&Npfits);

	vector<vector<double> > geoeff_targets = { { 0.510706, 0.509669, 0.508074, 0.506536, 0.50527, 0.504297, 0.50359, 0.503009, 0.502535, 0.502082,
        0.501583, 0.500996, 0.500422, 0.499869, 0.499222, 0.498433, 0.497585, 0.496783, 0.496044, 0.495208, 0.494273, 0.493272, 0.492179, 0.49114,
        0.489967, 0.488772, 0.487479, 0.486143, 0.484726, 0.483251, 0.481678, 0.48002, 0.478292, 0.476577, 0.47477, 0.472965, 0.471089, 0.469152,
        0.467128, 0.465008, 0.462828, 0.460588, 0.458341, 0.456033, 0.453635, 0.451174, 0.448623, 0.44602, 0.443424, 0.440772, 0.438094, 0.435327,
        0.432536, 0.429667, 0.426682, 0.42363, 0.420574, 0.417492, 0.414381, 0.411191, 0.407877, 0.404516, 0.401183, 0.397704, 0.394199, 0.390755,
        0.387221, 0.383588, 0.379894, 0.376017, 0.372109, 0.368071, 0.363972, 0.359813, 0.35556, 0.351178, 0.34673, 0.342182, 0.337569, 0.332832,
        0.327972, 0.323199, 0.318216, 0.3131, 0.307872, 0.302542, 0.2971, 0.291622, 0.285934, 0.280197, 0.274444, 0.26856, 0.262507, 0.256315,
        0.250046, 0.243737, 0.237307, 0.230723, 0.224063, 0.217292, 0.210484, 0.203652, 0.196653, 0.18959, 0.182469, 0.175417, 0.168371, 0.161247,
        0.154333, 0.147553, 0.14084, 0.134343, 0.128073, 0.122109, 0.116408, 0.110912, 0.10575, 0.100794, 0.0959735, 0.0913469, 0.0867613, 0.0823831, 
        0.0781501, 0.0740457, 0.0700232 }, 
        { 0.510704, 0.509667, 0.508074, 0.50654, 0.505279, 0.504313, 0.503618, 0.503052, 0.502591, 0.50215, 0.501665, 0.501095, 
		0.50054, 0.500006, 0.499376, 0.498604, 0.497775, 0.496993, 0.496272, 0.495453, 0.494535, 0.493552, 0.492475, 0.491453, 0.490296, 0.489116, 
		0.487839, 0.48652, 0.48512, 0.483661, 0.482105, 0.480464, 0.478752, 0.477056, 0.475268, 0.473483, 0.471628, 0.469712, 0.467709, 0.465609, 
		0.463452, 0.461237, 0.459014, 0.456731, 0.454359, 0.451926, 0.449403, 0.446828, 0.444264, 0.441645, 0.439001, 0.436273, 0.433521, 0.430691, 
		0.427746, 0.424734, 0.421723, 0.418687, 0.415626, 0.412486, 0.409226, 0.405919, 0.40264, 0.399216, 0.395774, 0.392394, 0.388929, 0.385365, 
		0.381736, 0.377926, 0.374088, 0.370126, 0.366106, 0.362027, 0.357853, 0.35355, 0.349186, 0.344725, 0.340191, 0.335537, 0.330762, 0.326072, 
		0.321174, 0.316144, 0.311007, 0.305759, 0.300407, 0.295017, 0.289411, 0.283756, 0.278083, 0.272278, 0.266306, 0.260187, 0.253988, 0.24775, 
		0.241384, 0.234856, 0.228249, 0.221512, 0.214741, 0.207942, 0.20095, 0.193889, 0.186745, 0.179647, 0.172515, 0.165283, 0.158234, 0.151289, 
		0.144374, 0.137659, 0.131165, 0.124978, 0.119042, 0.113315, 0.107939, 0.102772, 0.0977441, 0.0929101, 0.088125, 0.0835484, 0.0791276, 
		0.0748382, 0.0706401 } };

    ifstream infile;

    TString workdir(getenv("WORKDIR")!=NULL?getenv("WORKDIR"):"");
    TString index;
    //if (!uimanager.target()) index = workdir + Form("/tables/effrec%d", int(uimanager.isouter())) + uimanager.input_filename("mc") + "_si.dat";
    //else index = workdir + Form("/tables/effrec%d", int(uimanager.isouter())) + uimanager.input_filename("mc") + "_c12.dat";
    if (!uimanager.target()) index = workdir + Form("/tables/effrec%d_%.2f.dat", int(uimanager.isouter()), uimanager.get_lelas_cut());
    else index = workdir + Form("/tables/effrec%d_%.2f.dat", int(uimanager.isouter()), uimanager.get_lelas_cut());

    vector<double> geoeff;
    infile.open(index);
    if (!infile.is_open()) {
        cerr << "Couldn't open efficiency file " + index << ". Use default efficiency." << endl;
        geoeff = geoeff_targets[uimanager.target()];
    }
    else {
        cout << "Use " + index << endl;
        double temp;
        for (int i = 0; i < 125; i++) {
            infile >> temp;
            geoeff.push_back(temp);
        }
    }

	vector<vector<double> > cx_ilarin = { { 18.91780, 48.57532, 64.94920, 61.54921, 57.61427, 50.52202, 43.77200, 36.99528, 32.76882, 30.89266, 29.91695, 
		27.20740, 25.01229, 23.45575, 24.28930, 24.66808, 23.39777, 23.40127, 24.37760, 22.84130, 22.39708, 24.45628, 23.75674, 25.25551, 23.48691, 
		24.96391, 25.25141, 26.93284, 27.08244, 29.16905, 28.75519, 30.66874, 30.88640, 32.42720, 33.40725, 35.04254, 36.07894, 34.25146, 38.18644, 
		38.73247, 39.29086, 41.97525, 42.35574, 43.81289, 41.25137, 44.87543, 47.82605, 48.68939, 46.90195, 50.47281, 49.10114, 50.60778, 51.70198, 
		54.87564, 53.82585, 56.58145, 55.04273, 55.70100, 57.16637, 56.41164, 56.85794, 59.41360, 57.79117, 58.30198, 56.27025, 59.25443, 59.51877, 
		55.84207, 57.07999, 59.60442, 58.42272, 56.36616, 58.40206, 57.90853, 58.02875, 58.04807, 56.67309, 57.18576, 52.49699, 54.86701, 54.77655, 
		51.25727, 51.02219, 52.57930, 51.16555, 50.52422, 47.71101, 47.01061, 45.82554, 44.37566, 48.37138, 42.48201, 43.55245, 42.15264, 42.47041, 
		38.45493, 39.48478, 40.30909, 37.57168, 35.65630, 35.65940, 33.06170, 35.50649, 36.37715, 34.46681, 31.25928, 30.95740, 30.41431, 29.07522, 
		30.41549, 27.78221, 28.19196, 26.11155, 23.56880, 28.13475, 24.96551, 23.68807, 20.50472, 24.66277, 21.98660, 22.96906, 20.98963, 24.68886, 
		20.03215, 21.52606},
        { 3.442685, 8.998088, 11.985926, 11.760904, 10.105563, 8.290399, 7.701203, 6.641514, 5.660134, 5.283175, 5.496669, 
	    4.889200, 4.816358, 5.065584, 4.443582, 4.641646, 4.489931, 4.479843, 4.313957, 4.000865, 4.458523, 4.538793, 4.913742, 4.650720, 5.230009, 
	    5.642178, 5.910239, 5.584293, 5.464587, 6.033989, 6.853351, 7.087946, 7.801963, 7.683512, 7.769440, 8.982615, 8.890630, 9.270883, 8.629159, 
	    8.793224, 10.869985, 9.729795, 11.212657, 10.977981, 11.218017, 12.872381, 12.865654, 13.517530, 14.108892, 15.287876, 14.299283, 15.588189, 
	    14.587325, 16.746712, 16.536885, 16.043178, 17.874668, 18.428420, 18.264817, 18.803316, 18.956656, 19.268667, 20.332585, 20.253016, 21.118889, 
	    22.348028, 21.550949, 22.022933, 20.156333, 22.944878, 23.165395, 24.213000, 24.063725, 23.126102, 23.825259, 22.726505, 23.648168, 24.901193, 
	    24.718809, 25.138071, 25.000181, 24.816033, 25.807962, 25.293346, 25.378858, 25.516854, 25.781723, 25.042582, 26.038254, 26.538504, 25.174299, 
	    25.254950, 26.019678, 26.912126, 26.679864, 27.357095, 26.332152, 27.564930, 25.654239, 25.903864, 26.502541, 26.477217, 25.475794, 25.377203, 
	    24.881931, 25.061843, 25.884003, 25.783331, 25.186736, 25.155585, 23.160978, 23.005803, 22.527945, 24.448109, 24.469584, 23.130083, 22.703105, 
	    22.914310, 20.685884, 22.128808, 20.456825, 21.965954, 22.499514, 21.828200, 21.324870 } };

	vector<vector<double> > cx_ilarin_err = { { 0.68746, 1.12338, 1.37864, 1.28060, 1.37371, 1.18193, 1.12207, 1.03245, 1.00268, 0.86865, 0.83775, 0.95877, 0.88432, 
		0.84231, 0.86179, 0.86837, 0.79753, 0.80793, 0.85457, 0.83606, 0.83026, 0.85744, 0.85014, 0.86832, 0.85432, 0.85763, 0.87817, 0.88653, 0.90112, 
		0.91334, 1.03089, 0.94388, 0.93950, 0.96875, 0.98230, 1.00759, 1.02972, 1.01143, 1.04535, 1.07725, 1.07186, 1.16745, 1.10027, 1.13060, 1.07593, 
		1.12910, 1.17457, 1.18440, 1.12894, 1.20172, 1.12044, 1.21900, 1.34321, 1.28225, 1.26823, 1.23814, 1.25878, 1.29661, 1.32333, 1.32067, 1.32718, 
		1.34126, 1.48121, 1.34475, 1.33522, 1.37467, 1.38634, 1.34597, 1.37194, 1.38500, 1.37532, 1.37119, 1.40541, 1.41072, 1.42736, 1.41945, 1.40715, 
		1.40000, 1.50514, 1.43450, 1.40837, 1.39647, 1.39424, 1.42263, 1.64007, 1.42774, 1.38195, 1.37072, 1.26143, 1.51148, 1.42996, 1.37734, 1.40287, 
		1.40401, 1.42080, 1.45579, 1.38486, 1.46955, 1.30666, 1.34849, 1.54318, 1.37306, 1.36100, 1.46702, 1.46396, 1.48415, 1.52773, 1.35716, 1.50544, 
		1.48529, 1.35834, 1.50214, 1.48665, 1.30828, 1.58723, 1.54125, 1.49948, 1.52540, 1.67621, 1.61567, 1.75518, 1.65149, 1.74814, 1.68332, 1.64013},  
        {0.226512, 0.372092, 0.428863, 0.425156, 0.405393, 0.372563, 0.321133, 0.328444, 0.315022, 0.331758, 0.306464, 
		0.285716, 0.289402, 0.301439, 0.291478, 0.292641, 0.288822, 0.281666, 0.285904, 0.264578, 0.293364, 0.294336, 0.304978, 0.292056, 0.305397, 
		0.316210, 0.305021, 0.275813, 0.319757, 0.333264, 0.347352, 0.377616, 0.335346, 0.362802, 0.367790, 0.401407, 0.351244, 0.399750, 0.395700, 
		0.354994, 0.435035, 0.417638, 0.438833, 0.444297, 0.441750, 0.511580, 0.536677, 0.487552, 0.493685, 0.516638, 0.485513, 0.523724, 0.517968, 
		0.547759, 0.542830, 0.540399, 0.564760, 0.573425, 0.577735, 0.587894, 0.589473, 0.599569, 0.688651, 0.613750, 0.632143, 0.648278, 0.641109, 
		0.614920, 0.623455, 0.666199, 0.673376, 0.686777, 0.693581, 0.685115, 0.693833, 0.678328, 0.831913, 0.717060, 0.716146, 0.734276, 0.729967, 
		0.740515, 0.721553, 0.755665, 0.757691, 0.720484, 0.772406, 0.838873, 0.783602, 0.795628, 0.793007, 0.794527, 0.902065, 0.907811, 0.836117, 
		0.849892, 0.839673, 0.868976, 0.848518, 0.926670, 0.887261, 0.903089, 0.901933, 0.900230, 0.899196, 0.925610, 1.060280, 0.955547, 1.067688, 
		0.978620, 0.966547, 0.980669, 0.927845, 1.030947, 1.046858, 1.059511, 1.069190, 1.099772, 0.947902, 1.134340, 1.132490, 1.175472, 1.115583, 
		1.098927, 1.255856 } };


	vector<double> angle(nangle, 0), receff(nangle, 0), xerror(nangle, 0), cx(nangle, 0), cxerr(nangle, 0), yield(nangle, 0),
	yerror(nangle, 0), mean1(nangle, 0), sigma1(nangle, 0), fraction(nangle, 0), mean2(nangle, 0), sigma2(nangle, 0),
    emean1(nangle, 0), esigma1(nangle, 0), efraction(nangle, 0), emean2(nangle, 0), esigma2(nangle, 0),
    p0(nangle, 0), p1(nangle, 0), p2(nangle, 0), p3(nangle, 0), p0err(nangle, 0), p1err(nangle, 0),
	p2err(nangle, 0), p3err(nangle, 0), polyint(nangle, 0);

    vector<double> cx_per(nangle, 0.), cx_per_err(nangle, 0.);
	double total_cx = 0, err_cx = 0;
	TH1F *hchi2 = new TH1F("hchi2","rot-d mass fit #chi^{2}/NDF",100,0,0);
    hchi2->SetBit(TH1::kCanRebin);

	TH1F *diff1 = new TH1F("diff1","((d#sigma/d#theta)_{yang} - (d#sigma/d#theta)_{ilya})/ (d#sigma/d#theta)_{ilya} #theta < 0.5", 8, 0, 0);
    diff1->SetBit(TH1::kCanRebin);
	TH1F *diff2 = new TH1F("diff2","((d#sigma/d#theta)_{yang} - (d#sigma/d#theta)_{ilya})/ (d#sigma/d#theta)_{ilya} #theta > 0.5", 50, 0, 0);
    diff1->SetBit(TH1::kCanRebin);
	diff1->GetXaxis()->SetTitle("percentage");
	diff2->GetXaxis()->SetTitle("percentage");

	for(int i=0;i<nangle;i++){
		fitresult->GetEntry(i);
		angle[i] = angles;
		xerror[i] = 0;
        if (uimanager.use_Npi0()) {
            yield[i] = Npi0;
            yerror[i] = Npi0_err; 
        } else {
            yield[i] = parameters[0];
            yerror[i] = errors[0];
        }
		hchi2->Fill(chi2/NDF);
		//if(yerror[i]>3*sqrt(yield[i]) || yerror[i] < 0.1*sqrt(yield[i]))yerror[i] = 1.25*sqrt(yield[i]);
		double tocxmicrobarn = uimanager.flux()*uimanager.get_tsi()*uimanager.get_ta_corr()*bw*geoeff[i]*uimanager.get_adcerr()*uimanager.get_br_corr();
		cx[i] = yield[i]/tocxmicrobarn;
		cxerr[i] = yerror[i]/tocxmicrobarn;
        cout<< yield[i] << " " << yerror[i] << " " << cx[i] << " " << cxerr[i] << endl; 
		total_cx += cx[i];
		err_cx += cxerr[i]*cxerr[i];
        
        if (i < 125) {
            cx_per[i] = (cx[i]/cx_ilarin[uimanager.target()][i]-1)*100;
            cx_per_err[i] = 100*(max(cxerr[i],cx_ilarin_err[uimanager.target()][i])/cx_ilarin[uimanager.target()][i]);

            if (i < 25) diff1->Fill(cx_per[i]);
            else diff2->Fill(cx_per[i]);
        }

        if (uimanager.get_method() == 2) {
            mean1[i] = parameters[1];
            sigma1[i] = parameters[2];
            fraction[i] = parameters[3];
            mean2[i] = parameters[4];
            sigma2[i] = parameters[5];

            emean1[i] = errors[1];
            esigma1[i] = errors[2];
            efraction[i] = errors[3];
            emean2[i] = errors[4];
            esigma2[i] = errors[5];
        }

		p0[i] = parameters[7];
		p0err[i] = errors[7];
		p1[i] = parameters[8];
		p1err[i] = errors[8];
		p2[i] = parameters[9];
		p2err[i] = errors[9];
		p3[i] = parameters[10];
		p3err[i] = errors[10];

		polyint[i] = parameters[12];
	}
	total_cx = total_cx*nangle*bw;
	err_cx = sqrt(err_cx)*nangle*bw;
	cout<< endl << total_cx<<" "<<err_cx<<endl;

	TH1F *hyield = new TH1F("hyield","#pi^{0} yield", nangle, 0, max_angle);
	for(int i=0;i<nangle;i++){
		hyield->SetBinContent(i+1,yield[i]);
		hyield->SetBinError(i+1,yerror[i]);
	}
//draw yield
	TGraphErrors *gr = new TGraphErrors(nangle,&angle[0],&yield[0],&xerror[0],&yerror[0]);
	if(!uimanager.target()){
		gr->SetTitle("#pi^{0} silicon yield rotated m_{#gamma#gamma} method");
		gr->GetYaxis()->SetRangeUser(0,3500);
	}
    else {
		gr->SetTitle("#pi^{0} carbon yield rotated m_{#gamma#gamma} method");
		gr->GetYaxis()->SetRangeUser(0,1600);
	}
	gr->SetMarkerStyle(21);
	gr->SetMarkerColor(4);
	gr->GetXaxis()->SetTitle("#theta, (deg)");
	gr->GetYaxis()->SetTitle("dN/d#theta");
//cross section comparison
	TGraphErrors *gr0 = new TGraphErrors(nangle,&angle[0],&cx[0],&xerror[0],&cxerr[0]);
	gr0->SetTitle("#pi^{0} cross section rotated m_{#gamma#gamma} method");
	gr0->SetMarkerStyle(21);
	gr0->SetMarkerColor(4);
	gr0->SetFillColor(0);
	gr0->SetLineColor(kBlue);
	gr0->GetYaxis()->SetTitle("d#sigma/d#theta, (#mubarn/rad)");
	gr0->GetXaxis()->SetTitle("#theta, (deg)");
	gr0->SetMinimum(0);

	TGraphErrors *gr0_old = new TGraphErrors(125,&angle[0],&cx_ilarin[uimanager.target()][0],&xerror[0],&cx_ilarin_err[uimanager.target()][0]);
	gr0_old->SetTitle("#pi^{0} cross section constrict method");
	gr0_old->SetMarkerStyle(21);
	gr0_old->SetMarkerColor(kRed);
	gr0_old->SetFillColor(0);
	gr0_old->SetLineColor(kRed);
	gr0_old->GetYaxis()->SetTitle("d#sigma/d#theta, (#mubarn/rad)");
	gr0_old->GetXaxis()->SetTitle("#theta, (deg)");
	gr0_old->SetMinimum(0);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr0);
    mg->Add(gr0_old);
//cross section comparison persontage
	TGraphErrors *gr0_per = new TGraphErrors(125,&angle[0],&cx_per[0],&xerror[0],&cx_per_err[0]);
	gr0_per->SetMarkerStyle(21);
	gr0_per->SetMarkerColor(kBlue);
	gr0_per->SetLineColor(kBlue);
	gr0_per->SetFillColor(0);
	gr0_per->GetXaxis()->SetTitle("#theta (deg)");
	gr0_per->GetYaxis()->SetTitle("percentage");

    TGraphErrors *gmean1, *gsigma1, *gmean2, *gfraction, *gsigma2;

    if (uimanager.get_method() == 2) {
        gmean1 = new TGraphErrors(nangle, &angle[0], &mean1[0],&xerror[0],&emean1[0]);
        gmean1->SetTitle("mean1 vs. #theta");
        gmean1->SetMarkerStyle(21);
        gmean1->SetMarkerColor(4);

        gsigma1 = new TGraphErrors(nangle,&angle[0],&mean1[0],&xerror[0],&esigma1[0]);
        gsigma1->SetTitle("#sigma_{1} vs. #theta");
        gsigma1->SetMarkerStyle(21);
        gsigma1->SetMarkerColor(4);

        gfraction = new TGraphErrors(nangle,&angle[0],&fraction[0],&xerror[0],&efraction[0]);
        gfraction->SetTitle("fraction vs. #theta");
        gfraction->SetMarkerStyle(21);
        gfraction->SetMarkerColor(4);

        gmean2 = new TGraphErrors(nangle,&angle[0],&mean2[0],&xerror[0],&emean2[0]);
        gmean2->SetTitle("mean2 vs. #theta");
        gmean2->SetMarkerStyle(21);
        gmean2->SetMarkerColor(4);

        gsigma2 = new TGraphErrors(nangle,&angle[0],&mean2[0],&xerror[0],&esigma2[0]);
        gsigma2->SetTitle("#sigma_{2} vs. #theta");
        gsigma2->SetMarkerStyle(21);
        gsigma2->SetMarkerColor(4);
    }

	TGraphErrors *gr1 = new TGraphErrors(nangle,&angle[0],&p0[0],&xerror[0],&p0err[0]);
	gr1->SetTitle("poly. bkg p_{0} vs. #theta");
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerColor(4);

	TGraphErrors *gr2 = new TGraphErrors(nangle,&angle[0],&p1[0],&xerror[0],&p1err[0]);
	gr2->SetTitle("poly. bkg p_{1} vs. #theta");
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(4);

	TGraphErrors *gr3 = new TGraphErrors(nangle,&angle[0],&p2[0],&xerror[0],&p2err[0]);
	gr3->SetTitle("poly. bkg p_{2} vs. #theta");
	gr3->SetMarkerStyle(21);
	gr3->SetMarkerColor(4);
	gr3->GetYaxis()->SetRangeUser(-0.2,0.2);

	TGraphErrors *gr4 = new TGraphErrors(nangle,&angle[0],&p3[0],&xerror[0],&p3err[0]);
	gr4->SetTitle("poly. bkg p_{3} vs. #theta");
	gr4->SetMarkerStyle(21);
	gr4->SetMarkerColor(4);
	gr4->GetYaxis()->SetRangeUser(-0.008,0.008);

	TGraph *gr5 = new TGraph(nangle,&angle[0],&polyint[0]);
	gr5->SetTitle("integal under poly. bkg vs. #theta");
	gr5->SetMarkerStyle(21);
	gr5->SetMarkerColor(4);
//acceptance + rot-d elas. eff.
	//double angle1[135];
	//for(int i=0;i<135;i++)angle1[i] = i*0.02+0.01;
	//TGraph *gr6 = new TGraph(135,angle1,geoeff);
	//gr6->SetTitle("geometric acceptance (+basic energy cut) vs. #theta");
	//gr6->SetMarkerStyle(21);
	//gr6->SetMarkerColor(4);
	//gr6->SetFillColor(0);
	//gr6->SetLineColor(0);
	//gr6->GetYaxis()->SetRangeUser(0,0.65);


	TString filename("pi0alt" + uimanager.output_filename("fit"));
//dN/dt
	TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
	c1->SaveAs(filename+".pdf[");
	c1->cd();
	gr->Draw("ap");
	c1->SaveAs(filename+".pdf");
//cross section comparison
	c1->cd();
	gr0->Draw("ap");
	gr0->SetTitle("cross section Yang");
    gr0->GetXaxis()->SetTitle("#theta, (deg)");
    gr0->GetYaxis()->SetTitle("d#sigma/d#theta, (#mubarn/rad)");
    gr0->SetMinimum(0);
    c1->SaveAs(filename+".pdf");

    c1->cd();
	gr0_old->Draw("ap");
	gr0_old->SetTitle("cross section Ilya/LingLing");
    gr0_old->GetXaxis()->SetTitle("#theta, (deg)");
    gr0_old->GetYaxis()->SetTitle("d#sigma/d#theta, (#mubarn/rad)");
    gr0_old->SetMinimum(0);
    c1->SaveAs(filename+".pdf");

	c1->cd();
	mg->Draw("ap");
    mg->SetTitle("cross section comparison");
    mg->GetXaxis()->SetTitle("#theta, (deg)");
    mg->GetYaxis()->SetTitle("d#sigma/d#theta, (#mubarn/rad)");
    mg->SetMinimum(0);
    gPad->Modified();

	TLegend *leg1;
	leg1 = new TLegend(0.4,0.12,0.78,0.25);
	leg1->AddEntry(gr0,"rotated m_{#gamma#gamma} method");
	leg1->AddEntry(gr0_old,"ILya & LingLing");
	leg1->Draw();
	c1->Update();
	c1->SaveAs(filename+".pdf");

	TCanvas *c3 = new TCanvas("c3", "c3", 1600,1200);
	c3->Divide(2,1);
	c3->cd(1);
	diff1->Draw();
	c3->cd(2);
	diff2->Draw();
	c3->SaveAs(filename+".pdf");

	c1->cd();
	gr0_per->SetTitle("cross section comparison ((d#sigma/d#theta)_{yang} - (d#sigma/d#theta)_{ilya})/ (d#sigma/d#theta)_{ilya}");
    gr0_per->GetYaxis()->SetRangeUser(-50,50);
	gr0_per->Draw("ap");
	c1->Update();
	c1->SaveAs(filename+".pdf");
//acceptance
	TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
	c2->cd();
    hchi2->GetXaxis()->SetTitle("#chi^{2}/NDF");
	hchi2->Draw();
	c2->SaveAs(filename+".pdf");
    if (uimanager.get_method() == 2) {
        c1->cd();
        gmean1->Draw("ap");
        c1->SaveAs(filename+".pdf");
        c1->cd();
        gsigma1->Draw("ap");
        c1->SaveAs(filename+".pdf");
        c1->cd();
        gfraction->Draw("ap");
        c1->SaveAs(filename+".pdf");
        c1->cd();
        gmean2->Draw("ap");
        c1->SaveAs(filename+".pdf");
        c1->cd();
        gsigma2->Draw("ap");
        c1->SaveAs(filename+".pdf");
    }

	c1->cd();
	gr1->Draw("ap");
	c1->SaveAs(filename+".pdf");

	c1->cd();
	gr2->Draw("ap");
	c1->SaveAs(filename+".pdf");
	c1->cd();
	gr3->Draw("ap");
	c1->SaveAs(filename+".pdf");
	c1->cd();
	gr4->Draw("ap");
	c1->SaveAs(filename+".pdf");
	c1->cd();
	gr5->Draw("ap");
	c1->SaveAs(filename+".pdf");
	c1->SaveAs(filename+".pdf"+"]");

	TFile *outroot = new TFile(filename+".root","RECREATE");
	hyield->Write();
    gr->SetName("yield");
	gr->Write();
    gr0->SetName("cross_section");
	gr0->Write();
    gr0_old->SetName("cross_section_ilarin");
	gr0_old->Write();
	gr1->Write();
	gr2->Write();
	gr3->Write();
	gr4->Write();
	gr5->Write();
	//gr6->Write();
	outroot->Close();

	return 0;
}
