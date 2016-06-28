#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

using namespace std;

int main(int argc, char* argv[]){
    vector<double> eid = {0.0024208768, 0.0013084378, 0.0080143990, 0.0028857995, 0.0063417418, 0.0026269207, 0.0059270122, 0.0030106680, 0.0063838612, 0.0031747853, 0.0068374140, 0.0032393836, 0.0057597513, 0.0028512922, 0.0055873314, 0.0036027598, 0.0070705900, 0.0041495445, 0.0061263916, 0.0035052080, 0.0062953467, 0.0037479553, 0.0071386490, 0.0034913690, 0.0053078084, 0.0035455046, 0.0084149965, 0.0034658228, 0.0063346998, 0.0037505208, 0.0087057448, 0.0023106382, 0.0048839772, 0.0022146080, 0.0086812465, 0.0022498460, 0.0049359289, 0.0026215449, 0.0097910838, 0.0030694448, 0.0057011592, 0.0023492927, 0.0066916243, 0.0033049374, 0.0090496372, 0.0041227858, 0.0062285564, 0.0035706614, 0.0086107986, 0.0039306612, 0.0071202812, 0.0039015458, 0.0076017892, 0.0040982112, 0.0086700492, 0.0028680902, 0.0054940915, 0.0033979283, 0.0100612000, 0.0050620340, 0.0068545612, 0.0048460092, 0.0088548628, 0.0038326959, 0.0060979868, 0.0040158675, 0.0082774170, 0.0066623804, 0.0082180099, 0.0057133803, 0.0067992773, 0.0043726112, 0.0086418612, 0.0037498824, 0.0063548702, 0.0044739530, 0.0098256272, 0.0051580481, 0.0066690851, 0.0036691487, 0.0080661980, 0.0040915626, 0.0106631361, 0.0040991265, 0.0061937359, 0.0026620984, 0.0073910720, 0.0036334209, 0.0100982487, 0.0030275222, 0.0060826101, 0.0029543403, 0.0093005265, 0.0025709382, 0.0061023870, 0.0024706965, 0.0084049675, 0.0034990533, 0.0091763726, 0.0055113271, 0.0095589677, 0.0041280773, 0.0055169117, 0.0039466803, 0.0107112255, 0.0050447502, 0.0074413133, 0.0044523413, 0.0070669365, 0.0055834491, 0.0104544826, 0.0037378741, 0.0078883942, 0.0060580595, 0.0090152022, 0.0046068552, 0.0058981738, 0.0035763625, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0108579145, 0.0049906186, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0096771154, 0.0020747040, 0.0065395393, 0.0034850938, 0.0116043957, 0.0024547296, 0.0075223201, 0.0036390055, 0.0089209425, 0.0050220103, 0.0091729961, 0.0038423475, 0.0070865167, 0.0040228251, 0.0096624614, 0.0049473269, 0.0108678150, 0.0042415759, 0.0066914476, 0.0034631168, 0.0093252136, 0.0072740644, 0.0100549450, 0.0055336092, 0.0078818702, 0.0054673729, 0.0087701664, 0.0036696746, 0.0109633875, 0.0025023653, 0.0050732393, 0.0021195253, 0.0102841584, 0.0057396692, 0.0105363124, 0.0059488728, 0.0146054764, 0.0000000000, 0.0000000000, 0.0000000000, 0.0167939520, 0.0061681174, 0.0090204535, 0.0060559558, 0.0105734373, 0.0058989005, 0.0069586290, 0.0051160933, 0.0081325989, 0.0010565207, 0.0002536957, 0.0000457327};

    TString workdir(getenv("WORKDIR"));
    //float lelas = 1.25;
    //if(argc == 2) lelas = atof(argv[1]);
   
    int nmc = 540;
    int nrec = 135;

    const int nech = 180;
	ifstream events0, events1, stat, eflux, profile;

    vector<int> theta_stat(nmc, 0), count0(nmc, 0), count1(nmc,0);
    vector<double> theta_stat1(nmc, 0), count01(nmc, 0), count11(nmc,0);

	int echn, status, mc, count;
    float MC_angle, MC_E1, MC_E2, MC_x1, MC_y1, MC_x2, MC_y2, e1, e2, x_1, y_1, x_2, y_2;
    const float mc_mass_outer = 0.134344, mc_elas_outer = 0.995805;
    const float costheta = 0.70710678118, sintheta = 0.70710678118;

	const double DEGRAD = 0.0174532925199432958;
    TFile eff_mct_large("eff_mct_large.root", "RECREATE");
    TTree *ev = new TTree("ev", "M.C. accepted events");
    TTree *ge = new TTree("ge", "M.C. generated events");
    int _count = 0;
    float weid = 0, theta_mc = 0;
    ev->Branch("echn", &echn, "echn/I");
    ev->Branch("weid", &weid, "weid/F");
    ev->Branch("theta_mc", &MC_angle, "theta_mc/F");
    ge->Branch("echn", &echn, "echn/I");
    ge->Branch("weid", &weid, "weid/F");
    ge->Branch("theta_mc", &MC_angle, "theta_mc/F");
    ge->Branch("count", &_count, "count/I");
    	
	or(int i=1;i<=180;i++){
		for(int j=1;j<=20;j++){
            events0.open(Form(workdir + "efficiency2/pi0/ge%d.cout.mat1",(i+600)*100+j));
			stat.open(Form(workdir + "efficiency2/pi0/ge%d.cout.theta_stat",(i+600)*100+j));
			while(1){
				events0>>echn>>status>>MC_angle>>MC_E1>>MC_E2>>MC_x1>>MC_y1>>MC_x2>>MC_y2>>e1>>e2>>x_1>>y_1>>x_2>>y_2;
				if(events0.eof())break;
				int mc = int(MC_angle / 0.005) + 1;
				if (mc > nmc) continue;
                
                count1[mc - 1]++;
                count11[mc - 1] += eid[echn - 1];
                if(!status) {
                    count0[mc - 1]++;
                    count01[mc - 1] += eid[echn - 1];
                }
                theta_mc = MC_angle;
                weid = eid[echn - 1];
                ev->Fill();
            }
			while(1){
				stat>>mc>>count;
				if(stat.eof())break;
                if(mc>nmc)continue;
				theta_stat[mc-1] += count;
				theta_stat1[mc-1] += count * eid[echn - 1];
                _count = count;
                theta_mc = -0.005/2 + 0.005 * mc;
                echn = i;
                weid = eid[echn - 1];
                ge->Fill();
			}
			events0.close();
			stat.close();
            if(j==1)cout<<Form("ge%d.cout",(i+600)*100+j)<<endl;
		}
	}

    ofstream output("eff_mct_large.dat", ios::out);
    TH1D eff_mct("eff_mct", "eff_mct", nmc, 0, 2.7);
    TH1D eff_mct1("eff_mct1", "eff_mct normalized by eid", nmc, 0, 2.7);
    for(int i = 0; i < nmc; ++i)
    {
        eff_mct.SetBinContent(i + 1, count1[i]*1./theta_stat[i]);
        eff_mct.SetBinError(i + 1, sqrt(theta_stat[i])/theta_stat[i]*count1[i]*1./theta_stat[i]);
        eff_mct1.SetBinContent(i + 1, count11[i]*1./theta_stat1[i]);
        eff_mct1.SetBinError(i + 1, sqrt(theta_stat[i])/theta_stat[i]*count11[i]*1./theta_stat1[i]);
        output << setprecision(3) << fixed << setw(5) << i*0.005 << setprecision(7) << fixed << setw(10) << count1[i]*1./theta_stat[i] << setprecision(7) << fixed << setw(10) << sqrt(theta_stat[i])/theta_stat[i]*count1[i]*1./theta_stat[i] << endl;
    }
    eff_mct.Write();
    eff_mct1.Write();
    ev->Write();
    ge->Write();
    eff_mct_large.Close();
    output.close();
	return 0;
}
