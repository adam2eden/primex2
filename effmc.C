#include "primex.h"
#define mpi0 0.134976

using namespace std;

int main(int argc, char *argv[]){
	gROOT->ProcessLine("#include <vector>");

	if(argc < 3){cout<<"wrong number of inputs"<<endl;return 1;}

	int echn = atoi(argv[1]);
	int ifile = atoi(argv[2]);
	const int kMax = 100;
	int runnumber, eventid, nph, tid[kMax], eid[kMax], npi0, type1[kMax], type2[kMax], id1[kMax], id2[kMax];
	float Etot_cl, Etot, phE[kMax], hx1[kMax], hy1[kMax], hx2[kMax], hy2[kMax], E1[kMax], E2[kMax], invm[kMax], \
        angle[kMax], MC_angle, MC_z, MC_x1, MC_x2, MC_y1, MC_y2, MC_E1, MC_E2;

    TString filename(Form("ge%d.cout",(600+echn)*100+ifile));
	TFile *inroot = new TFile(("./mc/"+filename+".root"));
	TTree *event = (TTree*)inroot->Get("event");
	event->SetBranchAddress("runnumber", &runnumber);
	event->SetBranchAddress("eventid", &eventid);
	event->SetBranchAddress("Etot_cl",&Etot_cl);
	event->SetBranchAddress("Etot",&Etot);
	event->SetBranchAddress("nph", &nph);
	event->SetBranchAddress("phE", phE);
	event->SetBranchAddress("tid", tid);
	event->SetBranchAddress("eid", eid);
	event->SetBranchAddress("npi0", &npi0);
	event->SetBranchAddress("type1", type1);
	event->SetBranchAddress("type2", type2);
	event->SetBranchAddress("id1", id1);
	event->SetBranchAddress("id2", id2);
	event->SetBranchAddress("hx1", hx1);
	event->SetBranchAddress("hx2", hx2);
	event->SetBranchAddress("hy1", hy1);
	event->SetBranchAddress("hy2", hy2);
	event->SetBranchAddress("E1", E1);
	event->SetBranchAddress("E2", E2);
	event->SetBranchAddress("invm",invm);
	event->SetBranchAddress("angle",angle);
	event->SetBranchAddress("MC_angle",&MC_angle);
	event->SetBranchAddress("MC_z",&MC_z);
	event->SetBranchAddress("MC_x1",&MC_x1);
	event->SetBranchAddress("MC_x2",&MC_x2);
	event->SetBranchAddress("MC_y1",&MC_y1);
	event->SetBranchAddress("MC_y2",&MC_y2);
	event->SetBranchAddress("MC_E1",&MC_E1);
	event->SetBranchAddress("MC_E2",&MC_E2);

    TString dirname("./pi0/");
	ofstream outfile0;
    outfile0.open(dirname+filename+".mat");//inv. mass method
	ofstream outfile1;
    outfile1.open(dirname+filename+".mat1");//hybrid mass method
	ofstream theta_stat;
    theta_stat.open(dirname+filename+".theta_stat");//generator stat.
    if(!outfile0||!theta_stat){
        cerr << "can not open output file"<<endl;
        exit(1);
    }

	vector<int> MC_Npi0(640,0);
    const double DEGRAD = 0.0174532925199432958;
	int nevents = event->GetEntries();
	cout<<nevents<<endl;
	for(int i=0;i<nevents;i++){
		//if(i%100000==0)cout<<i<<endl;
		event->GetEntry(i);
        int cc = int(MC_angle/0.005);
        if(cc<0||cc>=640)continue;
		MC_Npi0[cc]++;
        for(int j=0;j<npi0;j++){
            // loose MC matching: | rec-d coord - projected MC coordinate |   < 4cm && | 1 - (rec-d energy / original MC energy) | < 0.3

            int match1 = 0, match2 = 0;
            if(fabs(hx1[j]-MC_x1) < 4. && fabs(hx2[j]-MC_x2) < 4. && fabs(hy1[j]-MC_y1) < 4. && fabs(hy2[j]-MC_y2) < 4. && fabs(1-E1[j]/MC_E1) < 0.3 && fabs(1-E2[j]/MC_E2) < 0.3) match1 = 1;
            if(fabs(hx2[j]-MC_x1) < 4. && fabs(hx1[j]-MC_x2) < 4. && fabs(hy2[j]-MC_y1) < 4. && fabs(hy1[j]-MC_y2) < 4. && fabs(1-E2[j]/MC_E1) < 0.3 && fabs(1-E1[j]/MC_E2) < 0.3) match2 = 1;
            
            if (!match1 && !match2) continue;
            //energy cuts
            if(E1[j]<0.5||E2[j]<0.5)continue;
            if(E1[j]+E2[j]<3.5 || E1[j]+E2[j]>8)continue;
    
            //fiducial cuts: center hole area
            if(id1[j]<1000||id2[j]<1000)continue; //glass cut off
            int irow1, irow2, icol1, icol2;
            irow1 = (id1[j]-1001)/34+1;
            icol1 = (id1[j]-1001)%34+1;
            if(icol1>=16 && icol1<=19 && irow1>=16 && irow1<=19) continue;
            irow2 = (id2[j]-1001)/34+1;
            icol2 = (id2[j]-1001)%34+1;
            if(icol2>=16 && icol2<=19 && irow2>=16 && irow2<=19) continue;
            //mark transitional area
            int isouterlayer = (icol1 == 1) || (icol1 == 34) || (irow1 == 1) || (irow1 == 34) || (icol2 == 1) || (icol2 == 34) || (irow2 == 1) || (irow2 == 34);
            
            //reconstructed angular bin            
            int cc1 = int(angle[j]/0.02);
            if(cc1<0||cc1>=160)continue;

            float de1 = 0., de2 = 0., dx1 = 0., dy1 = 0., dx2 = 0., dy2 = 0.;
            float e1 = 0., e2 = 0., x_1 = 0., x_2 = 0., y_1 = 0., y_2 = 0.;
            if(match1) {
                e1 = E1[j]; e2 = E2[j]; x_1 = hx1[j]; y_1 = hy1[j]; x_2 = hx2[j]; y_2 = hy2[j];
                de1 = E1[j] - MC_E1;
                de2 = E2[j] - MC_E2;
                dx1 = hx1[j] - MC_x1;
                dy1 = hy1[j] - MC_y1;
                dx2 = hx2[j] - MC_x2;
                dy2 = hy2[j] - MC_y2;
            }
            else if (match2) {
                e1 = E2[j]; e2 = E1[j]; x_1 = hx2[j]; y_1 = hy2[j]; x_2 = hx1[j]; y_2 = hy1[j];
                de1 = E2[j] - MC_E2;
                de2 = E1[j] - MC_E1;
                dx1 = hx2[j] - MC_x2;
                dy1 = hy2[j] - MC_y2;
                dx2 = hx1[j] - MC_x1;
                dy2 = hy1[j] - MC_y1;
            }

            float cos45 = 0.70710678119;
            float sin45 = 0.70710678118;

           /* e1 = E1[j] - alpha * de1;
            e2 = E2[j] - alpha * de2;
            x_1 = hx1[j] - beta * dx1;
            y_1 = hy1[j] - beta * dy1;
            x_2 = hx2[j] - beta * dx2;
            y_2 = hy2[j] - beta * dy2;
*/
            float Epi0 = e1 + e2;
            float elas = Epi0 / (MC_E1 + MC_E2);
            float ct = (x_1 * x_2 + y_1 * y_2 + 702 * 702) / sqrt(x_1 * x_1 + y_1 * y_1 + 702 * 702) / sqrt(x_2 * x_2 + y_2 * y_2 + 702 * 702 );
            if(ct >= 1) continue;
            float m0 = sqrt(2 * e1 * e2 * (1 - ct));
            float r1 = sqrt(x_1 * x_1 + y_1 * y_1 + 702 * 702);
            float r2 = sqrt(x_2 * x_2 + y_2 * y_2 + 702 * 702);
            float px = e1 * x_1 / r1 + e2 * x_2 / r2; 
            float py = e1 * y_1 / r1 + e2 * y_2 / r2; 
            float pz = e1 * 702 / r1 + e2 * 702 / r2; 
            float angle_c = atan(sqrt(px * px + py * py) / pz) / DEGRAD;

            float altinvm = m0/mpi0*cos45 - elas*sin45;
            float altelas = m0/mpi0*sin45 + elas*cos45; 
            
            outfile0 << echn << " " << cc+1 << " " << cc1+1 << " " << isouterlayer << " " << MC_angle << " " << angle[j] << " " << angle_c << " " << altinvm << " " << altelas << " " << m0 << " " << elas << endl;
            outfile1 << echn << " " << isouterlayer << " " << MC_angle << " " << MC_E1 << " " << MC_E2 << " " << MC_x1 << " " << MC_y1 << " " << MC_x2 << " " << MC_y2 << " " << e1 << " " << e2 << " " << x_1 << " " << y_1 << " " << x_2 << " " << y_2 << endl;
		}
    }

	for(int i=1;i<=640;i++) theta_stat << i << " " << MC_Npi0[i-1] << endl;
    
    return 0;
}
