#include <iostream>
#include <iomanip>
#include <vector>

#include <TFile>
#include <TH1D>

using namespace std;

int main()
{
	TFile outroot("eff_mc.root", "RECRATE");
	ifstream input("weight.dat", ios::in);
	vector<double> weight(180, vector<double>(540, 0));
	for(int i = 0; i < 180; ++i)
	{
		for(int j = 0; j < 540; ++j)
		{
			int ech, mc;
			input >> ech >> mc >> weight[i][j];
			if(ech != i || mc != j)
			{
				cout << "incorrect weighting file" << endl;
				return 1;
			}
		}
	}
	input.close();
	
	vector<vector<double>> eff(180, vector<double>(540, 0));
	vector<double> norm(540, 0);
	vector<double> eff_mc(540, 0);
	input.open("Nmc_Nacc.dat", ios::in);
	for(int i = 0; i < 180; ++i)
	{
		for(int j = 0; j < 540; ++j)
		{
			int ech, mc, Nmc, Nacc;
			input >> ech >> mc >> Nmc >> Nacc;
			if(ech != i + 1 || mc != j + 1)
			{
				cout << "incorrect eff file" << endl;
				return 1;
			}
			eff[i][j] = Nacc * 1./ Nmc;
		}
	}
	
	for(int i = 0; i < 540; ++i)
	{
		for(int j = 0; j < 180; ++j)
		{
			eff_mc[i] += eff[i][j] * weight[i][j];
			norm[i] += weight[i][j];
		}
		eff_mc[i] /= norm[i];
	}
	
	TH1D *eff_mc = new TH1D("eff_mc", "eff_mc", 540, 0, 2.7);
	for(int i = 0; i < 540; ++i) eff_mc->SetBinContent(i + 1, eff_mc[i]);
	outroot.cd();
	eff_mc->Write();
	outroot.Close();
	
	return 0;
}