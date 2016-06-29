#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

int main()
{
	string workdir(getenv("WORKDIR"));
    
	ifstream feid(workdir + "/tables/eflux.dat", ios::in);
	vector<double> eid(180, 0);
	for(int i = 0; i < 180; ++i) feid >> eid[i];
	feid.close();
	
	vector<vector<vector<double>>> dfp_mc(180, vector<vector<double>>(540, vector<double>(5, 0.)));
	vector<vector<double>> weight(180, vector<double>(540, 0));
	
	ifstream fdfp_mc(workdir + "/tables/dfpmc_si.dat", ios::in);
	ofstream fweight("weight.dat", ios::out);
    double norm = 0;

	for(int i = 0; i < 180; ++i)
	{
		for(int j = 0; j < 540; ++j)
		{
			for(int k = 0; k < 5; ++k)
			{
				int ech, theta, proc;
                double val;
				fdfp_mc >> ech >> theta >> proc >> val;
				if(ech - 1 != i || theta - 1 != j || proc - 1 != k) 
                {
                    cout << "incorrect cross section table" << endl;
                    return 1;
                }
				dfp_mc[i][j][k] = val;
			}
			weight[i][j] = 7.7 * dfp_mc[i][j][0] + sqrt(7.7 * 1.0)
							* (cos(0.8) * dfp_mc[i][j][1] + sin(0.8) * dfp_mc[i][j][4])
							+ 1.0 * dfp_mc[i][j][2] + 0.5 * dfp_mc[i][j][3];
			weight[i][j] *= eid[i];
            norm += weight[i][j];
		}
	}
	for(int i = 0; i < 180; ++i)
	{
		for(int j = 0; j < 540; ++j)
		{
			fweight << setw(3) << i << setw(4) << j << setprecision(5) << scientific << setw(15) << weight[i][j]/norm << endl;
        }
    }
	fweight.close();
	return 0;
}
