#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#define MAX_HHITS 1728

using namespace std;

extern "C" struct {int mcrun, mcevent; float hegen[64];
                int nout, iout[MAX_HHITS], aout[MAX_HHITS];} read_mcfile_com_;
extern struct {int mcread_flag, mc_eof_flag, nreadmc;} mcread_stat_com_;
extern "C" int mcfile_read(int iopt);
    
int main(){ 
    int mcstat = mcfile_read(0);
    if(mcstat == 2) {
        cout << "hycal error" << endl;
        exit(2);
    }
    else if(mcstat == 1) {
        cout << "error reading mc file" << endl;
        exit(1);
    }
    char fname[256];
    sprintf(fname, "ge%d.cout.dat", read_mcfile_com_.mcrun);
    ofstream output(fname, ofstream::out);

    for(int i = 0; i < 32; ++i) {
        output << setprecision(6) << fixed << setw(12) << read_mcfile_com_.hegen[i];
    }
    output << endl;

    while(!mcread_stat_com_.mc_eof_flag) {
        int mcstat = mcfile_read(0);
        if(mcstat == 2) {
            cout << "hycal error" << endl;
            exit(2);
        }
        else if(mcstat == 1) {
            cout << "error reading mc file" << endl;
            exit(1);
        }
        for(int i = 0; i < 32; ++i) {
            output << setprecision(6) << fixed << setw(12) << read_mcfile_com_.hegen[i];
        }
        output << endl;
    }
    output.close();
}
