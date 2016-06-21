#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#define MAX_HHITS 1728
#define HEGEN(N) read_mcfile_com_.hegen[N-1]

static FILE *fp;

struct {int mcrun, mcevent; float hegen[64];
        int nout, iout[MAX_HHITS], aout[MAX_HHITS];} read_mcfile_com_;
struct {int mcread_flag, mc_eof_flag, nreadmc;} mcread_stat_com_;

static int first_call = 1;
void mcdone();

int mcfile_read(int iopt) {

  if(first_call) {
    first_call = 0;
    fp = fopen("mc.out", "rb");
    mcread_stat_com_.nreadmc = 0;
    mcread_stat_com_.mc_eof_flag = 0;
  }

  float ftmp;
  int   itmp, nelements;
  int i; size_t result;
/*
      hegen(3)  = current_tid         ! T-counter
      hegen(4)  = current_eid         ! E-channel
      hegen(5)  = bm_mom              ! Beam energy
      hegen(16) = gkin(1,2) ! gamma from pi0 decay with greater energy
      hegen(17) = gkin(2,2)
      hegen(18) = gkin(3,2)
      hegen(25) =  th_x*RADDEG        ! Theta pi0 [deg.]
      hegen(26) = gkin(1,1) ! gamma from pi0 decay with lesser energy
      hegen(27) = gkin(2,1)
      hegen(28) = gkin(3,1)
      hegen(31) = type_generated
*/
    result = fread(&itmp, sizeof(itmp), 1, fp); if(result != 1) { mcdone(); return 2;}
    read_mcfile_com_.mcrun = itmp;


//  for some reason it helps:
    fseek (fp, -4, SEEK_CUR);

    result = fread(&itmp, sizeof(itmp), 1, fp); if(result != 1) { printf("p2\n"); return 1;}
    read_mcfile_com_.mcevent = itmp;

 //   for(i = i; i <= 64; ++i) {  // for omega
    for(i = i; i <= 32; ++i) {  // for pi0
      result = fread(&ftmp, sizeof(ftmp), 1, fp); if(result != 1) { printf("p3\n"); return 1;}
      HEGEN(i) = ftmp;
    }
    result = fread(&nelements, sizeof(nelements), 1, fp); if(result != 1) { printf("p1\n"); return 1;}
    read_mcfile_com_.nout = nelements;

    if(nelements)
      for(i = 0; i < nelements; ++i) {
        result = fread(&itmp, sizeof(itmp), 1, fp); if(result != 1) { printf("p4\n"); return 1;}
        read_mcfile_com_.iout[i] = itmp;
        result = fread(&itmp, sizeof(itmp), 1, fp); if(result != 1) { printf("p5\n"); return 1;}
        read_mcfile_com_.aout[i] = itmp;
      }

  mcread_stat_com_.nreadmc += 1;
  if(mcread_stat_com_.mcread_flag == 1) {
    printf("ERR: read_mcfile: mcevent has not been analyzed\n"); exit(1);
  }
  mcread_stat_com_.mcread_flag = 1;
  return 0;
}

void mcdone() {
    FILE *fpw;
    fpw = fopen("mc.done", "w");
    int nmceventsread = mcread_stat_com_.nreadmc;
    fprintf(fpw,"%i",nmceventsread);
    fclose(fpw);
    mcread_stat_com_.mc_eof_flag = 1;
    fclose(fp);
}

int main(){ 
    cout << mcread_stat_com_.mcread_flag << endl; 
    int mcstat = mcfile_read(0);
    cout << mcread_stat_com_.mcread_flag << endl; 
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

    for(int i = 1; i <= 32; ++i) {
        output << setprecision(6) << fixed << setw(12) << HEGEN(i);
    }
    output << endl;

    while(!mcread_stat_com_.mc_eof_flag) {
        if(mcread_stat_com_.mcread_flag == 0) exit(0);
        mcread_stat_com_.mcread_flag = 0;
        int mcstat = mcfile_read(0);
        if(mcstat == 2) {
            //cout << "hycal error" << endl;
            exit(2);
        }
        else if(mcstat == 1) {
            cout << "error reading mc file" << endl;
            exit(1);
        }
        for(int i = 1; i <= 32; ++i) {
            output << setprecision(6) << fixed << setw(12) << HEGEN(i);
        }
        output << endl;
    }
    output.close();
}
