/*
RBS:
$ g++ -g -pthread -m64 -Wno-deprecated -I${ROOTSYS}/include -o oscill_bin oscill_bin.C `root-config --cflags --libs` -lMinuit

if sampling is high allocate memory on heap and delete


*/

/*
  Name:           read_binary.cpp
  Created by:     Stefan Ritt <stefan.ritt@psi.ch>
  Date:           July 30th, 2014
  Purpose:        Example file to read binary data saved by DRSOsc.
 
  Compile and run it with:
 
  gcc -o read_binary read_binary.cpp
 
  ./read_binary <filename>
  This program assumes that a pulse from a signal generator is split
  and fed into channels #1 and #2. It then calculates the time difference
  between these two pulses to show the performance of the DRS board
  for time measurements.
  $Id: read_binary.cpp 21495 2014-09-26 14:20:49Z ritt $
*/
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdint.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <TMath.h>
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPostScript.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TStyle.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <TF2.h>
#include <TH1.h>
#include "TStyle.h"
#include <string.h>
#include <sstream>
#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fcntl.h>
#include <unistd.h>
typedef struct {
  char           unknown[11];
  char           blockname[16];
  char			  temp[15];
  unsigned short comtype;
  unsigned short comorder;
}THEADER;
//ext
//int		      wv_des[10];.
typedef struct {
  char           inst_name[16];//
  char 		  inst_numb[4];
  char           tracelabel[16];//
  char           words1[4];
}UHEADER;
//int           wv_array_count[9];
//char words3[4];  
//float          vert_gain;
//float		  vert_off;
//float		  max;
//float		  min;
//char words4[4];
//float		  h_interval;
//double         offset1;
//double         pix_off;
typedef struct {
  char           defn[100];
  char		  time_stmp_sec[16];
  char          acq_dur[4];
  unsigned	short recordtype;
  unsigned	short processing;
  char           words4[4];
  unsigned	short timebase;
  unsigned	short vert_coup;
  char		  prob_atten[4];
  unsigned	short vert1_gain;
  unsigned	short bandwidth;
  char          vert_vernier[8];
  unsigned	short wavesource;
} BHEADER;

using namespace std;
int main()
{
  THEADER th;
  BHEADER bh;
  UHEADER uh;
  char words3[4];
  char words4[4];
  char ext;
  int		      wv_des[10];
  int           wv_array_count[9];
  float		  h_interval;
  double         hor_off;
  double         pix_off;  
  float          vert_gain;
  float		  vert_off;
  float		  max;
  float		  min;
  const int n_chnls = 1;   
  const int samples = 502;// 500MS/s window: 500ns 1ms: 10GS/S
   
   int16_t voltage[samples];
   double waveform[samples], time[samples];

  // //allocate memory on heap since sample size is more
  // int16_t* voltage = new int16_t[samples];
  // double* waveform = new double[samples];
  // double* time = new double[samples];

   int i;
   char filename[256];
   char name[300];
   const int N_events = 9999;
  
   const double tenp3 = 1000.;
   const double tenp9 = pow(10,9);
   const int nhv=1;
   int HV1[nhv] = {0};
   int HV2[nhv] = {0};



  double peak_ht[N_events];
  double time_val[N_events];
  double chrg_val[N_events]; 
  TH1F *raw_evt[N_events];


  char infile11[200];  
  char infile[200];
  int folder_number=0;
  char title[100];
  char foldername[100];
  char logfilename[100];
  char textfilename[100];

  ifstream voltage_folder("folders.log");
  int evt = 0;
  int ch = 0;
  int hv_no=0;
  
  while(!(voltage_folder.eof())){
    cout <<"mmmmmmmm : "<<endl;
    voltage_folder >> foldername >> logfilename;
    if(voltage_folder.eof()) break;
    folder_number++;
    ifstream comfile_db;
    sprintf(infile11,"%s",foldername);
    //  cout <<"First : "<< infile11<< endl;      
    comfile_db.open(logfilename);
    evt = 0;
	  
    sprintf(name, "%s.root", foldername);
    TFile *newfile = new TFile(name,"recreate");
    //  TH1F *chnl_raw[N_events];

   
    while(!(comfile_db.eof())){
      comfile_db >> textfilename;
      if (strstr(textfilename,"#")) continue;

      if(comfile_db.eof()) break;
      //sprintf(infile,"%s/%s",FOLDER, textfilename);
      sprintf(infile,"%s/%s",foldername, textfilename);
      if(1){cout <<folder_number<<"   "<<evt<<"       "<< infile <<"   "<< endl;}
      FILE *f = fopen(infile, "rb");
      //f.open(infile, "r");
      if (f == NULL) {
	printf("Cannot find file \'%s\'\n", infile);
	continue;
	//return 0;
      }
      fread(&th, sizeof(th), 1, f);
      fread(&ext, sizeof(ext), 1, f);  
      fread(wv_des, sizeof(int), 10, f);  
      fread(&uh, sizeof(uh), 1, f);
      fread(wv_array_count, sizeof(int), 9, f);  
      fread(words3, sizeof(words3), 1, f);  
      if(hv_no==0 && ch==0 && evt ==0){
	cout<<"SIZE th : "<<sizeof(th)<<endl;
	cout<<"SIZE uh : "<<sizeof(uh)<<endl;
	cout<<"SIZE bh : "<<sizeof(bh)<<endl;
      }
      fread(&vert_gain, sizeof(float), 1, f);
      fread(&vert_off, sizeof(float), 1, f);
      fread(&max, sizeof(float), 1, f);
      fread(&min, sizeof(float), 1, f);
      fread(words4, sizeof(words4), 1, f);
      fread(&h_interval, sizeof(float), 1, f);
      fread(&hor_off, sizeof(double), 1, f);
      fread(&pix_off, sizeof(double), 1, f);
      fread(&bh, sizeof(bh), 1, f);
      if(hv_no==0 && ch==0 && evt ==0){
	for(int ij=0;ij<11; ij++){printf("%c", th.unknown[ij]);}
	printf("\n");
	for(int ij=0;ij<16; ij++){printf("%c", th.blockname[ij]);}
	printf("\n");
	for (i=0 ; i<16 ; i++) {printf("%c", th.temp[i]);}
	printf("\n");
	printf("comtype : %X\n", th.comtype);
	printf("comorder : %X\n", th.comorder);
	for (i=0 ; i<16 ; i++) {printf("%c", uh.inst_name[i]);}
	printf("\n");
	printf("Length in bytes of descriptor : %i\n", wv_des[0]);
	printf("WS 1 : %i\n", wv_des[1]);
	printf("WS 2 : %i\n", wv_des[2]);
	printf("Length in bytes of data array : %i\n", wv_des[6]);
	for (i=0 ; i<9 ; i++) {printf("Number of data points : %i\n", wv_array_count[i]);}
	printf("vert_gain : %f\n", vert_gain);
	printf("vert_off : %f\n", vert_off);
	printf("max : %f\n", max);
	printf("min : %f\n", min);
	printf("h_interval : %g\n", h_interval);
	printf("hor_off : %g\n", hor_off);
	printf("pix_off : %g\n", pix_off); 
	printf("bandwidth : %d\n", bh.bandwidth); 
	printf("vert1_gain : %d\n", bh.vert1_gain);
	printf("wavesource : %d\n", bh.wavesource);
      }
      fread(voltage, sizeof(short), samples, f);
      for (i=0 ; i<samples ; i++) {
	//printf("%i  %d\n  ", i, voltage[i]);
	waveform[i] = ((voltage[i]*vert_gain)-vert_off)*tenp3;
	time[i] = ((h_interval*i) + hor_off)*tenp9;	
      }
      sprintf(name,"raw_%i", evt);			  
      TH1F* chnl_raw=new TH1F(name,name,wv_array_count[0],time[0],time[samples-1]); 
	      
      for (i=0 ; i<samples ; i++) {
	chnl_raw->SetBinContent(i+1, waveform[i]);
	//chnl_raw->Fill(time[i], waveform[i]);	
      }
      //cout<<evt<<endl;
      newfile->cd();
      if(chnl_raw){
      chnl_raw->Write(0,TObject::kOverwrite);
      delete chnl_raw;
      chnl_raw = 0;
      }
      
      evt++;
      if(evt>N_events-1){break;}
      fclose(f);
    }
    comfile_db.close();
    
    ch++;
    
    cout<<"Writing to output file"<<endl;
    for(int ij=0;ij<N_events;ij++){
      //  cout<<ij<<endl;
      // chnl_raw[ij]->Write(0,TObject::kOverwrite);
      // delete chnl_raw[ij];
      // chnl_raw[ij] = 0;
    }
    cout<<"check0"<<endl;

    if(newfile){
      newfile->Close();
      delete newfile; newfile=0;
    }
    cout<<"check1"<<endl;
  }
  cout<<"check2"<<endl;

  // delete[] voltage;
  // delete[] waveform;
  // delete[] time;
  cout<<"check3"<<endl;

  
}
	
