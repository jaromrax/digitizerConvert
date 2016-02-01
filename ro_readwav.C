#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TPad.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TList.h"
//==============  rotate
#include<algorithm>
//===============
#include "TTree.h"
#include "TFile.h"

#include "coinonline.h"


#include "TVirtualFFT.h"
/*

scp aaron:/home/ojr/root/lib/root/libFFT*  /home/ojr/root/lib/root/


tree->BuildIndex("SNR");
TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();
for( int i = index->GetN() - 1; i >=0 ; --i ) {
    Long64_t local = tree->LoadTree( index->GetIndex()[i] );
    branch->GetEntry(local);
    ....
}
*/
//
//  #include <stdint.h>Short_t
// READ WAVEFORM   FOLLOWING 1st int16  size
//=============================
//  mixed mode: 47471
//    81582 generators in chan0, total   163164 in ch0,1
//    total 47471 samples ... ~23735 per channel  3.4x less
//
// 0 0  7   0 200 61447 0 0 41431  8        1 
// b ch ev  fmt    ?              
// 0 0  6   0 200 61447 0 0 4506   65524    0      nxt   e extr???
//                 ?    0 0 fast   slower   upbit
// 1 2  3   4 5    6    7 8  9      10      11     12    13 14                                                 
// ######################################
//##   #  fitting the  line-slope-line
//#######################################
  Double_t fitf1(Double_t *v, Double_t *par) {
    //     Double_t low=par[0];
    Double_t low=0.0;
    Double_t crs1=par[0];
    Double_t crs2=par[1];
    Double_t hig=par[2];
    double penalty=1.0;
    double fitval;

    
     if (hig<=low){ penalty=1.0+0.1*(low-hig);hig=low;}
     if (crs2<=crs1+2) {penalty=1.0+10; }
     if (v[0] < crs1 ){       fitval=low;     }
     if ( (v[0] >= crs1)&&(v[0]<= crs2) ){
       fitval= (hig-low)/(crs2-crs1) * (v[0]-crs1) +  low;
     }
     if (v[0] >crs2 ){       fitval=hig;     }
     return fitval* penalty;
   } //--------------------------DEFINITION END --------------

//########################################
//########     FIT - WRAP FUNCTION
//########################################
double* fit1( int max, TH1F *k1, char *OPT){
 
  float avgend=200.0;
  int csini1,csini2,csini3;
  int i;
  int newmax=max;
  // prepare initial value for amplitude by averaging:
  //  for (i=250;i<250+105;i++){avgend=avgend+k1->GetBinContent(i);}
  //  avgend=avgend/(105.);
  avgend=k1->GetMaximum();
  
  //now look for 1st avgend/2 value
  csini1=50; csini2=105;
  for (i=0;i<max;i++){
    if (k1->GetBinContent(i) < avgend/5.){ csini1=i;}
    if (k1->GetBinContent(i) > avgend/2.){ csini2=i;break;}
  }
  
  newmax=csini2+ (csini2-csini1)+100;
  //  printf("new maximum =%d\n", newmax );
  TF1 *func1 = new TF1("fit1", fitf1, 0, 1.0*newmax,  3); //3 params
  //k1->GetBinContent(0),  k1->GetBinContent(max) 
  func1->SetParameters(   csini1, csini2+ (csini2-csini1) , avgend);
  //  func1->FixParameter( 2, avgend );
  //  func1->SetParameters(  avgend , 80,120 );
  func1->SetParNames(    "crs1","crs2","hig");
  
  k1->Fit("fit1", OPT ); // .... WW error=1; N nodraw; Q quiet
      //      k1->Fit("fit1", "WWN0QR" ); // .... WW error=1; N nodraw; Q quiet
      //      k1->Fit("fit1", "WWN0Q"); // .... WW error=1; N nodraw; Q quiet
  //   k1->Fit("fit1", "WW"); // 
  // WWL no
  double rlow=0;
  //func1->GetParameter(0);
  double rhig=func1->GetParameter(2);
  rhig=k1->GetMaximum()-k1->GetMinimum();
  double rc1=func1->GetParameter(0);
  double rc2=func1->GetParameter(1);
  double xi=func1->GetChisquare();

  double* d=new double[3];
  d[0]=rhig-rlow;
  d[1]=rc2-rc1;
  d[2]=xi;
   if (gPad!=NULL){
         k1->SetStats(kFALSE);  gPad->Modified(); gPad->Update();
   }
  return d;
}//........................................








/**********************************************************
 * #######################################################
 *        MAIN PROGRAM 
 *
 *
 *
 *
 *   start stop ... important for parallel!
 *
 *
 *   clerun=  "run"      run number unknown
 *        new filename can be local OR remote
 *        important for parallel run
 *
 *    options WW ... weight==1 ! important
 *             0Q ... quiet, ...
 *
 *
 *ro_readwav("run0050_20151212_090602.dat",0,9999999,"WW")
 *
 *   without fit=  361770/2 sec.     1568 seconds
 ***********************************************************/
void ro_readwav(char* fname , int start=0, int stop=0, char *OPT="WWN0Q", char *clerun="run"){

  FILE *f;
  int i,j=0;
  Short_t ssample[10];
  UShort_t buffer[2200]; // usually 600

  //============ OPEN FILE========== READ 1st word
  f=fopen( fname , "r" );
  fread( ssample, sizeof(Short_t), 1 , f );
  printf("i ... sample size = %d ... one integer is %d  bytes\n",
	 ssample[0] , sizeof(Short_t)  );
  // read the rest =================== 1st sample 
  fread( buffer, sizeof(Short_t), ssample[0]-1, f );
  for (i=0;i<ssample[0];i++){buffer[i]=0;} // clean buffer
  printf("%s\n", "buffer clean");


  //----------------- DECLARATIONS #2
  TH1F *hch[8]; // channel histograms: 0-7 FW
  for (i=0;i<8;i++){hch[i]=NULL;}
  TH1F *hchf[8]; // chan fit 
  for (i=0;i<8;i++){hchf[i]=NULL;}
  TH1F *htime[8]; // I cannot Fill time histo !?
  for (i=0;i<8;i++){htime[i]=NULL;}
  TH2F* bidim[8];  // slope vs ene
  for (i=0;i<8;i++){bidim[i]=NULL;}
  //  TH2F* ede=new TH2F("ede","E x dE matrix",400,0,8000,400,0,8000);



  
  TH1F *h; // WAVE HISTOGRAM
  
  char hname[20];   // histo channels
  int result;    
  int sizeOfStr=16; // SIZE OF struct with b,ch,t,E
  // ---- structure data:
  Short_t b,ch, evnt, ene ;  // board channel energy
  unsigned long long time=0;
  // tets differences
  unsigned long long timelast=0;
  Short_t extr;
  float avg=0.0;
  // fit parameters:
  double* pars;
  

    //==============================FILE  TTREE ================
    struct Tdata{
      ULong64_t t;   //l  time 
      Int_t n;       // i  event number
      Int_t ch;      // i   channel
      Int_t e;       //i    energy
      Float_t efit;  //F    energy from fit
      Float_t xi;    //F    Xi2 from  fit
      Float_t s;     //F    s - slope from fit
                     //-------------------COINC from NOW
      Int_t e0;      // COINC CASE    E0
      Int_t e1;      // COINC CASE    E0
      Int_t dt;      //I   +-

      Int_t efit0;
      Int_t efit1;
      Int_t s0;
      Int_t s1;
      
    } event;
      
    char outname[300];   // OUT =   infile.root
    char wc[100];int wci=0;
    //    sprintf(outname,"%s_%06d_%06d.root", fname , start, stop );
    sprintf(outname,"%s_%06d_%06d.root" , clerun, start, stop );
    
    TFile *fo=new TFile(  outname ,"recreate");
    TTree *to=new TTree("pha","dpp_pha_data");
    to->Branch("time",&event.t,"t/l");
    to->Branch("data",&event.n,"n/i:ch/i:e/i:efit/F:xi/F:s/F");
    to->Branch("coin",&event.e0,"e0/i:e1/i:dt/I:efit0/I:efit1/I:s0/I:s1/I");


    MyCoinc *mc=new MyCoinc(0,1);
//
    //
    //

    //#####################################  LOOP #####################
    int JMIN=start;
    int JMAX=stop;
    if (stop==0){JMAX=999999999;}
    //  JMAX=111;    // limit for test
    for (j=0;j<=JMAX;j++){
    
    event.n=j;
    event.e=0;
    event.efit=0;
    event.xi=0.0;
    event.s=0.0;


    // ch  t 
    
    result=fread( buffer,sizeof(Short_t) , sizeOfStr, f );
    if (result<sizeOfStr) break;
    //---------------- FW header decode
    b=buffer[0];
    ch=buffer[1];
    evnt=buffer[2];
    //[8]-[11] ... time
    //--------------- Time decode
    time=0;
    time=buffer[11];      time=time<<16;
    time=time+buffer[10]; time=time<<16;
    time=time+buffer[9]; time=time<<16;
    time=time+buffer[8];
    timelast=time;
    
    ene=buffer[12];
    extr=buffer[13];
    //    printf( "ch/evnt/ene   %3d   %3d  %3d  E=%3d  %3d\n", ch, evnt, ene, extr );
    //    printf("%lld  =======  %3d  \n", time, ene);

    event.t=time;
    event.e=ene;
    event.ch=ch;

    
    //================ FILL ENERGY HISTO
    if (hch[ ch ]==NULL){
      sprintf(hname,"hifw_%02d", ch);
      hch[ ch ]= new TH1F( hname,hname,16000,0,16000);
    }
    //  printf("%s\n", "J loop hifw ok");
    hch[ch]->Fill( ene );

    //================ FILL TIME HISTO
    if (htime[ch]==NULL){
      sprintf(hname,"hitime_%02d", ch); 
      //      htime[ch]= new TH1F( hname,"timeprofile;[seconds]",1000000,0,100);
      int seconds=1000;
      htime[ch]= new TH1F( hname,"hname",seconds*100,0,seconds);
      printf("defined  %d. htime\n",ch);
    }
    htime[ch]->Fill( float(time/100000000.) );

    if (j%1000==0){printf( "j=%d ...\n", j );}

    
    //    if ( j>= start){ //=======PARALLEL LOOP===============
    //===============================================LOAD WAVEFORM
    result = fread( buffer,sizeof(Short_t) , ssample[0], f );
    if (result< ssample[0]) break;

    
   


    
    // ============ FIT / ANALYZE ================
    int makefit=0;
    if ((j>=start)&&(j<stop) ){ makefit=1;}
    
    if (makefit==1){
      //====CREATE HISTOGRAM wave ==>>> VARIANT:create/analyze/root it/delete
      sprintf(hname,"wave_%06d", j); // here I keep the channel
      h=new TH1F( hname,hname, ssample[0],0,ssample[0] );
      //====== prepare ZERO (one analyze parameter now)
      avg=0.0;
      for (i=0;i<50;i++){avg=avg+buffer[i];}
      avg=avg/(50.);


      // ### ======= Fill SHIFTED TO ZERO histogram
      for (i=0;i<ssample[0];i++){    // i add ch 0/1 to the end!
	h->SetBinContent( i+1, buffer[i] - avg ); // 0-1 is #1
      }

      // ### FFT HERE //
      int n=ssample[0];
      if (1 == 1){
	TH1 *hmsm =0;
	TH1 *hpsm =0;
	hmsm=h->FFT(hmsm,"MAG");
	hpsm=h->FFT(hpsm,"PH");
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	Double_t *re_full = new Double_t[ssample[0]];
	Double_t *im_full = new Double_t[ssample[0]];
	fft->GetPointsComplex(re_full,im_full);
	for (int j=0;j<ssample[0]/2; j++){
	  hmsm->SetBinContent( j,  re_full[j]   );
	  hpsm->SetBinContent( j,  im_full[j]   );
	}
	//	for (int j=16;j>11;j--){
	hmsm->Smooth(2);
	hpsm->Smooth(2);
	for (int j=16;j>5;j--){
	  hmsm->SetBinContent( j,  0.8*re_full[j]   );
	  hpsm->SetBinContent( j,  0.8*im_full[j]   );
	}
	hmsm->Smooth(2);
	hpsm->Smooth(2);

	for (int j=0;j<ssample[0]/2;j++){
	  re_full[j]=hmsm->GetBinContent(j);
	  im_full[j]=hpsm->GetBinContent(j);
	}
	TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
	fft_back->SetPointsComplex(re_full,im_full);
	fft_back->Transform();
	TH1 *hb = 0;
	hb = TH1::TransformHisto(fft_back,hb,"Re");
	double mx= mx=hb->GetMaximum();
	double omax= h->GetMaximum();
	hb->Scale( omax/mx );
	delete fft_back;
	delete [] re_full;
	delete [] im_full;
	if (gPad!=NULL){
	  h->Draw();
	  hb->Draw("same");hb->SetLineColor(3);
	  gPad->Modified(); gPad->Update();
	}
	for (int j=0;j<ssample[0];j++){
	  h->SetBinContent( j, hb->GetBinContent(j) );
	}
	delete hb;
	delete hmsm;
	delete hpsm;
	//	sleep(1);
	//	h->Draw(); gPad->Modified(); gPad->Update();

      }// ========== FFT HERE  END
      
      // ### ===========     fit 1 
      pars=fit1( ssample[0],  h , OPT ); // FIT 400 points in a sample
      printf("%06d/ %6.1f %5.2f    %5.2f  \n", j, pars[0], pars[1] ,pars[2]  );
      if ( bidim[ch]==NULL ){
	sprintf(hname,"bishape_%02d", ch);
	bidim[ch]=new TH2F( hname,"shapes;Energy;slope in channels", 2500,0,2500,240,2,25);
      }
      bidim[ch]->Fill( pars[0], pars[1] );
      //      sleep(1);
      
      //##### DEBUG HERE".....................
      // # pars[1]===S
      // # pars[0]==E
      //      if ( (ch==0)&&(pars[0]>960.) ){
      /*
      if ( (ch==0)&&(pars[1]>100.) ){
	printf("S problem @ %d chan %d S=%f\n",j,ch, pars[1]);
	sleep(5);
      }
      */
      //================ FILL ENERGY FIT HISTO
      if (hchf[ ch ]==NULL){
	sprintf(hname,"hifit_%02d", ch);
	hchf[ ch ]= new TH1F( hname,hname,16000,0,16000);
      }
      hchf[ch]->Fill( pars[0] );

      event.efit=pars[0];
      event.s=pars[1];
      event.xi=pars[2];
      if ( strcmp(OPT,"WWN0Q")==0){ delete h;}else{h->Write();} //write to tree if not N0Q
    } // fit
    // =========== FIT / ANALYZE  END  ================


    event.e0=0;
    event.e1=0;
    event.dt=0;
    
    //=================COINCIDENCES ==========push
    if (1==mc->add( ch, event.t, event.e, event.efit, event.xi,event.s) ){
                                    // event.efit
      event.e0=mc->getlastC0();
      event.e1=mc->getlastC1();
      event.efit0=mc->getlastC0fit();
      event.efit1=mc->getlastC1fit();
      event.s0=mc->getlastC0s();
      event.s1=mc->getlastC1s();
      event.dt=-mc->getlastDT();
      if ( (event.dt>400)||(event.dt<-400)){
	printf("%4d mc ============================\n", event.dt);
      }
    }
    //=================COINCIDENCES END ==========
    


    
    if (makefit==1){ delete pars; } // no pars if no fit
    //    if (h!=NULL){delete h;}   // delete  WAVE HISTOGRAM

    //==========SAVE TTREE======
    //    if (makefit==1){
      to->Fill();
      //    }
    //    } // ===========PARALLEL LOOP ===========
    
  } //------- j < JMAX ------=====================================END




    fclose(f);
    
    if (j>=JMAX){  // i dont create thousands of th1f anymore
      printf("!... limit in number of histograms reached....%d.  STOP \n", j);
    }else{
      printf("!... total number of histograms  == %d \n", j);
      }     
    if ( (gPad!=NULL)&&(bidim[0]!=NULL)){    bidim[0]->Draw("col"); }
    double dlt=(timelast/100/1000);
    printf("%7.3f s =last time\n", dlt/1000.);

    //====SAVE==============================
    for (ch=0;ch<8;ch++){
      //      printf("ch == %d\n", ch);
      if (hch[ch]!=NULL){ hch[ch]->Write();}
      if (htime[ch]!=NULL){ htime[ch]->Write();}
      if (bidim[ch]!=NULL){ bidim[ch]->Write();}
    }

    //  ede->Write();
    // ttree:=====================
    //    mc->print();
    TH2F* mcbi=(TH2F*) mc->getbidim();
    TH1F* mchidt=(TH1F*) mc->gethidt();
    mcbi->Write();
    mchidt->Write();
 
    to->Print();
    to->Write();
    fo->Close();

} //######################################################## END ##############
