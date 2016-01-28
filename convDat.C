/*
 *
 * TO get properly TIME organized data from the file - 
 *  -  possible sort -n should be done on the data file first:
 *
 *  The buffer is read by chunk with chan0, then chunk with chan1
 *    time to sort ~ time to convert to root
 *===================================================
 * I have to rework:   1st  TIME
 *                     2nd  CHANNEL
 *  to properly sortout!!!
 *
 */

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

const double ROOT_offset=7.889112e+08;

const double COINWINDOW=3.0;


int convDat( const char* filename, int sortme=1 ){
  FILE *f;
  char line[500];
  int res;

  // basic structure for   "time ene sig chan"
  struct Tdata{
    ULong64_t time; // l
    //    UInt_t ene; // i
    Int_t ene; // i
    Int_t ch; // i        // b  1byte
    Int_t sign; // i      // b  1byte
    Int_t coE0; // i
    Int_t coE1; // i
    
  } event;

  
  //  ULong64_t lastt0=0;
  //  ULong64_t lastt1=0;
  //  ULong64_t lastt2=0;

  ULong64_t last_01_0=0;
  ULong64_t last_01_1=0;

  ULong64_t lastEvTime=0;


  int enelast0=0, enelast1=0;
  
  int last_flipflop=0; // channel 0  or  1

  ULong64_t evts=0;
  UInt_t d2;
  char outname[300];   // OUT =   infile.root
  char wc[100];int wci=0;
  sprintf(outname,"%sA.root", filename );

  
  TFile *fo=new TFile(  outname ,"recreate");
  TTree *to=new TTree("pha","dpp_pha_data");
  to->Branch("b1",&event.time,"time/l:ene/i:ch/i:sign/i:coE0/i:coE1/i");
   

  // HISTOGRAMS:
  TH1F *Hsaturation=new TH1F("saturation","Saturation times;time [us]",1000,-10000,100000);
  TH1F *erlang0=new TH1F("erlang0","erlang0;time [us]",100000,-10000,100000);
  TH1F *erlang1=new TH1F("erlang1","erlang1;time [us]",100000,-10000,100000);

  TH1F *coinc01=new TH1F("coinc01","Coincidences01;time [us]",100000,-50000,50000);
  TH2F *coinc01B=new TH2F("coinc01B","Coincidences01 bidim",500,0,5000,500,0,5000);

  TH1F *chan0=new TH1F("chan0","vme channel0",18000,0,18000);
  TH1F *chan1=new TH1F("chan1","vme channel1",18000,0,18000);

  TH1F *dchan0=new TH1F("dchan0","Saturation marked vme channel0 : 1 event after ",46000,0,46000);
  TH1F *dchan1=new TH1F("dchan1","Saturation marked vme channel1:  1 event after",46000,0,46000);
  //LATER!!  TH1F *timeline=new TH1F("timeline","timeline;time",100000,0,100000);


  
  //================== DIRTY LINES WIT BASH SCRIPTING =================
  // =========== to recover basic info on data file =================
  printf("counting lines.....\n%s", "");
  sprintf( line, "cat  %s | wc -l > LINES" , filename );
  system( line );

  printf("getting last line.....\n%s", "");
  sprintf( line, "cat  %s |  tail -1 | cut -d \" \" -f 1 > LASTTIME" , filename );
  system( line );

  printf("getting last line.....\n%s", "");
  sprintf( line, " stat -c \"%%Y\" %s > LASTMOD" , filename );
  system( line );
  // stat -c "%Y" run0001_20151113_165732.dat

  
  f=fopen( "LINES","r");
  fread( wc, 1, sizeof(wc), f); sscanf(wc,"%d",&wci); fclose(f);
  ULong64_t lasttime=0;
  f=fopen( "LASTTIME","r");
  fread( wc, 1, sizeof(wc), f); sscanf(wc,"%lld",&lasttime); fclose(f);
  f=fopen( "LASTMOD","r"); ULong64_t lmod;
  fread( wc, 1, sizeof(wc), f); sscanf(wc,"%lld",&lmod); fclose(f);

  ULong64_t starttime= lmod- lasttime/100./1e+6 ;
  
  printf("# LINES= %d....\n",         wci);
  printf("# lastLine =%lld x10ns\n",  lasttime);

  ULong64_t duration=lasttime/100/1000/1000;
  printf("# duration =%lld sec \n"  , duration);
  printf("# lastMOD = %lld sec.\n",  lmod);
  printf("# start   = %lld\n ", starttime);



 printf("checking creation date.....=lastMOD-duration\n%s", "");
 sprintf( line, "  date -d @%lld > LASTCREATE" , starttime );
 system( line );
 sprintf( line, " cat LASTCREATE; head -1  %s" , filename );
 system( line );

 char fsorted[300];
 if (sortme==1){
   sprintf(fsorted, "%s.sorted", filename );
   printf("SORTING TO NEW FILE.....\n%s", "");
   sprintf( line, "time cat %s | sort -n > %s" ,  filename , fsorted);
   system( line );
 }
 printf("COINCIDENCE WINDOW IS.....%.1f us\n", COINWINDOW);

 //=================================================== BASH END
 //==========================================================================
 
 TH1F *timeline0=new TH1F("timeline00","timeline00;time",
			    (lmod-starttime) ,
			    1.0*starttime - ROOT_offset,
			    1.0*lmod- ROOT_offset  );
 TH1F *timeline1=new TH1F("timeline01","timeline01;time",
			    (lmod-starttime) ,
			    1.0*starttime - ROOT_offset,
			    1.0*lmod- ROOT_offset  );


 int  sum0=0, sumnega=0, sum1637=0, sumpileup=0,sumtot=0, sumSatu=0, sumRoll=0, sumTTres=0,sumFake=0;
 int FORMAT=0; //  caen = 0   HEADER;     gregory= 1    


 printf("==================================================\n%s","");
 int  test, test2, test3, chan,  energy;
  //  f=fopen( filename ,"r");
 if (sortme==1){
   f=fopen( fsorted ,"r");
 }else{
   f=fopen( filename ,"r");
 }   
  int linenow=0;
  if (f!=NULL){
    //    while (not(eof(f) )){
    while ( fgets( line, 300, f ) != NULL ){
      //      res=fread( line, 1,  sprintf( title, "read from %s", fname );

      
      while ( (line[0]=='#')||(line[0]=='H') ){
	if (line[0]=='#'){ printf("changed to Gregory format \n",""); FORMAT=1;}
	if (line[0]=='H'){ printf("changed to Caen mc2 format\n%s","");FORMAT=0;}
	fgets( line, 300, f );
	// printf("%s\n", "comment line # or H");
      }

      
      // DT5780 :   not ch but pilup
      //ORIGO      sscanf( line , "%lld %d %d %d", &event.time, &event.ene, &test, &chan );  
      // DT5780 -------->>>>>>>   3 columns.  Energy sometimes < 0

      //      FORMAT=2;


      if (FORMAT==0)//  CAEN
	sscanf( line , "%lld %d %d %d", &event.time, &energy, &test );
      if (FORMAT==1)// GREGORY	
	sscanf( line , "%lld %d %d %d", &event.time, &energy, &test, &chan );  

      if (FORMAT==2)
	sscanf( line , "%lld %d %d %d %d %d", &event.time, &energy, &test, &chan,
		&test2, &test3);  
      sumtot++;

      
      if (energy>16370)sum1637++;
      if (energy<0)sumnega++;
      if (energy==0)sumpileup++;

      if (test==0)sum0++;
      if ( (test&0b1)!=0){
	sumSatu++; //Saturation before
	Hsaturation->Fill(  0.01*(event.time - lastEvTime) );
      }
      lastEvTime=event.time;
      if ( (test&0b10)!=0) sumRoll++; //RolloverTime before
      if ( (test&0b100)!=0) sumTTres++; //externat TTreset
      if ( (test&0b1000)!=0) sumFake++; //fake time sync

      //      sscanf( line , "%lld %d %d %d", &event.time, &test2, &event.ene, &test );   //test2=ch; run0002
      //       event.ch=2;
      if (linenow<=5){
       printf("%14lld %5d  /%d/  chan=%d \n",event.time,event.ene,test,chan );
       linenow++;
      }

      event.sign=test;
      event.ch=chan;
      event.ene=energy;
	event.coE0=0;
	event.coE1=0; //i dont know if it will be coincidence or single


      //===================  chan 0  or  chan 1   sections
      if (event.ch==0){
	event.coE0=0;
	event.coE1=0;
	timeline0->Fill(  1.0*event.time/100./1.e+6 + 1.0*starttime- ROOT_offset );
	erlang0->Fill(  0.01*(event.time-last_01_0) ); // 0.01...in us
	chan0->Fill(event.ene);
	if (( (event.sign&0b1)!=0)||(event.ene<=0) ){dchan0->Fill(event.ene);}
	
	coinc01->Fill( 0.01*(event.time-last_01_1) );

	if ( abs( 0.01*(event.time-last_01_1))<COINWINDOW )
	  {
	    coinc01B->Fill( enelast1, event.ene );
	    event.coE0=event.ene;
	    event.coE1=enelast1;
	  }
 	last_01_0=event.time;
	enelast0=event.ene;
     }
      if (event.ch==1){
	event.coE0=0;
	event.coE1=0;
	timeline1->Fill(  1.0*event.time/100./1.e+6 + 1.0*starttime- ROOT_offset );
	erlang1->Fill( -0.01*(event.time-last_01_1) ); // 0.01...in us
	chan1->Fill(event.ene);
	if (( (event.sign&0b1)!=0)||(event.ene<=0)){dchan1->Fill(event.ene);}
	
	coinc01->Fill( 0.01*(event.time-last_01_0)  );
	if ( abs( 0.01*(event.time-last_01_0))<COINWINDOW ){
	  coinc01B->Fill( event.ene, enelast0 );
	  event.coE0=enelast0;
	  event.coE1=event.ene;
	}
	last_01_1=event.time;
	enelast1=event.ene;
      }

      //===========FILL EVENT
      if ( (test&0b1000)==0){  //not fake
	to->Fill();
      }

      evts++;
      //      if (event.ch<2){last_01=event.time;} // either 0 or 1 ... last time for coin
      
      if (evts%50000==0){
	printf("\r   %9.3f / %9.3f M events processed",1.0e-6*evts, 1e-6*wci  );
	fflush(stdout);    
	//to->Write();    fo->Close();	break;
      } // printout
    }//  whil not eof
    fclose(f);
    printf("\n%s","");
    to->Print();
    to->Write();
    erlang0->Write();
    erlang1->Write();
    coinc01->Write();
    coinc01B->Write();
    chan0->Write();
    chan1->Write();
    dchan0->Write();
    dchan1->Write();
    timeline0->Write();
    timeline1->Write();
    Hsaturation->Write();
    
    fo->Close();
    printf( "TFile *ff=new TFile(\"%s\"); \n", outname );
  }else{
    printf("cannot open <%s>\n",  filename );
  } // f NOT NULL


  if (FORMAT==1){ printf("Gregory format was used\n",""); }
  if (FORMAT==0){ printf("Caen mc2 format was used\n%s","");}
  if (FORMAT==2){ printf("Special format was used - check source code\n%s","");}

  printf("# duration = %lld s \n",lasttime/100/1e+6 );
  printf("TOTAL brut.. %6d  \n",  sumtot );
  printf("TOTAL-FAKE.. %6d  \n",  sumtot-sumFake );

  double countrate=1.0*(sumtot-sumFake);
  //   ULong64_t duration=lasttime/100/1.0e+6;
   //  printf("COUNTRATE .. %6.1f  cps\n",  countrate );
   //  printf("COUNTRATE .. %6.1f  cps\n",  duration );
  countrate=countrate/duration;
  printf("COUNTRATE .. %6.1f  cps\n",  countrate );



  printf("E>16370  ... %6d\n",  sum1637);
  printf("E==0     ... %6d  ... PileUP   %.2f %% DT\n",sumpileup,100.0*sumpileup/(sumtot-sumFake));
  printf("E<0      ... %6d  ... after <0, PileUP is frequent\n\n",  sumnega);


  printf("xtra==0  ... %6d\n",  sum0);
  printf("Saturation.. %6d       %.2f++ %% DT\n",   sumSatu ,  100.0*sumSatu/sumtot );
  printf("Rollover  .. %6d      \n",   sumRoll  );
  printf("extTTreset.. %6d   \n",   sumTTres );
  printf("Fake     ... %6d\n"  ,   sumFake );




  
}// ro_pha2tree
