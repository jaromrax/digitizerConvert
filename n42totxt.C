// simplexml -v RadInstrumentData/RadMeasurement/Spectrum/ChannelData  20160121-bg001.n42  > z.dat
//
//simplexml -v RadInstrumentData/RadMeasurement/Spectrum/ChannelData  20160121-bg001.n42 | sed ':a;N;$!ba;s/\n/ /g'  |  tr " " "\n" > z.dat


int getvalue(const char *ch, const char *ext){
  int i=0;
  char ba[500],name[500];
  sprintf(name, "%s%s", ch, ext );
  FILE* temf=fopen( name ,"r");
  if (temf!=NULL){  fread( ba, 1, sizeof(ba), temf ); fclose(temf);}
  i=atoi( ba);
  return i;
}


float getreal(const char *ch, const char *ext){
  float i=0;
  char ba[500],name[500];
  sprintf(name, "%s%s", ch, ext );
  FILE* temf=fopen( name ,"r");
  if (temf!=NULL){  fread( ba, 1, sizeof(ba), temf ); fclose(temf);}
  i=atof( ba);
  return i;
}

const char* getchar(const char *ch, const char *ext){
  char name[500];
  char *ba=new char[100];
  sprintf(name, "%s%s", ch, ext );
  FILE* temf=fopen( name ,"r");
  if (temf!=NULL){  fread( ba, 1, 100, temf ); fclose(temf);}
  printf("BA: %s\n", ba);
  return ba;
}


//============================================================================
//============================================================================
//============================================================================

void n42totxt(const char* fname){
  TH1F* h=new TH1F("n42","title", 16000,0,16000);
  TString oneline, title=fname, token;
  int MAXLINES=64000;
  int i,j;

//======================================Conversion to cmd
  char cmd[500], newname[500], basename[500], newtitle[500];
  //  sprintf(cmd, "simplexml -v RadInstrumentData/RadMeasurement/Spectrum/ChannelData  20160121-bg001.n42 | sed ':a;N;$!ba;s/\n/ /g'  |  tr \" \" \"\n\" > z.dat", fname);
  sprintf( cmd, "basename -z  %s .n42 | tr \"\\n\" \" \" > BASENAME", fname );
  //  printf("%s\n", cmd );
  system(cmd);
  FILE* temf=fopen("BASENAME","r");
  if (temf!=NULL){  fread( basename, 1, sizeof(basename), temf ); fclose(temf);}

  // calib
  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/EnergyCalibration/CoefficientValues %s | cut -d \" \" -f 2 > %s.cal1 ", fname, basename );
  //  printf( "%s\n", cmd );
  system(cmd);
  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/EnergyCalibration/CoefficientValues %s | cut -d \" \" -f 1 > %s.cal0 ", fname, basename );
  //  printf( "%s\n", cmd );
  system(cmd);


  // realtime
  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/RealTimeDuration %s | cut -b 8-9 > %s.durD ", fname, basename );
  system(cmd);

  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/RealTimeDuration %s | cut -b 12-13 > %s.durH ", fname, basename );
  system(cmd);

  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/RealTimeDuration %s | cut -b 15-16 > %s.durM ", fname, basename );
  system(cmd);

  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/RealTimeDuration %s | cut -b 18-19 > %s.durS ", fname, basename );
  system(cmd);
  int totalsec=  getvalue(basename,".durD")*24*3600+
  getvalue(basename,".durH")*3600 +
  getvalue(basename,".durM")*60 +
  getvalue(basename,".durS");
  printf("%d seconds REAL\n",totalsec);


    // livetime
  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/Spectrum/LiveTimeDuration %s | cut -b 8-9 > %s.livD ", fname, basename );
  system(cmd);

  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/Spectrum/LiveTimeDuration %s | cut -b 12-13 > %s.livH ", fname, basename );
  system(cmd);

  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/Spectrum/LiveTimeDuration %s | cut -b 15-16 > %s.livM ", fname, basename );
  system(cmd);

  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/Spectrum/LiveTimeDuration %s | cut -b 18-19 > %s.livS ", fname, basename );
  system(cmd);
  int livesec=  getvalue(basename,".livD")*24*3600+
  getvalue(basename,".livH")*3600 +
  getvalue(basename,".livM")*60 +
  getvalue(basename,".livS");
  printf("%d seconds LIVE      .....  %.2f %% deadtime\n",livesec,  100.0*(totalsec-livesec)/totalsec );

  


  // START TIME
  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/StartDateTime %s > %s.start ", fname, basename );
  system(cmd);



     // DATA ITSELF ===================================================
  sprintf( newname, "%s.dat", basename );
  sprintf( newtitle, "h%s", basename );
  h->SetName( newtitle );
  sprintf( newtitle, "%s  %.2f%% DT  cal: %.4f,%.4f", getchar(basename,".start") ,  100.0*(totalsec-livesec)/totalsec,   getreal(basename,".cal1"),  getreal(basename,".cal0")  );
  h->SetTitle( newtitle );

  sprintf(cmd, "/usr/bin/simplexml -v RadInstrumentData/RadMeasurement/Spectrum/ChannelData %s | tr \"\\n\" \" \"   |  tr \" \" \"\\n\" > %s ", fname, newname );
  //  printf( "%s\n", cmd );
  system(cmd);
  
  FILE * pFile;
  pFile=fopen( newname  ,"r" ); 
  if (pFile==NULL) {
    printf("cannot open %s,STOPping\n", fname ); 
    return 0;
  } // error
  
    i=0;
  int lastlen;// remove spaces
  while ((i<MAXLINES)&&( feof(pFile)==0) ){ 
    if ( oneline.Gets (pFile, kTRUE) ){//chop true ---------
      /*
      do {
	lastlen=oneline.Length();
          if (oneline.Index(" ")==0){oneline=oneline(1,oneline.Length()-1);}
      lastlen=oneline.Length();
      if (oneline.Index(" ")==lastlen){oneline=oneline(0,oneline.Length()-1);}

       oneline.ReplaceAll("\t"," ");
       oneline.ReplaceAll("  "," ");
       oneline.ReplaceAll("  "," ");
      }while( lastlen!=oneline.Length());
      */      
      if (i<10)printf("%d  %s\n", i, oneline.Data() );
      for (int j=0;j<oneline.Atoi();j++){
	//	h->SetBinContent( i ,  oneline.Atof()  );
	h->Fill(  i  );
      }
      i++;
    }// if gets ------------------------------------------
    
  }// while end

  h->Draw();
  h->Print();
  printf( "CALIBRATION:\n%.4f,%.4f\n", getreal(basename,".cal1"),  getreal(basename,".cal0")  );
  sprintf( newname, "%s_histo.root", basename );
  TFile* f=new TFile(newname,"update");
  h->Write();
  f->Close();
} // FUNCTION =============================================
// simplexml -v RadInstrumentData/RadMeasurement/Spectrum/ChannelData  20160121-bg001.n42  > z.dat
/*

 simplexml -v RadInstrumentData/RadMeasurement/Spectrum/ChannelData  20160121-bg001.n42 | sed ':a;N;$!ba;s/\n/ /g'


*/
