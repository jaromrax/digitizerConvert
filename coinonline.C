#include "coinonline.h"

#include<stdio.h>
// ClassImp(MyCoin);



OneChan::OneChan(int ch){
  printf("created OneChan %d\n",ch);
  chan=ch;
  max=0;
  for (int i=0;i<BMAX;i++){time[i]=0;value[i]=0;efit[i]=0;xi[i]=0;s[i]=0;}
}



OneChan::~OneChan(){
  printf("destroing OneChan %d\n",chan);
}


int OneChan::push(unsigned long long int t,int val,int vefit,int vxi,int vs){
  std::rotate(time,  time+BSIZE-1, time+BSIZE);
  std::rotate(value, value+BSIZE-1, value+BSIZE);
  std::rotate(efit, efit+BSIZE-1, efit+BSIZE);
  std::rotate(xi, xi+BSIZE-1, xi+BSIZE);
  std::rotate(s, s+BSIZE-1, s+BSIZE);
  
  time[0]=t;
  value[0]=val;
  efit[0]=vefit;
  xi[0]=vxi;
  s[0]=vs;
  if (max<BSIZE){max++;}
}

int OneChan::print(){
  char a[100];
  sprintf(a,"%*s",chan*15," ");
  for (int i=0;i<max;i++){
    printf("%s %6d %9lld \n", a, i, time[i]);
  }

} // print



int OneChan::getmax(){
  return max;
}

unsigned long long int OneChan::gett(int i){
  return time[i];
}


int OneChan::getv(int i){
  return value[i];
}

int OneChan::getvf(int i){
  return efit[i];
}
int OneChan::getvs(int i){
  return s[i];
}


int OneChan::remove(int i){
   std::rotate(time+i,  time+i+1, time+max);
   std::rotate(value+i, value+i+1, value+max);

   std::rotate(efit+i, efit+i+1, efit+max);
   std::rotate(xi+i, xi+i+1, xi+max);
   std::rotate(s+i, s+i+1, s+max);
   max--;
}

//===================================== COINC ===========
//===================================== COINC ===========

MyCoinc::MyCoinc(int ch0,int ch1){
  max=0;
  totalsum=0;
  printf("creating MyCoinc\n%s","");
  OC0=new OneChan(ch0);
  OC1=new OneChan(ch1);
  char ss[100];
  sprintf(ss,"coin_%d_%d", ch0,ch1);
  bidim=new TH2F(ss,ss,500,0,5000,500,0,5000);  // the 1st BIDIM = coincidences
  sprintf(ss,"coindt_%d_%d", ch0,ch1);
  hidt=new TH1F(ss,ss,WIN*2+1,-WIN,WIN);

  sprintf(ss,"fitcoin_%d_%d", ch0,ch1);
  fitbidim=new TH2F(ss,ss,500,0,2000,500,0,2000);  // the 2st BIDIM = coincidences

}


MyCoinc::~MyCoinc(){
  printf("destroing MyCoinc\n%s","");
  delete OC0;
  delete OC1;
}

Int_t MyCoinc::getlastC0(){return ch0last;}
Int_t MyCoinc::getlastC1(){return ch1last;}
Int_t MyCoinc::getlastDT(){return dtlast;}


Int_t MyCoinc::getlastC0fit(){return  ch0flast;}
Int_t MyCoinc::getlastC1fit(){return  ch1flast;}
Int_t MyCoinc::getlastC0s(){  return  ch0slast;}
Int_t MyCoinc::getlastC1s(){  return  ch1slast;}



int MyCoinc::add(int ch,unsigned long long int t,int val){
  return add(ch,t,val,   0,0,0);
}

int MyCoinc::add(int ch,unsigned long long int t,int val,
		 int vefit,int vxi,int vs){
  int i, sig;
  int found=0;
  unsigned long long int t1,t0, dt;
  ch1last=0;ch0last=0;dtlast=0;
  ch1flast=0;ch1slast=0;
  ch0flast=0;ch0slast=0;

  if (ch==0){
    // search in ch==1 first
    for (i=0;i<OC1->getmax();i++){
      t1=OC1->gett(i);
      if (t1>t){ dt=t1-t;sig=1;}else{dt=t-t1;sig=-1;}
      if ( dt<WIN){
	if (DEBUG>=1){printf("%lld ~ %lld from ch=%d\n", t, t1, 1);}
	bidim->Fill( OC1->getv(i) ,val); totalsum++;
	fitbidim->Fill( OC1->getvf(i) , vefit);// totalsum++;
	ch1last=OC1->getv(i);ch0last=val;
	ch1flast=OC1->getvf(i);ch0flast=vefit;
	ch1slast=OC1->getvs(i);ch0slast=vs;

	dtlast=dt*sig;
	hidt->Fill(dtlast);
	found++;
	OC1->remove(i);
	break;
      } // if
    } //all i
    if (found ==0){ OC0->push(t, val, vefit,vxi,vs);}
  } //ch==0
  
  if (ch==1){
    // search in ch==1 first
    for (i=0;i<OC0->getmax();i++){
      t0=OC0->gett(i);
      if (t0>t){ dt=t0-t;sig=-1;}else{dt=t-t0;sig=1;};
      if ( dt<WIN){
	if (DEBUG>=1){printf("%lld ~ %lld from ch=%d\n", t, t0, 1);}
	bidim->Fill(val, OC0->getv(i));totalsum++;
	fitbidim->Fill(vefit, OC0->getvf(i)); // fit
	ch1last=val;ch0last=OC0->getv(i);
	ch1flast=vefit;ch0flast=OC0->getvf(i);
	ch1slast=vs;ch0slast=OC0->getvs(i);

	
	dtlast=dt*sig;
	hidt->Fill(dtlast);
	//printf("%4d  co\n", dtlast);

	found++;
	OC0->remove(i);
	break;
      } // if
    } //all i
    if (found ==0){ OC1->push(t, val, vefit,vxi,vs);}   
  } //ch ==1 

  if (found!=0){ return 1;}
  if (DEBUG>=2){
  print();
  }
} //======== add =====





int MyCoinc::print(){
  printf("---------------------%s\n","");
  OC0->print();
  OC1->print();
  printf("Total coincidences: %d\n",totalsum );
}



TH2F* MyCoinc::getbidim(){
    return bidim;
}
TH2F* MyCoinc::getfitbidim(){
    return fitbidim;
}
TH1F* MyCoinc::gethidt(){
    return hidt;
}


int coinonline(){
  
  MyCoinc *mc=new MyCoinc(0,1);

 mc->add(0,10,100);
 mc->add(0,20,200);
 mc->add(0,30,300);
 mc->add(0,40,400);
 mc->add(0,50,500);
 

 mc->add(1,32,800);
 mc->add(1,41,900);
 mc->add(1,11,600);
 mc->add(1,15,700);

 mc->add(0,60,1000);
 mc->add(0,70,1100);
 mc->add(0,80,1200);
 mc->add(0,90,1300);
 mc->add(0,100,1400);

 mc->add(1,85,1500);
 mc->add(1,89,1600);
 

 mc->print();

delete mc;

}
