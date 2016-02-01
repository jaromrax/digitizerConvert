
#if !defined( __CINT__) || defined (__MAKECINT__)
#endif

#if !defined(COINONLINE_H)
#define COINONLINE_H

// 0   1 just ~    2 allstatus
#define DEBUG 0
//======= rotate ===
#include<algorithm>
#include "TH2F.h"


static const int BMAX=1000;
static const int BSIZE=56;  // buffer size can be 3*9 +
static const int WIN=100;   // TAC window -100 .. +100

using namespace std;



class OneChan{
  Int_t chan;
  int max;
  //------------ BUFFERS-----------
  unsigned long long int time[BMAX];
  Int_t value[BMAX];  //energy
  Int_t efit[BMAX];  
  Int_t xi[BMAX];  
  Int_t s[BMAX];  


  
 public:

  OneChan(int ch );
  ~OneChan();

  //  int push(unsigned long long int t, int val);
  int push(unsigned long long int t,int val,
	   int vefit,int vxi,int vs);
  int print();
  int remove(int i);


  int getmax();
  unsigned long long int gett(int i);
  int getv(int i);
  int getvf(int i);
  int getvs(int i);


  
  //  ClassDef(OneChan,0);
};




//======================== COINC
//======================== COINC

class MyCoinc{
 public:
  // there are 4 queues:
  
  int max;
  int totalsum;
  OneChan *OC0;
  OneChan *OC1;
  TH2F *bidim;
  TH1F *hidt;
  Int_t ch0last,ch1last,dtlast,  ch0flast,ch1flast, ch0slast,ch1slast;
  
  MyCoinc(int ch0, int ch1);
  ~MyCoinc();
  int add(int ch, unsigned long long int t, int val);
  int add(int ch, unsigned long long int t, int val,int vefit,int vxi,int vs);
  int print();
  TH2F* getbidim();
  TH1F* gethidt();
  
  Int_t getlastC0();
  Int_t getlastC1();
  Int_t getlastDT();
  Int_t getlastC0fit();
  Int_t getlastC1fit();
  Int_t getlastC0s();
  Int_t getlastC1s();
  //  int getmax();
  //  unsigned long long int gett(int i);

  
  //  ClassDef(MyCoinc,0)
};

//MyCoinc *mc;
#endif
