# digitizerConvert

here are **root** scripts to convert V1724 and DT5780 data:

**convDat** ... convert ASCII files from listmode.  DT5780 (from MC2 software) is recognized by the word(s) HEADER in the file,  V1724 is recognized from the initial comment **\#** with starting time. 

First parameter is filename

Second parameter is 0/1 :  0=do not sort the file before;  1=sort the file before

For a good performance, load with

*.L convDat.C+*

and run like - if you do not care about coincidences -

convDat("run001.dat",0)



**convWav.C** ... conversion of wav format created by gregory using V1724 data.


**n42totxt.C** ... conversion of n42 format of MC2. *simplexml* code must be installed.

