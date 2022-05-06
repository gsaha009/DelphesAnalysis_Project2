# Delphes-Analysis FrameWork 

C++ and ROOT based code setup.

## Instructions : 
 - Check the setup.sh file. Use correct paths of Delphes and ROOT
 - `source setup.sh`
 - Change the paths in `DelphesClasses.h/cc` and `ExoAnalysis.h`
 - ... Ready for compilation

## How to compile?

```
make clean -f Makefile.exoanalysis
make cling -f Makefile.exoanalysis
make -f Makefile.exoanalysis
``` 

## How to run?
```
./FCNC.exe test.job
```