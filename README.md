# Delphes-Analysis FrameWork 

C++, ROOT and Delphes based code setup. It is mandatory to install `Delphes-3.4.2` before. 

## Instructions : 
 - Check the setup.sh file. Use correct paths of Delphes and ROOT
 - `source setup.sh`
 - Use the correct `Delphes` path in `Makefile.exoanalysis`
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