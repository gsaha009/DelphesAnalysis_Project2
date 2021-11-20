echo "Setting Delphes-3.4.2 environment ... "
echo "W A R N I N G : >>------->  Do not forget to use the alias 'gcc9' [to enable c++14] before source this environment"
#conda deactivate
APPDIR=/home/gsaha
WORKDIR=$APPDIR/Work/DelphesAnalysis_2
ROOTPATH=/home/gsaha/Packages/buildroot/bin/thisroot.sh
cd $WORKDIR
export LD_LIBRARY_PATH=/home/gsaha/Packages/Delphes-3.4.2/:$LD_LIBRARY_PATH
echo libPath=$LD_LIBRARY_PATH
echo "Getting ROOT-6.16"
source $ROOTPATH
echo "... Done"
