echo "Setting Delphes-3.4.2 environment ... "
echo "W A R N I N G : >>------->  Do not forget to use the alias 'gcc9' [to enable c++14] before source this environment"
DELPHES_PATH=/home/gsaha/Packages/Delphes-3.4.2
ROOT_SH=/home/gsaha/Packages/buildroot/bin/thisroot.sh
ROOT_INCLUDE_PATH=$DELPHES_PATH/external/ExRootAnalysis
echo "setting exroot-analysis"
export PATH=ROOT_INCLUDE_PATH:$PATH
echo "adding delphes in LD_LIBRARY_PATH"
export LD_LIBRARY_PATH=$DELPHES_PATH:$LD_LIBRARY_PATH
echo libPath=$LD_LIBRARY_PATH
echo "Getting ROOT-6.16"
source $ROOT_SH
echo "... Done"
