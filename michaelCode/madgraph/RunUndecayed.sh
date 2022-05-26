#!/bin/bash

#usage ./Run.sh mvGeV ebeamGeV runNumber outSubDir nevents

CWD=`pwd`
RANDOM=`date +%N|sed s/...$//`
SCRATCH=$CWD/scratch.$3
#SCRATCH=$HOME/Local

echo "=== starting job in directory $CWD ==="

OUTDIR=$CWD/Repository
mv=$1
energy=$((2*$2+10))
run=$3
outSubDir=$4
nevnts=$5
mkdir $SCRATCH
cp AprimeWUndecayed.tar.gz $SCRATCH
cd $SCRATCH
tar -xzf AprimeWUndecayed.tar.gz
cd AprimeWUndecayed

line=622" "$mv" ""# APMASS"
sed -in "s/622.*# APMASS/$line/" Cards/param_card.dat

line="Anuc = 207.0" 
sed -in "s/Anuc \= 184.0/$line/" Source/MODEL/couplings.f 

line="2  82.0000000e+00 # Znuc    ,nuclear charge" 
sed -in "s/2.*# Znuc    ,nuclear charge/$line/" Cards/param_card.dat 

line=$nevnts" = nevents ! Number of unweighted events requested"
sed -in "s/.*nevents.*/$line/" Cards/run_card.dat

line=$run" = iseed ! rnd seed (0=assigned automatically=default))"
sed -in "s/.*iseed.*/$line/" Cards/run_card.dat

line=$energy" = ebeam1  ! beam 1 energy in GeV"
sed -in "s/.*ebeam1.*/$line/" Cards/run_card.dat

line="192.51 = ebeam2 ! beam 2 energy in GeV"
sed -in "s/.*ebeam2.*/$line/" Cards/run_card.dat

line="192.51 = mbeam2 ! beam 2 energy in GeV"
sed -in "s/.*mbeam2.*/$line/" Cards/run_card.dat

line="623 192.51 #HPMASS"
sed -in "s/623.*# HPMASS (tungsten)/$line/" Cards/param_card.dat

line="F77 = f95 -std=legacy"                                                                    
sed -in "s/.*f77.*/$line/" Source/makefile  

line="FC       = f95 -std=legacy"  
sed -in "s/.*f77.*/$line/" Source/DHELAS/Makefile   

line="F77 = f95 -std=legacy" 
sed -in "s/.*f77.*/$line/" Source/PDF/makefile 

line="F77 = f95 -std=legacy" 
sed -in "s/.*f77.*/$line/" Source/MODEL/makefile  

line="F77 = f95 -std=legacy" 
sed -in "s/.*f77.*/$line/" SubProcesses/makefile

line="F77 = f95 -std=legacy" 
sed -in "s/.*f77.*/$line/" SubProcesses/P0_e-n+_e-n+x/makefile   

line="F77 = f95 -std=legacy" 
sed -in "s/.*f77.*/$line/" Source/makefile_dynamic 

./bin/generate_events 0 LDMX_Pb_UndecayedAP."$energy"GeV.W.mA.$mv.$run
mkdir -p $OUTDIR/$outSubDir
mv Events/*unweighted*gz $OUTDIR/$outSubDir
cd $CWD
rm -rf $SCRATCH
cd Repository/$outSubDir
gunzip LDMX_Pb_UndecayedAP."$energy"GeV.W.mA."$mv"."$run"_unweighted_events.lhe.gz
#python $CWD/LHE_Root.py LDMX_W_UndecayedAP."$energy"GeV.W.mA."$mv"."$run"_unweighted_events.lhe
#rm LDMX_W_UndecayedAP."$energy"GeV.W.mA."$mv"."$run"_unweighted_events.lhe



# -lh /local/cms/user/revering/madgraph/scratch  | grep revering
