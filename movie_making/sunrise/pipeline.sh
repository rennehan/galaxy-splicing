#! /bin/bash -l

snapDir=$1
snap=$2

snap=$(printf "%03d" $snap)

templatedir=/home/rennehan/configs/sunrise
#outputdir=/home/rennehan/analysis/sunrise/protocluster_lowres
#outputdir=/home/rennehan/analysis/sunrise/protocluster
outputdir=/home/rennehan/analysis/sunrise/protocluster_highres

mkdir -p $outputdir/nofilter
mkdir -p $outputdir/psf
mkdir -p $outputdir/background
mkdir -p $outputdir/rebin

#simparfile=/home/rennehan/params/EddyDiffusion/MFM/TurbOff/TurbOff_QuickDiceTest.tex
#simparfile=/home/rennehan/params/EddyDiffusion/MFM/TurbOff/TurbOff_Protocluster.tex
simparfile=/home/rennehan/params/EddyDiffusion/MFM/TurbOff/TurbOff_Protocluster_HighRes.tex

sfrhistparfile=$snapDir/sfrhist/sfrhist_$snap.params
mcrxparfile=$snapDir/mcrx/mcrx_$snap.params
broadbandparfile=$snapDir/broadband/broadband_$snap.params

mkdir -p $snapDir/restructured
mkdir -p $snapDir/sfrhist
mkdir -p $snapDir/mcrx
mkdir -p $snapDir/broadband

# We need to be inside the sunrise bin directory because we need
# units.dat which is in there.
cd /home/rennehan/software/sunrize-1.0/bin

###############################################
# Convert snapshot to SUNRISE compatible format
python /home/rennehan/analysis/sunrise/scripts/snap2sunrise.py $snapDir $snap $snap

# This is the output from snap2sunrise
snapShot=$snapDir/restructured/snapshot_$snap.hdf5

###############################################
#Run sfrhist
cp $templatedir/sfrhist.params $sfrhistparfile

sfrhistoutputfile=$snapDir/sfrhist/grid_$snap.fits
rm $sfrhistoutputfile

sed -i "s@SNAPSHOTFILE@$snapShot@g" $sfrhistparfile
sed -i "s@OUTPUTFILE@$sfrhistoutputfile@g" $sfrhistparfile
sed -i "s@SIMPARFILE@$simparfile@g" $sfrhistparfile

./sfrhist $sfrhistparfile 

###############################################
# Run MCRX
cp $templatedir/mcrx.params $mcrxparfile

mcrxinputfile=$sfrhistoutputfile
mcrxoutputfile=$snapDir/mcrx/mcrx_$snap.fits
rm $mcrxoutputfile

sed -i "s@MCRXINPUT@$mcrxinputfile@g" $mcrxparfile
sed -i "s@MCRXOUTPUT@$mcrxoutputfile@g" $mcrxparfile

./mcrx $mcrxparfile

###############################################
# Run broadband
cp $templatedir/broadband.params $broadbandparfile

broadbandinputfile=$mcrxoutputfile
broadbandoutputfile=$snapDir/broadband/broadband_$snap.fits
rm $broadbandoutputfile

sed -i "s@BROADBANDINPUT@$broadbandinputfile@g" $broadbandparfile
sed -i "s@BROADBANDOUTPUT@$broadbandoutputfile@g" $broadbandparfile

./broadband $broadbandparfile

###############################################
# Run sunpy
# Must set python path to include sunpy directory
export PYTHONPATH=$PYTHONPATH:/home/rennehan/software:/home/rennehan/software/sunpy

python /home/rennehan/software/sunpy/produce_image.py $broadbandoutputfile $outputdir $snap

exit
