#!/usr/bin/env bash

# Number of CPU cores
N=4

#-------------------------------------------------------------------------------
# Extract data
#-------------------------------------------------------------------------------
mkdir -p all
tar -xvf ../PCM/PCM_data_all.tgz -C ./all

mkdir -p unitary
tar -xvf ../PCM/PCM_data_unitary.tgz -C ./unitary

cp -R unitary unitary_tau

#-------------------------------------------------------------------------------
# Average over Langevin history (ignoring data before \tau=250)
#-------------------------------------------------------------------------------

for folder in all unitary
do
    cd $folder

    find . -name "N*1.0000e-02*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 250 -b -f
    find . -name "N*8.7500e-03*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 285 -b -f
    find . -name "N*7.5000e-03*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 333 -b -f
    find . -name "N*6.2500e-03*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 400 -b -f
    find . -name "N*5.0000e-03*s100*.dat" -print | xargs -P $N -n 1 ../../../Scripts/data_analyser.py -u -L 500 -b -f

    cd -
done

cd unitary_tau

# V = 16x16 (\tau: 100 -- 500)
########################################
find . -name "N*16x16*5.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 200 -b -f
find . -name "N*16x16*7.5000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 133 -b -f
find . -name "N*16x16*8.7500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 114 -b -f


# V = 20x20 (\tau: 75 -- 300)
########################################
find . -name "N*20x20*5.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 150 -K 600 -b -f
find . -name "N*20x20*6.2500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 120 -K 480 -b -f
find . -name "N*20x20*7.5000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 100 -K 400 -b -f
find . -name "N*20x20*8.7500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L  86 -K 343 -b -f
find . -name "N*20x20*1.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L  75 -K 300 -b -f


# V = 24x24 (\tau: 150 -- 400)
########################################
find . -name "N*24x24*5.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 300 -K 800 -b -f
find . -name "N*24x24*6.2500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 240 -K 640 -b -f
find . -name "N*24x24*7.5000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 200 -K 533 -b -f
find . -name "N*24x24*8.7500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 172 -K 458 -b -f
find . -name "N*24x24*1.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 150 -K 400 -b -f

# V = 32x32 (\tau: 170 -- 500)
########################################
find . -name "N*32x32*5.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 340 -b -f
find . -name "N*32x32*6.2500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 272 -b -f
find . -name "N*32x32*7.5000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 227 -b -f
find . -name "N*32x32*8.7500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 195 -b -f
find . -name "N*32x32*1.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 170 -b -f


# V = 48x48 (\tau: 220 -- 500)
########################################
find . -name "N*48x48*5.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 440 -b -f
#find . -name "N*48x48*6.2500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 352 -b -f
find . -name "N*48x48*7.5000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 293 -b -f
#find . -name "N*48x48*8.7500e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 252 -b -f
find . -name "N*48x48*1.0000e*s100*.dat" -print | xargs -P $N -n 1 ../../../scripts/data_analyser.py -u -L 220 -b -f


cd -

#-------------------------------------------------------------------------------
# Extrapolation to vanishing \epsilon
#-------------------------------------------------------------------------------

for folder in all unitary unitary_tau
do
    cd $folder

    mkdir -p tmpN{03..06}
    mkdir -p tmpN12

    echo {03..06} 12 | xargs  -I {} -P $N -n 1 sh -c 'mv -v N{}*.stat tmpN{}'
    echo {03..06} 12 | xargs  -I {} -P $N -n 1 sh -c '../../../scripts/pcm_analyser.py -d tmpN{}'
    echo {03..06} 12 | xargs  -I {} -P $N -n 1 sh -c '../../../scripts/make_eps0.00_all_v.py -d tmpN{}'

    cd -
done

#-------------------------------------------------------------------------------
# Extrapolation to infinite volume
#-------------------------------------------------------------------------------

#printf "03 04 05 06 12 "
for folder in all unitary unitary_tau
do
    cd $folder

    echo {03..06} 12 | xargs  -I {} -P $N -n 1 sh -c 'cd ./tmpN{} &&  ../../../../scripts/fit.py -z'

    cd -
done

#-------------------------------------------------------------------------------
# Calculate ratios
#-------------------------------------------------------------------------------

for folder in all unitary unitary_tau
do
    cd $folder
    # Collect data for all N
    mkdir -p results
    echo {03..06} 12 | xargs  -I {} -P $N -n 1 sh -c 'mv ./tmpN{}/*.*at ./results/'
    cd results
    ../../../../scripts/pcm_ratios.py
    ../../../../scripts/renormalon_fit.py

    cd ../../
done
