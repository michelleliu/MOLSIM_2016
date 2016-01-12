#! /bin/bash
ttype=$1
cp ../w.${ttype}.dat w.dat
srcdir=/Users/michelle/MOLSIM_2016/Exercises/ExercisesC/11Jan_Chapter07/3-IsingModel
N=32
Beta=0.5
Ncycle=100
echo $N $Beta $Ncycle | time $srcdir/Source/ising
exit
