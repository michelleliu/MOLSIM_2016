#! /bin/bash
bins=(-5 -4 -3 -2 -1 0 1 2 3 4 5)

cat bounds | while read B
do
    b=( ${B} )
    echo ${b[0]} ${b[1]}
    echo 1000 0.1 ${b[0]} ${b[1]} | ./run
    echo writing results to ${b[0]}.dat
    cp Umbrella.dat umbrella${b[0]}.dat
done
