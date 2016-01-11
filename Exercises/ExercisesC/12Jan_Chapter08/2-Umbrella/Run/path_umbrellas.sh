#! /bin/bash
bins=(-10 -8 -6 -4 -2 0 2 4 6 8)
for i in "${bins[@]}"
do
    let j=i+1
    echo "1000 1 "$i" "$j | ./run 
    cp Umbrella.dat "L_"$i"_R_"$j".dat"
done 
