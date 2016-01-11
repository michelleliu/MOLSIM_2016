#! /bin/bash
bins=(-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9)
count=1
for i in "${bins[@]}"
do
    let j=i+2
    echo "("$i" , "$j")"
    echo "1000 0.1 "$i" "$j | ./run 
    echo "Writing results to "$count".dat"
    cp Umbrella.dat $count".dat"
    let count=$count+1
done 
