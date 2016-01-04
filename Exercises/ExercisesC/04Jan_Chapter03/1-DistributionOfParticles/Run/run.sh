for N in 100 300 500 700 900 1000; do
    mkdir $N
    cd $N
    time ../../Source/distribution $N 2 10
    cd ..
done
