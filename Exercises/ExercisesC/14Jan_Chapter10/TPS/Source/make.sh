for x in 1.75 1.85 1.9 1.95 2.0 2.05 2.10 2.15 2.25
do
make clean; export angle=$x; make
cp tps.x tps.x-$x

done
