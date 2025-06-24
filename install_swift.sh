#!/bin/bash
 
echo " " 
echo " " 
echo "Installing the swift N-body lib"
echo " " 
echo " " 

cd ./swift_j/;
awk -v a="$PWD" '{ if (NR == 3) print "set SWIFT_DIR="a; else print $0}' @make_temp > @make;
csh @makeall;
#cp libswift.a ../;
cd ../;


echo " " 
echo " " 
echo "Compiling symba5_j_migrate6 and other N-body routines!"
echo " " 
echo " " 

gfortran -O3 ./source/swift_symba5_j_migrate6.f -o swift_symba5_j_migrate6 ./swift_j/libswift.a;
gfortran -O3 ./source/geninit_j3.f -o geninit_j3 ./swift_j/libswift.a;
gfortran -O3 ./source/follow_symba2.f -o follow_symba2 ./swift_j/libswift.a;
 




