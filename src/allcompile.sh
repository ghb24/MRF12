#!/bin/bash -x

touch *.f90
cd ../
make
touch src/*.f90
make
