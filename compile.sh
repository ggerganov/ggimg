#!/bin/bash

cur="ggimg"
echo "Compiling ${cur} ... "
g++ -std=c++11 -O3 -I. examples/${cur}.cpp -o ${cur}
