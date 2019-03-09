#!/bin/bash

cur="ggimg"
echo "Compiling ${cur} ... "
g++ -std=c++11 -Wall -O3 -I. examples/${cur}.cpp -o ${cur}

cur="ggimg-mt"
echo "Compiling ${cur} ... "
g++ -std=c++11 -Wall -O3 -I. examples/${cur}.cpp -o ${cur}
