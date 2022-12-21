#!/bin/bash
g++ -std=c++11  -c ../pcscl_impl.cpp valgrind.cpp
g++ valgrind.o pcscl_impl.o
valgrind ./a.out
rm *.o *.out
