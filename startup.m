cd cpp_impl;
% cd pcscl_noperm-master;
mex pcscl_jsc.cpp
% mex *.cpp
cd ..;
addpath('cpp_impl');
% addpath('pcscl_noperm-master');
addpath('tools');
fprintf('Startup routines completed.\n');