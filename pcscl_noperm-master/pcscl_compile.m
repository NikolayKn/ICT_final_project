% Select compiler options
if ismac
    % No compiler options required for MAC platform
    compiler_option = '';
elseif isunix
    % Code to run on Linux platform
    compiler_option = '-march=native -Ofast -Wl,--stack,4194304';
elseif ispc
    % Code to run on Windows platform
    compiler_option = '/Ox /GL /fp:fast /arch:AVX2 /F 4194304';
else
    compiler_option = '';
    error('Platform not supported');
end
mex( ...
    'pcscl_impl.cpp',  'pcscl_noperm.cpp', ...
    ['COMPFLAGS="$COMPFLAGS ' compiler_option '"'] ...
    );