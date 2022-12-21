// Copyright 2020 Skoltech
#include "pcscl_impl.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Read decoder parameters
    unsigned int blockLenLog = (unsigned int)*(mxGetPr(prhs[0]));
    unsigned int blockLen = 1 << blockLenLog;
    unsigned int listSize = (int)*(mxGetPr(prhs[1]));
    bool use_pathMetricApprox = (bool)*(mxGetPr(prhs[2]));
    bool useMinSum = (bool)*(mxGetPr(prhs[3]));
    // Instantiate decoder
    PolarSCL polar_scl = PolarSCL(blockLenLog,
                                  listSize,
                                  useMinSum,
                                  use_pathMetricApprox);
    // Read frozen bits
    double* input = mxGetPr(prhs[4]);
    std::vector<PolarPath::bit_t> f(blockLen);
    for (unsigned int i = 0; i < blockLen; ++i) {
        f[i] = (PolarPath::bit_t) input[i];
    }
    // Read LLR vector
    input = mxGetPr(prhs[5]);

    std::vector<PolarPath::llr_t> llr_in(blockLen);
    for (unsigned int i = 0; i < blockLen; ++i) {
        llr_in[i] = (PolarPath::llr_t) input[i];
    }
    PolarSCL::DecoderOutput decoderOutput = polar_scl.Decode(llr_in, f);

    // Return results to MATLAB
    // bits list
    plhs[0] = mxCreateDoubleMatrix(listSize, blockLen, mxREAL);
    double* output = mxGetPr(plhs[0]);
      for (unsigned int i = 0; i < blockLen; ++i) {
        for (unsigned int j = 0; j < listSize; ++j) {
        output[i * listSize + j] = (double) decoderOutput.bits[j][i];
      }
    }
    // LLR
    plhs[1] = mxCreateDoubleMatrix(listSize, blockLen, mxREAL);
    double* output1 = mxGetPr(plhs[1]);
    for (int i = 0; i < blockLen; ++i) {
        for (int j = 0; j < listSize; ++j) {
            output1[i * listSize + j] = (double) decoderOutput.llrs[j][i];
        }
    }
    // path metrics
    plhs[2] = mxCreateDoubleMatrix(listSize, 1, mxREAL);
    double* output3 = mxGetPr(plhs[2]);
    for (unsigned int i = 0; i < listSize; i++) {
      output3[i] = (double) decoderOutput.pMetrics[i];
    }
}
