// Copyright 2020 Skoltech
#include "pcscl_impl.h"

extern "C"
void decode(double* llrChannel,
        uint8_t* isFrozen,
        uint8_t* decodedBits,
        double* outLLRs, double* pathMetrics, unsigned int n, unsigned int L) {

    unsigned int N = 1ULL << n;
    std::vector<double> llrChannelVec (llrChannel, llrChannel + N);
    std::vector<uint8_t> isFrozenVec (isFrozen, isFrozen + N);
    
    PolarSCL polar_scl = PolarSCL(n, L, true, true);
    PolarSCL::DecoderOutput decoderOutput = polar_scl.Decode(llrChannelVec, isFrozenVec);

    for(unsigned int i = 0; i < L; i ++) {
        for(unsigned int j = 0; j < N; j ++) {
            // Copy output LLRs
            outLLRs[i * N + j] = decoderOutput.llrs[i][j];
            // Copy decoded bits
            decodedBits[i * N + j] = decoderOutput.bits[i][j];
        }
        // Copy path metrics
        pathMetrics[i] = decoderOutput.pMetrics[i];
    }
}

