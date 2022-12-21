#include<stdio.h>
#include "pcscl_impl.h"
int main()
{
    // Create k = 21, n = 64 frozen pattern
    unsigned int N = 64;
    unsigned int m = 6; // log2(N)
    std::vector<PolarPath::bit_t> f = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    unsigned int L = 4096; // List size
    // Hard-coded input LLR 
    std::vector<PolarPath::llr_t> llr_in = {-7.48, -0.92, 2.26, 11.57, 6.11, -5.65, 5.34, -1.57, 2.69, 4.84, 6.42, -5.17, -12.61, 6.62, 6.21, -4.97, -0.85, -14.72, 1.87, -11.98, -8.02, -6.28, 15.80, 5.65, 3.47, 9.54, -5.22, -8.16, 3.42, 2.40, -3.50, -1.52, -0.45, 3.32, -3.49, -5.60, 1.53, -9.96, -9.95, -1.41, -1.62, -5.98, -3.92, -6.96, 2.78, 5.88, 10.74, 12.17, -1.33, -1.83, 1.98, 7.22, -0.58, 1.13, 0.15, 4.19, -2.71, 1.21, -14.30, -3.58, 2.43, 7.24, 4.44, 17.99};
    PolarSCL polar_scl = PolarSCL(m, L, true, true);
    PolarSCL::DecoderOutput decoder_output = polar_scl.Decode(llr_in, f);
    return 0;
}
