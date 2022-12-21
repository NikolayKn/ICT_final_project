// Copyright 2020 Skoltech
#ifndef PCSCL_IMPL_H_
#define PCSCL_IMPL_H_

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cstdint>
#ifndef _WIN32
#define __forceinline __attribute__((always_inline))
#endif

class PolarPath {
 public:
    // Define types for indices, LLRs and bits
    typedef uint32_t index_t;
    typedef uint8_t bit_t;
    typedef double llr_t;
    typedef double pm_t;
    // Initialize path with the blocklength and node operation parameters
    PolarPath(index_t blockLenLog, bool useMinSum, bool usePmApprox);
    // Initialize a path with channel LLR vector
    void setChannelLLR(const std::vector<llr_t> &llrChannel);
    void propagateLLR(index_t bitNum);
    void propagateBits(index_t bitNum);
    void makeBitDecision(bit_t bit, index_t bitNum);
    void copyFrom(const PolarPath& other);
    // get decoding results
    std::vector<llr_t> getOutputLLRs() const;
    std::vector<bit_t> getDecodedBits() const;
    pm_t getPathMetric() const;
    // Max path metric corresponds to the best (minimum) path
    bool operator < (const PolarPath &other) const
        {return m_pathMetric > other.m_pathMetric;};

 private:
    // Min-sum approximation for check node operation (cnop)
    static llr_t cnopMS(const llr_t l1, const llr_t l2);
    // Sum-product approximation for check node operation (cnop)
    static llr_t cnopSP(const llr_t l1, const llr_t l2);
    // Variable node operation: (-1)^bit * l1 + l2
    static llr_t vnop(const bit_t bit, const llr_t l1, const llr_t l2);

    llr_t cnop(const llr_t, const llr_t);

    index_t m_blockLen;
    // The depth of recusrive polar transform block length = 2 ^ depth;
    index_t m_blockLenLog;
    // If true, use min-sum approximation for Check node operation
    bool m_useMinSum;
    bool m_usePmApprox;
    // Path metric for a given path
    pm_t m_pathMetric;
    // Array of length 2 * blocklength for LLR propagation procedure
    std::vector<llr_t> m_llrInternal;  // Used by propagateLLR only
    // Array of length blocklength representingthe information word LLRs
    std::vector<llr_t> m_llrFinal;
    // Array of length 2 * blocklength for bit propagation procedure
    std::vector<bit_t> m_bitsInternal;  // Used by propagateBits only
    // Array of length blocklength representing the information word bits
    std::vector<bit_t> m_bitsFinal;
};

class PolarSCL {
 public:
    // Define the decoder output. First dimension of each output filed
    // has <listSize> entities
    typedef struct {
        // 2D array of decoded bits corresponding to the information word
        std::vector<std::vector<PolarPath::bit_t> > bits;
        // 2D array of LLRs corresponding to information bits
        std::vector<std::vector<PolarPath::llr_t> > llrs;
        // Path metrics for each decoded information word
        std::vector<PolarPath::pm_t> pMetrics;
    } DecoderOutput;
    // Use the logarithm of the block length and the list size for init.
    PolarSCL(
            unsigned int blockLenLog,    // Log2 of the block length
            unsigned int listSize,       // Output list size
            bool useMinSum,              // Min-Sum / Sum-Prod check node
            bool use_pathMetricApprox);  // Exact or approx path metric

    // Decoder: pass channel LLR vector and frozen bit pattern.
    DecoderOutput Decode(
            const std::vector<PolarPath::llr_t> &llrChannel,
            const std::vector<PolarPath::bit_t> &isFrozen);

 private:
    // Succesive cancellation for every bit consists of four steps:
    void propagateLLR(unsigned int bitNum);
    void processInfromationBit(unsigned int bitNum);
    void processFrozenBit(unsigned int bitNum);
    void propagateBits(unsigned int bitNum);
    // Return information word bits, corresponding LLRs and path metrics
    PolarSCL::DecoderOutput getOutput() const;

 private:
    unsigned int m_blockLen;
    // Output list size
    unsigned int m_listSize;
    // Current list size
    unsigned int m_listSizeCur;
    std::vector<PolarPath> m_paths;
};
#endif  //  PCSCL_IMPL_H_
