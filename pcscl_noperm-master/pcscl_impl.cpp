// Copyright 2020 Skoltech
#include "pcscl_impl.h"
/*
 * ******************** Aux. function implementation ********************
 */

inline uint32_t bitReverse(uint32_t x, uint32_t b) {
    // Fast bit reversal by D. Knuth
    x = x - 1;
    uint32_t t;
    x = (x << 15) | (x >> 17);
    t = (x ^ (x >> 10)) & 0x003f801f;
    x = (t + (t << 10)) ^ x;
    t = (x ^ (x >> 4)) & 0x0e038421;
    x = (t + (t << 4)) ^ x;
    t = (x ^ (x >> 2)) & 0x22488842;
    x = (t + (t << 2)) ^ x;
    x = x >> (32 - b);
    return x + 1;
}

inline PolarPath::index_t logInteger(PolarPath::index_t x) {
    // return floor(log2(x))
    PolarPath::index_t xLog = 0;
    while (x) {
        x >>= 1;
        xLog = xLog + 1;
    }
    return xLog - 1;
}

inline PolarPath::llr_t sign(const PolarPath::llr_t a) {
    return (a > 0) ? 1.0 : -1.0;
}

inline PolarPath::pm_t jacobLog(PolarPath::llr_t x) {
    if (x < -17.0) {
        return 0;
    }
    if (x > 15.0) {
        return x;
    } else {
        return log(1 + exp(x));
    }
}

inline PolarPath::pm_t pmUpdateApprox(const PolarPath::bit_t bit,
                               const PolarPath::llr_t llr) {
    PolarPath::bit_t bitFromLLR = llr < 0 ? 1 : 0;
    if (bitFromLLR != bit) {
        return -std::abs(llr);
    }
    return 0.0;
}

inline PolarPath::pm_t pmUpdateExact(const PolarPath::bit_t bit,
                              const PolarPath::llr_t llr) {
    PolarPath::pm_t pathMetric = - jacobLog(llr);
    if (bit == 0) {
        pathMetric += llr;
    }
    return pathMetric;
}

/*
 * ******************** Polar path implementation ********************
 */
PolarPath::llr_t PolarPath::cnopMS(const PolarPath::llr_t l1,
                                   const PolarPath::llr_t l2) {
    return (sign(l1) * sign(l2)) * std::min(std::abs(l1), std::abs(l2));
}

PolarPath::llr_t PolarPath::cnopSP(const PolarPath::llr_t l1,
                                   const PolarPath::llr_t l2) {
    llr_t a1 = std::tanh(l1 / 2);
    llr_t a2 = std::tanh(l2 / 2);

    llr_t ath;
    if (std::abs(a1 * a2) < 1) {
        ath  = 2 * std::atanh(a1 * a2);
    } else {
        ath = sign(l1) * sign(l2) * std::min(std::abs(l1), std::abs(l2));
    }
    return ath;
}

PolarPath::llr_t PolarPath::cnop(const PolarPath::llr_t l1,
                                 const PolarPath::llr_t l2) {
    return m_useMinSum ? cnopMS(l1, l2) : cnopSP(l1, l2);
}

PolarPath::llr_t PolarPath::vnop(const PolarPath::bit_t bit,
                                 const PolarPath::llr_t l1,
                                 const PolarPath::llr_t l2) {
    //  (-1) ^ bit * l1 + l2;
    return (bit & 1) ? (-l1 + l2) : (l1 + l2);
}

PolarPath::PolarPath(PolarPath::index_t blockLenLog,
                     bool useMinSum, bool usePmApprox) :
    m_blockLen((index_t) 1 << blockLenLog),
    m_blockLenLog(blockLenLog),
    m_useMinSum(useMinSum),
    m_usePmApprox(usePmApprox),
    m_pathMetric(0.0),
    m_llrInternal(std::vector<llr_t> (2 * m_blockLen)),
    m_llrFinal(std::vector<llr_t> (m_blockLen)),
    m_bitsInternal(std::vector<bit_t> (2 * m_blockLen)),
    m_bitsFinal(std::vector<bit_t> (m_blockLen))
{};


void
PolarPath::setChannelLLR(const std::vector<PolarPath::llr_t> &llrChannel) {
    for (unsigned int i = 0; i < m_blockLen; i++) {
        index_t a = bitReverse(i + 1, m_blockLenLog);
        m_llrInternal[i + m_blockLen] = llrChannel[a - 1];
    }
}

void PolarPath::copyFrom(const PolarPath& other) {
    std::copy(other.m_llrInternal.begin(), other.m_llrInternal.end(),
              m_llrInternal.begin());
    std::copy(other.m_llrFinal.begin(), other.m_llrFinal.end(),
              m_llrFinal.begin());
    std::copy(other.m_bitsInternal.begin(), other.m_bitsInternal.end(),
              m_bitsInternal.begin());
    std::copy(other.m_bitsFinal.begin(), other.m_bitsFinal.end(),
              m_bitsFinal.begin());
    m_pathMetric = other.m_pathMetric;
}

void PolarPath::propagateBits(PolarPath::index_t bitNum) {
    m_bitsInternal[2 + (bitNum & (index_t)1)] = m_bitsFinal[bitNum];
    if (!(bitNum & (index_t)1)) {
        // Nothing to update for even bits
        return;
    }
    index_t bitLevel = 0;
    for (index_t i = 0; i <= m_blockLenLog; i++) {
        if ((bitNum & ((index_t)1 << i)) == 0) {
            break;
        }
        bitLevel += 1;
    }
    if (bitLevel == m_blockLenLog) {
        bitLevel = 1;
    }
    index_t indexBound = (index_t)1 << bitLevel;
    for (index_t i = 2; i < indexBound; i += 2) {
        m_bitsInternal[i * 2 + 1] = m_bitsInternal[i] ^ m_bitsInternal[i + 1];
        m_bitsInternal[(i + 1) * 2 + 1] = m_bitsInternal[i + 1];
    }
    for (index_t i = indexBound; i < 2 * indexBound; i += 2) {
        m_bitsInternal[i * 2] = m_bitsInternal[i] ^ m_bitsInternal[i + 1];
        m_bitsInternal[(i + 1) * 2] = m_bitsInternal[i + 1];
    }
}

void PolarPath::propagateLLR(PolarPath::index_t bitNum) {
    index_t id = bitNum ^ (bitNum > 0 ? bitNum - 1 : m_blockLen - 1);
    // Find the maximum power of two smaller than id
    index_t nextLevel = (index_t)1 << logInteger(id);
    bool lGfunc = bitNum & nextLevel;
    for (index_t i = (nextLevel * 2 - 1); i > 0; i--) {
        if (lGfunc) {
            m_llrInternal[i] = vnop(m_bitsInternal[2 * i],
                                    m_llrInternal[2 * i],
                                    m_llrInternal[2 * i + 1]);
        } else {
            m_llrInternal[i] = cnop(m_llrInternal[2 * i],
                                    m_llrInternal[2 * i + 1]);
        }
        if (i == nextLevel) {
            nextLevel >>= 1;
            lGfunc = bitNum & nextLevel;
        }
    }
    m_llrFinal[bitNum] = m_llrInternal[1];
}

void
PolarPath::makeBitDecision(PolarPath::bit_t bit,
                           PolarPath::index_t bitNum) {
    // Update bit values
    m_bitsFinal[bitNum] = bit;
    // Update path metrics
    if (m_usePmApprox) {
        m_pathMetric += pmUpdateApprox(bit, m_llrFinal[bitNum]);
    } else {
        m_pathMetric += pmUpdateExact(bit, m_llrFinal[bitNum]);
    }
}

std::vector<PolarPath::llr_t> PolarPath::getOutputLLRs() const {
    return m_llrFinal;
}
std::vector<PolarPath::bit_t> PolarPath::getDecodedBits() const {
    return m_bitsFinal;
}
PolarPath::pm_t PolarPath::getPathMetric() const {
    return m_pathMetric;
}


/*
 * ******************** Polar SCL implementation ********************
 */

PolarSCL::PolarSCL(unsigned int blockLenLog,    // Block length log
                   unsigned int listSize,       // Output list size
                   bool useMinSum,              // Check node operation
                   bool usePathMetricApprox) :  // Path metric type
        m_blockLen(1 << blockLenLog),
        m_listSize(listSize),
        m_listSizeCur(1) {
    // Note that maximum lst size is the next power-of-two number greater
    // than doubled list size, initial list size values are powers of two
    m_paths = std::vector<PolarPath> (
            2 * std::pow(2, std::ceil(std::log2(listSize))),
            PolarPath(blockLenLog, useMinSum, usePathMetricApprox));
}

PolarSCL::DecoderOutput PolarSCL::Decode(
        const std::vector <PolarPath::llr_t> &llrChannel,
        const std::vector <PolarPath::bit_t> &isFrozen
        ) {
    m_listSizeCur = 1;
    // Set LLRs to the first path
    m_paths[0].setChannelLLR(llrChannel);

    for (unsigned int bitNum = 0; bitNum < m_blockLen; bitNum++) {
        // Propagate LLR to make bit decision
        propagateLLR(bitNum);
        // Perform single SC step
        if (isFrozen[bitNum] == 1) {
            processFrozenBit(bitNum);
        } else {
            processInfromationBit(bitNum);
        }
        // Propagate bits to construct information word
        propagateBits(bitNum);
    }
    // Construct the output
    return getOutput();
}

void PolarSCL::propagateLLR(unsigned int bitNum) {
    for (unsigned int j = 0; j < m_listSizeCur; j++) {
        m_paths[j].propagateLLR(bitNum);
    }
}

void PolarSCL::processInfromationBit(unsigned int bitNum) {
    for (unsigned int j = 0; j < m_listSizeCur; j++) {
        // clone path arrays
        m_paths[j + m_listSizeCur].copyFrom(m_paths[j]);
        // Make bit decisions (two hypotheses)
        m_paths[j].makeBitDecision(0, bitNum);
        m_paths[j + m_listSizeCur].makeBitDecision(1, bitNum);
    }
    m_listSizeCur *= 2;
    // Perform list pruning if required
    if (m_listSizeCur >= (2 * m_listSize)) {
        std::sort(m_paths.begin(), m_paths.begin() + m_listSizeCur);
        m_listSizeCur = m_listSize;
    }
}

void PolarSCL::processFrozenBit(unsigned int bitNum) {
    for (unsigned int j = 0; j < m_listSizeCur; j++) {
        m_paths[j].makeBitDecision(0, bitNum);
    }
}

void PolarSCL::propagateBits(unsigned int bitNum) {
    for (unsigned int j = 0; j < m_listSizeCur; j++) {
        m_paths[j].propagateBits(bitNum);
    }
}

PolarSCL::DecoderOutput PolarSCL::getOutput() const {
    DecoderOutput output;
    output.bits = std::vector<std::vector<PolarPath::bit_t> > (m_listSize);
    output.llrs = std::vector<std::vector<PolarPath::llr_t> > (m_listSize);
    output.pMetrics = std::vector<PolarPath::pm_t> (m_listSize);
    for (unsigned int i = 0; i < m_listSize; i++) {
        output.bits[i]     = m_paths[i].getDecodedBits();
        output.llrs[i]     = m_paths[i].getOutputLLRs();
        output.pMetrics[i] = m_paths[i].getPathMetric();
    }
    return output;
}
