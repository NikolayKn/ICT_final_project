//#include "mex_commands.h"
//#include "profiler.h"
#include "mex.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>
#include <cmath>

using namespace std;

#define ABS(A) (((A) >= 0) ? (A) : -(A))
void ntt(vector<double>& p); // Multi-dimensional Fourier Transform

//=======================================================================================================================================
//===================================================== Log distribution ================================================================
//=======================================================================================================================================

class LogDistr
{
public:
    int m_K;
    int m_Q; // 2^K (number of messages)
    vector<double> m_logProbs;
    vector<bool> m_isnegative;
public:
    LogDistr()
        : MAX_EXP(exp((double) 700))
        , MIN_EXP(exp((double)-700))
    {
        m_K = 0;
        m_Q = 0;
    }
    LogDistr(int K)
        : MAX_EXP(exp((double) 700))
        , MIN_EXP(exp((double)-700))
    {
        m_K = K;
        m_Q = 1 << K;
        m_logProbs = vector<double>(m_Q, 0);
        m_isnegative = vector<bool>(m_Q, false);
    }

    LogDistr(const LogDistr& logDistr)
        : MAX_EXP(exp((double) 700))
        , MIN_EXP(exp((double)-700))
    {
        m_K = logDistr.m_K;
        m_Q = logDistr.m_Q;
        m_logProbs = logDistr.m_logProbs;
        m_isnegative = logDistr.m_isnegative;
    }
    const LogDistr& operator=(const LogDistr& logDistr)
    {
        m_K = logDistr.m_K;
        m_Q = logDistr.m_Q;
        m_logProbs = logDistr.m_logProbs;
        m_isnegative = logDistr.m_isnegative;
        return *this;
    }
public:
    int getHard() const
    {
        double max = m_logProbs[0];
        int    index = 0;

        for (int i = 1; i < m_Q; ++i)
        {
            if (max <= m_logProbs[i])
            {
                max = m_logProbs[i];
                index = i;
            }
        }

        return index;
    }
    
    /* permutation */
    LogDistr operator+(int value) const
    {
        value %= m_Q;
        LogDistr result(m_K);

        for(int i = 0; i < m_Q; ++i)
        {
            result.m_logProbs[i^value] = m_logProbs[i];
        }

        return result;
    }

    /* vnop operation */
    LogDistr operator+(const LogDistr& logDistr) const
    {
        LogDistr result(m_K);

        for(int i = 0; i < m_Q; ++i)
        {
            result.m_logProbs[i] = m_logProbs[i] + logDistr.m_logProbs[i];
            result.m_isnegative[i] = (m_isnegative[i] != logDistr.m_isnegative[i]);
        }

        return result;
    }
    const LogDistr& operator+=(const LogDistr& logDistr)
    {
        for(int i = 0; i < m_Q; ++i)
        {
            m_logProbs[i] += logDistr.m_logProbs[i];
            m_isnegative[i] = (m_isnegative[i] != logDistr.m_isnegative[i]);
        }

        return *this;
    }
    LogDistr operator-(double value) const
    {
        LogDistr result(m_K);

        for(int i = 0; i < m_Q; ++i)
        {
            result.m_logProbs[i] = m_logProbs[i] - value;
        }

        return result;
    }
    const LogDistr& operator-=(double value)
    {
        for(int i = 0; i < m_Q; ++i)
        {
            m_logProbs[i] -= value;
        }

        return *this;
    }
    const LogDistr& normalize()
    {
        double sum = 0;

        for (int i = 0; i < m_Q; ++i)
        {
            sum += safe_exp(m_logProbs[i]);
        }
        if (sum <= 0)
        {
            sum = m_Q;
        }
        sum = safe_log(sum);
        for (int i = 0; i < m_Q; ++i)
        {
            m_logProbs[i] -= sum;
        }
        return *this;
    }
    const LogDistr& sub_max()
    {
        int index = getHard();
        *this -= m_logProbs[index];
        return *this;
    }
    const LogDistr& sub_zero()
    {
        *this -= m_logProbs[0];
        return *this;
    }
    /* cnop operation */
    const LogDistr& fft() // Multi-dimensional Fourier Transform
    {
        for (int i = 0; i < m_Q; ++i)
        {
            m_logProbs[i] = safe_exp(m_logProbs[i]);
            if (m_isnegative[i])
            {
                m_logProbs[i] = -m_logProbs[i];
            }
        }
        ntt(m_logProbs);
        for (int i = 0; i < m_Q; ++i)
        {
            m_isnegative[i] = false;
            if (m_logProbs[i] < 0)
            {
                m_isnegative[i] = true;
            }
            m_logProbs[i] = safe_log(ABS(m_logProbs[i]));
        }
        return *this;
    }
    void print()
    {
        for (int i = 0; i < m_Q; ++i)
        {
            mexPrintf("%e ", m_logProbs[i]);
        }
        mexPrintf("\n");
    }
    void printExp()
    {
        for (int i = 0; i < m_Q; ++i)
        {
            mexPrintf("%e ", m_isnegative[i]? -safe_exp(m_logProbs[i]) : safe_exp(m_logProbs[i]));
        }
        mexPrintf("\n");
    }
private:
    const double MAX_EXP;
    const double MIN_EXP;
private:
    double safe_exp(double x)
    {
        if (x > 700)
        {
            return MAX_EXP;
        }
        if (x < -700)
        {
            return MIN_EXP;
        }
        return exp(x);
    }
    double safe_log(double x)
    {
        if (x < MIN_EXP)
        {
            return -700;
        }
        return log(x);
    }
};

//=======================================================================================================================================
//===================================================== Helper functions ================================================================
//=======================================================================================================================================

void ntt(vector<double>& p) // Multi-dimensional Fourier Transform
{
    int Q = p.size();
    int K = (int) log2(Q);
    int factor = 1;
    for (int b = 0; b < K; b++)
    {
        for (int rest = 0; rest < Q/2; rest++)
        {
            int restH = rest >> b;
            int restL = rest & (factor-1);
            int rest0 = (restH << (b+1)) + restL;
            int rest1 = rest0 + factor;
            double prest0 = p[rest0];
            p[rest0] += p[rest1];
            p[rest1] = prest0 - p[rest1];
        }
        factor += factor;
    }
}

//=======================================================================================================================================
//===================================================== Joint SC list decoder ===========================================================
//=======================================================================================================================================


struct pcscl_list 
{
    vector<int> x;
    int idx;
};

vector<vector<float> > prob;
vector<vector<int> > u;
unsigned int idx; /* current bit index */
unsigned int L; /* list size */
vector<bool> used;
vector<tuple<float,int,int> > sorted;

// Joint Successive Cancellation algorithm
vector<pcscl_list> pcscl_jsc(const vector<vector<LogDistr> >& y, vector<vector<int> >::const_iterator f_it)
{
    const unsigned int N = y[0].size();
    const unsigned int L0 = y.size();
    const unsigned int K = y[0][0].m_K;
    const unsigned int Q = 1 << K;

    //mexPrintf("%d %d %d\n", N, L0, K);
    
    if (N == 1) 
    { // End recursion
        vector<pcscl_list> res;
        
        int f_mask = 0;
        int f_value = 0;
        const vector<int>& f_vec = *f_it;

/*
        for (int i = 0; i < K; ++i)
        {
            mexPrintf("%d ", f_vec[i]);
        }

        return res;
*/
        for (int i = 0; i < K; ++i)
        {
            if (f_vec[i] != 2)
            {
                f_mask += (1 << i);
                f_value += (f_vec[i] << i);
            }
        }

        //mexPrintf("%d %d\n", f_mask, f_value);    
        
        sorted.clear();
        sorted.reserve(Q*L0);

        for (int j = 0; j < Q; ++j)
        {
            if ((j & f_mask) != f_value)
            {
                // skip if frozen bits are incorrect
                continue;
            }
            for (int i = 0; i < L0; ++i)
            {
                sorted.emplace_back(prob[i][idx] + y[i][0].m_logProbs[j] /* log prob */, i /* path num */, j /* value */);
            }

        }

        if (L < sorted.size()) 
        {
            sort(sorted.begin(), sorted.end(), greater<decltype(sorted[0])>());
            sorted.resize(L);
        }
        res.resize(sorted.size());

        vector<int> last(L0,-1);
        for (unsigned i=0; i<sorted.size(); i++)
        {
            last[get<1>(sorted[i])] = i; // Find the last use
        }
        for (unsigned i=0; i<L0; i++)
        {
            if (last[i] < 0) used[i] = false;
        }
        for (unsigned i=0; i<sorted.size(); i++)
        {
            if (last[get<1>(sorted[i])] == (int)i) 
            {
                unsigned int j = get<1>(sorted[i]);
                u[j][idx] = get<2>(sorted[i]);
                prob[j][idx+1] = get<0>(sorted[i]);
                res[j].x.push_back(get<2>(sorted[i]));
                res[j].idx = j;
            } 
            else if (last[get<1>(sorted[i])] >= 0) 
            {
                unsigned int j=0;
                while (j<L && used[j]) j++; // Find empty space
                copy_n(u[get<1>(sorted[i])].begin(), idx, u[j].begin()); // CLONE PATH
                u[j][idx] = get<2>(sorted[i]);
                prob[j][idx+1] = get<0>(sorted[i]);
                res[j].x.push_back(get<2>(sorted[i]));
                res[j].idx = get<1>(sorted[i]);
                used[j] = true;
            }
        }
        

        ++idx;
        //mexPrintf("%d\n", idx);
        return res;
    } 
    else 
    {   // Recursion
        vector<vector<LogDistr> > u_est(L0, vector<LogDistr>(N/2, LogDistr(K)));
        LogDistr temp;
        for (int i=0; i<L0; i++)
        {
            for (int j=0; j<N/2; j++) 
            {
                // FFT-based convolution
                temp = y[i][2*j];
                u_est[i][j] = temp.fft();
                temp = y[i][2*j+1];
                u_est[i][j] += temp.fft();
                u_est[i][j].fft();
                //u_est[i][j].normalize();
                u_est[i][j] -= log(Q);
            }
        }

        vector<pcscl_list> res1 = pcscl_jsc(u_est, f_it);
        
        u_est.resize(res1.size(), vector<LogDistr>(N/2, LogDistr(K)));
        for (int i=0; i<res1.size(); i++) 
        {
            const vector<int>& u1hardprev = res1[i].x;
            const vector<LogDistr>& yy = y[res1[i].idx];
            for (int j=0; j<N/2; j++)
            {
                u_est[i][j] = (yy[2*j]+u1hardprev[j]) + yy[2*j+1];
                u_est[i][j].normalize();
            }
        }
        vector<pcscl_list> res2 = pcscl_jsc(u_est, f_it+N/2);
        for (int i=0; i<res2.size(); i++) 
        {
            vector<int> x(N);
            const vector<int> &x1 = res1[res2[i].idx].x;
            const vector<int> &x2 = res2[i].x;
            for (unsigned j=0; j<N/2; j++) 
            {
                x[2*j] = x1[j] ^ x2[j];
                x[2*j+1] = x2[j];
            }
            res2[i].x = move(x);
            res2[i].idx = res1[res2[i].idx].idx;
        }
        return res2;
    }
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Read arguments */
    
    int K = (int)*(mxGetPr(prhs[0])); // User count
    int m = (int)*(mxGetPr(prhs[1])); // Length
    int N = 1 << m;
    int Q = 1 << K;
    L = (int)*(mxGetPr(prhs[2])); // List size
    
    // f matrix
    double* input = mxGetPr(prhs[3]);
    
    vector<vector<int>> f_matrix(N, vector<int>(K, 0));
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < K; ++j)
        {
            f_matrix[i][j] = (int) input[i*K + j];
        }
    }
    

    // Log distributions
    input = mxGetPr(prhs[4]);
    
    vector<LogDistr> in_logDistrs(N, K);
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < Q; ++j)
        {
            in_logDistrs[i].m_logProbs[j] = input[i*Q + j];
        }
        in_logDistrs[i].normalize();
    }

    /* Initialization */
    used.assign(L,false);
    used[0] = true; // ?????
    u.assign(L, vector<int>(N));
    prob.assign(L, vector<float>(N+1, 0));
    idx = 0;
 
    /* Decode */
    pcscl_jsc({in_logDistrs}, f_matrix.begin());

    /* Return results to MATLAB */
    // iwd list
    plhs[0] = mxCreateDoubleMatrix(L, N, mxREAL);
    double* output = mxGetPr(plhs[0]);
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < L; ++j)
        {
            output[i*L + j] = u[j][i];
        }
    }


    // probabilities
    plhs[1] = mxCreateDoubleMatrix(L, N+1, mxREAL);
    output = mxGetPr(plhs[1]);
    for(int i = 0; i < N+1; ++i)
    {
        for(int j = 0; j < L; ++j)
        {
            output[i*L + j] = prob[j][i];
        }
    }
}









