/******************************************************************************
The contents of this source file are a lightly edited copy of parts of
six-subgraph (https://github.com/pkaski/six-subgraph), which is subject to 
the following license:

The MIT License (MIT)

Copyright (c) 2017-2019 P. Kaski

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

******************************************************************************/

#ifndef __IDIGM_POLY_H__
#define __IDIGM_POLY_H__

#include <memory>
#include "util.cpp"



/* Miscellanea ***************************************************************/

index_t ceilLog2(index_t n)
{
    assert (n > 0);
    index_t k {0};
    index_t nn {n};
    while (nn > 1) {
        nn = nn / 2;
        k += 1;
    }
    return k + (n > (1L << k) ? 1 : 0);
}

void randperm(index_t n, index_t *p)
{
    for(index_t i = 0; i < n; i++)
        p[i] = i;
    for(index_t j = 0; j < n-1; j++) {
        index_t u = j+rand()%(n-j);
        index_t t = p[u];
        p[u] = p[j];
        p[j] = t;
    }
    index_t *f = arrayAllocate<index_t>(n);
    for(index_t i = 0; i < n; i++)
        f[i] = 0;
    for(index_t i = 0; i < n; i++) {
        index_t u = p[i];
        assert(u >= 0 && u < n);
        f[u]++;
    }
    for(index_t i = 0; i < n; i++)
        assert(f[i] == 1);
    arrayDelete<index_t>(f);
}



/* Functions for working with arrays *****************************************/

template <typename T>
void arrayAdd(index_t na, const T* a, index_t nb, const T* b,
              index_t nc, T* c)
{
    //TODO: parallelize
    assert ( (nc >= na) && (nc >= nb) );
    if (na != nc || nb != nc) {
        for (index_t i = 0; i < nc; ++i)
            c[i] = (i < na ? a[i] : T::zero) + (i < nb ? b[i] : T::zero);
    }
    else {
        for (index_t i = 0; i < nc; ++i)
            c[i] = a[i] + b[i];
    }
}

template <typename T>
void arraySub(index_t na, const T* a, index_t nb, const T* b,
              index_t nc, T* c)
{
    // TODO: parallelize
    assert ( (nc >= na) && (nc >= nb) );
    if (na != nc || nb != nc) {
        for (index_t i = 0; i < nc; ++i)
            c[i] = (i < na ? a[i] : T::zero) - (i < nb ? b[i] : T::zero);
    }
    else {
        for (index_t i = 0; i < nc; ++i)
            c[i] = a[i] - b[i];
    }
}

template <typename T>
void arrayNeg(index_t n, const T* a, T* b)
{
    //TODO: parallelize
    for (index_t i = 0; i < n; ++i)
        b[i] = -a[i];
}

/* Gather N-length to n-length at 
 * offset, offset+1, ..., offset+n-1 over the D-bundle. */
template <typename F>
void arrayGather(index_t D, index_t N, index_t n, 
                 const F* in, F* out, index_t offset = 0)
{
    assert (N >= n);
    //#pragma omp parallel for if(D*n >= par_threshold)
    for(index_t v = 0; v < D*n; v++) {
        // in a CUDA kernel need to check v < D*n here 
        index_t vLow = v % n;
        index_t vHigh = v / n;
        index_t u = vHigh*N + vLow + offset;
        out[v] = in[u];
    }
}

/* Power-of-two version of the above. */
template <typename F>
void arrayGather2(index_t d, index_t N, index_t n, 
                  const F* in, F* out, index_t offset = -1)
{
    arrayGather(1L << d, 1L << N, 1L << n, in, out, 
                offset < 0 ? 0L : 1L << offset);
}

template <typename F>
void arrayScatter(index_t D, index_t n, index_t N, 
                  const F* in, F* out, index_t offset = 0)
{
    assert(N >= n);
    for(index_t v = 0; v < D*N; v++) {
        // in a CUDA kernel need to check v < D*N here 
        index_t vLow = v % N;
        index_t vHigh = v / N;
        index_t uLow = vLow - offset;       
        index_t u = vHigh * n + uLow;
        out[v] = (uLow >= 0 && uLow < n) ? in[u] : F::zero;
    }
}

/* Power-of-two version of the above. */
template <typename F>
void arrayScatter2(index_t d, index_t n, index_t N,
                   const F* in, F* out, index_t offset = -1)
{
    arrayScatter(1L << d, 1L << n, 1L << N, in, out, 
                 offset < 0 ? 0L : 1L << offset);
}

/* Assign ones to components of a power-of-two bundle. */
template <typename F>
void arrayMonic2(index_t d, index_t n, index_t l, F* f)
{
    assert(l < n && l >= 0);
    for (index_t i = 0; i < (1L << d); i++)
        f[i*(1L << n) + (1L << l)] = F::one;
}

/* Add one to each component of a power-of-two bundle. */
template <typename F>
void arrayAddOne2(index_t d, index_t n, F* f) 
{   
    for(index_t i = 0; i < (1L << d); i++)
        f[i*(1L << n)] = f[i*(1L << n)] + F::one;
}

template <typename F>
void arrayOne2(index_t d, index_t n, F* out)
{
//#pragma omp parallel for if((1L << (d+n)) >= par_threshold)
    for(index_t v = 0; v < (1L << (d+n)); v++)
        out[v] = ((v & ((1L << n)-1)) == 0) ? F::one : F::zero;
}

template <typename T>
void arrayScalarMul(index_t n, const T* in, T s, T* out)
{
    // TODO: parallelize
    for(index_t i = 0; i < n; ++i)
        out[i] = in[i] * s;
}

template <typename T>
void arrayZero(index_t n, T* f)
{
    // TODO: parallelize
    for(index_t i = 0; i < n; ++i)
        f[i] = T::zero;
}

template <typename F>
void arrayRand(index_t n, F* arr)
{
    for(index_t i = 0; i < n; i++)
        arr[i] = F::rand();
}

/* Reverse components of the N-length D-bundle. */
template <typename F>
void arrayRev(index_t D, index_t N, const F* f, F* g)
{
    for(index_t j = 0; j < D; j++)
        for(index_t i = 0; i < N; i++)
            g[j*N+N-1-i] = f[j*N+i];
}

template <typename T>
void arrayInterleave(index_t D, index_t N, const T* in0, const T* in1, T* out)
{
    //TODO: parallelize
    for (index_t v = 0; v < 2*D*N; ++v) {
        // in a CUDA kernel need to check v < 2*D*N here 
        index_t vLow = v % N;
        index_t vHigh = v / N;
        index_t vSelect = v % 2;
        vHigh = vHigh/2;        
        index_t u = vHigh*N + vLow;
        out[v] = vSelect == 0 ? in0[u] : in1[u];
    }
}

/* Subroutines for inverse modulo x^{L}. */

template <typename F>
void arrayAssignTo2(index_t d, index_t N, index_t k, 
                    const F* in, F* out,
                    const bool negate = false,
                    const F constantOffset = F::zero)
{
    // in  = a (d,N)-vector
    // out = a (d,k+1)-vector
    // out = (constantOffset*x^0 + (-1)^{negate}*in) mod x^{2^k} 
//#pragma omp parallel for if((1L << (d+k+1)) >= par_threshold)
    for(index_t v = 0; v < (1L << (d+k+1)); v++) {
        index_t vLow  = v&((1L << k)-1);
        index_t vPad = (v >> k)&1;
        index_t vHigh  = v >> (k+1);
        index_t cConstant = (vLow == 0) && (vPad == 0);
        index_t uHigh = vHigh * N;
        F s = cConstant ? constantOffset : F::zero;
        F t = (vPad == 0 && vLow < N) ? in[uHigh + vLow] : F::zero;
        out[v] = negate ? s - t : s + t;
    }
}

template <typename F>
void arrayAssignFrom2(index_t d, index_t k, index_t L, F* in, F* out)
{
    // in  = a (d,k+1)-vector
    // out = a (d,L)-vector
//#pragma omp parallel for if((1L << d)*L >= par_threshold)
    for(index_t v = 0; v < (1L << d)*L; v++) {
        // need to check v < (1L << d)*L here in a CUDA kernel
        index_t vLow = v % L;
        index_t vHigh = v / L;
        index_t u = vHigh*(1 << (k+1)) + vLow; 
        out[v] = in[u];
    }
}

// hostSSMul is defined further below
template <typename F>
void hostSSMul(index_t d, index_t n, const F* x, const F* y, F* z);

template <typename F>
void arrayInverseModXPowL(index_t d, index_t N, const F* f, index_t L, F* g)
{
    // f = a (d,N)-vector
    // g = a (d,L)-vector

    index_t l = ceilLog2(L);
    F* gg = arrayAllocate<F>(1L << (d+l+1));
    F* t  = arrayAllocate<F>(1L << (d+l+1));

    arrayOne2(d, 1, t); 
    index_t k = 0;
    while(k < l) {
        // Invariant: t is a (d,k+1)-array that holds the inverse mod x^{2^k}
        //            in the least 2^k coefficients of each component

        k++;
        // Now make the invariant hold for k, assuming it holds for k-1  ...
        //metric_push(m_quorem_layers_detail);
        arrayAssignTo2(d, 1L << k, k, t, gg); 
            // gg = t                    (as a (d,k+1)-array)
        arrayAssignTo2(d, N, k, f, t); // could-reverse-assign here
            // t = f mod x^{2^k}         (as a (d,k+1)-array)
        //ss_mul(d, k+1, t, gg, t); // could halve the length here to k
        hostSSMul(d, k+1, t, gg, t); // could halve the length here to k
            // t = fg                    (as a (d,k+1)-array)
        arrayAssignTo2(d, 1L << (k+1), k, t, t, 1, F::one + F::one);
            // t = (2-fg) mod x^{2^k}    (as a (d,k+1)-array)
            // Caveat: reading and writing to the same array here       
        //ss_mul(d, k+1, t, gg, gg); // need length k+1 due to cyclic wrap
        hostSSMul(d, k+1, t, gg, gg); // need length k+1 due to cyclic wrap
            // g = (2-fg)g               (as a (d,k+1)-array)
        arrayAssignTo2(d, 1L << (k+1), k, gg, t);
            // t = g mod x^{2^k}         (as a (d,k+1)-array)
        // The invariant now holds
        /*metric_pop(m_quorem_layers_detail,
                   "invlayer:  "
                   "d = %2ld,                              "
                   "l = %2ld, L = %10ld, N = %10ld, k = %2ld",
                   d, l, L, N, k);*/
    }
    arrayAssignFrom2(d, k, L, t, g);
    arrayDelete(t);
    arrayDelete(gg);
}



/* Schoenhage-Strassen multiplication ****************************************/

#define SIZE(a,b,c)         (1L << ((a)+(b)+(c)))
#define DIGIT2(u,a,b,c)     ((u) >> ((b)+(c)))
#define DIGIT1(u,a,b,c)     (((u) >> (c)) & ((1L << (b))-1))
#define DIGIT0(u,a,b,c)     ((u) & ((1L << (c))-1))
#define BUILD(x,y,z,a,b,c)  ((x) << ((b)+(c)))+((y) << (c))+(z)
#define MOD(v,n)            ((v) & ((1L << (n)) - 1))
#define NEG(j,s,n)          (MOD(MOD((j)-(s),n)+(s),(n)+1) >= (1L << (n)))

/* Expand from (d,t,m) to (d,t,m+1). */

template <typename F>
void ssExpand(index_t d, index_t t, index_t m, const F* in, F* out)
{   
    // metric_push(m_ss_layers);
    for(index_t v = 0; v < SIZE(d,t,m+1); v++) { 
        index_t p  = DIGIT2(v,d,t,m+1);
        index_t i  = DIGIT1(v,d,t,m+1);
        index_t j  = DIGIT0(v,d,t,m+1);
        index_t s  = i << (m+1-t);       
        index_t js = MOD(j-s, m+1);
        index_t u  = BUILD(p,i,js,d,t,m);        
        F t = (js < (1L << m)) ? in[u] : F::zero;
        out[v] = NEG(j, s, m+1) ? -t : t;
    }
    /*
    double time = metric_time(m_ss_layers);
    double trans_bytes = 3*sizeof(F)*(1L << (d+t+m));   
    metric_pop(m_ss_layers,
               "expand:   "
               "d = %2ld, t = %2ld, m = %2ld         (%6.2lfGiB/s)", 
               d, t, m, trans_bytes/((double)(1L << 30))/(time/1000.0));
    */
}


/* Base-case (d,m+1) multiply in F[x] / <x^{2M}+1>. */

template <typename F>
void ssMulBase(index_t d, index_t m, const F* x, const F* y, F* z)
{   
    // metric_push(m_ss_layers);
    for(index_t v = 0; v < (1L << (d+m+1)); ++v) {
        index_t L = 1L << (m+1);
        F s = F::zero;
        index_t uBase = v & ~(L-1);
        index_t k = v & (L-1);
        for(index_t a = 0; a < L; a++) {
            index_t u0 = uBase + a;
            index_t u1 = uBase + ((k-a)&(L-1));
            F t0 = x[u0];
            F t1 = y[u1];
            F p = t0*t1;
            index_t negate = (k < a) && (a <= k + L);
            if(negate)
                p = -p;
            s = s + p;
        }
        z[v] = s;
    }
    /*
    double time = metric_time(m_ss_layers);
    index_t mul_count = 1L << (d+2*(m+1));
    metric_pop(m_ss_layers, 
               "base:     "
               "d = %2ld, n = %2ld                (%6.2lfGmul/s)", 
                d, m+1, mul_count/1e9/(time/1000.0));
    */
}


/* Forward (d,t,m+1) butterfly at level w = t-1,...,1,0. */

template <typename F>
void ssButterflyForward(index_t d, index_t t, index_t m, index_t w, 
                        const F* in, F* out)
{
    //TODO: parallelize
    for(index_t v = 0; v < SIZE(d,t,m+1); v++) {
        index_t p  = DIGIT2(v,d,t,m+1);
        index_t i  = DIGIT1(v,d,t,m+1);
        index_t j  = DIGIT0(v,d,t,m+1);
        index_t iw = (i >> w)&1;
        index_t W  = 1L << w;
        index_t i0 = i & ~W;
        index_t i1 = i | W;
        index_t s  = (iw != 0) ? (i & (W-1)) << (m+1-w) : 0L; 
        index_t js = MOD(j-s,m+1);
        index_t u0 = BUILD(p,i0,js,d,t,m+1);
        index_t u1 = BUILD(p,i1,js,d,t,m+1); 
        F t0 = in[u0];      
        F t1 = in[u1];
        t0 = NEG(j, s, m+1) ? -t0 : t0;
        t1 = NEG(j, s, m+1) ? -t1 : t1;
        out[v] = (iw != 0) ? t0-t1 : t0+t1;
    }
}


/* Forward (d,t,m+1) FFT. Output is in bit-reversed order. */

template <typename F>
F* ssFFTForward(index_t d, index_t t, index_t m, 
                F* in, F* scratch)
{
    for(index_t w = t-1; w >= 0; w--) {
        //metric_push(m_ss_layers);
        ssButterflyForward(d, t, m, w, in, scratch);
        /*
        double time = metric_time(m_ss_layers);
        double trans_bytes = 3*sizeof(F)*(1L << (d+t+m+1));
        metric_pop(m_ss_layers,
                   "forward:  "
                   "d = %2ld, t = %2ld, m = %2ld, w = %2ld (%6.2lfGiB/s)", 
                    d, t, m, w, trans_bytes/((double)(1L << 30))/(time
                    1000.0));
        */
        F* temp = in;
        in = scratch;
        scratch = temp;
    }
    return in;
}


/* Inverse (d,t,m+1) butterfly at level w = 0,1,...,t-1. */

template <typename F>
void ssButterflyInverse(index_t d, index_t t, index_t m, index_t w, 
                        const F* in, F* out)
{
//#pragma omp parallel for if(SIZE(d,t,m+1) >= par_threshold)
    for(index_t v = 0; v < SIZE(d,t,m+1); v++) {
        index_t p  = DIGIT2(v,d,t,m+1);
        index_t i  = DIGIT1(v,d,t,m+1);
        index_t j  = DIGIT0(v,d,t,m+1);
        index_t iw = (i >> w)&1;
        index_t W  = 1L << w;       
        index_t i0 = i & ~W;
        index_t i1 = i | W;
        index_t s1 = -((i1 & (W-1)) << (m+1-w));
        index_t j0 = j;
        index_t j1 = MOD(j-s1,m+1);     
        index_t u0 = BUILD(p,i0,j0,d,t,m+1);
        index_t u1 = BUILD(p,i1,j1,d,t,m+1);        
        F t0 = in[u0];      
        F t1 = in[u1];
        t1 = NEG(j, s1, m+1) ? -t1 : t1;
        out[v] = (iw != 0) ? t0-t1 : t0+t1;
    }
}

/* Inverse (d,t,m+1) FFT. Assumes input is in bit-reversed order. */

template <typename F>
F* ssFFTInverse(index_t d, index_t t, index_t m, F* in, F* scratch)
{
    F z = F::inv_2_pow_k(t);
    arrayScalarMul(1L << (d+t+m+1), in, z, in);
    for(index_t w = 0; w < t; w++) {
        //metric_push(m_ss_layers);
        ssButterflyInverse(d, t, m, w, in, scratch);
        /*
        double time = metric_time(m_ss_layers);
        double trans_bytes = 3*sizeof(F)*(1L << (d+t+m+1));
        metric_pop(m_ss_layers,
                   "inverse:  "
                   "d = %2ld, t = %2ld, m = %2ld, w = %2ld (%6.2lfGiB/s)", 
                    d, t, m, w, trans_bytes/((double)(1L << 30))/(time
                    1000.0));
        */
        F* temp = in;
        in = scratch;
        scratch = temp;
    }
    return in;
}

/* Compress from (l,t,m+1) to (l,t,m). */

template <typename F>
void ssCompress(index_t d, index_t t, index_t m, const F* in, F* out)
{   
    //metric_push(m_ss_layers);
//#pragma omp parallel for if(SIZE(d,t,m) >= par_threshold)
    for(index_t v = 0; v < SIZE(d,t,m); v++) {
        index_t M    = 1L << m;
        index_t T    = 1L << t;
        index_t p    = DIGIT2(v,d,t,m);
        index_t i    = DIGIT1(v,d,t,m);
        index_t j    = DIGIT0(v,d,t,m);   
        index_t i1   = (i-1) & (T-1);
        index_t jM   = j+M;
        index_t s    = -(i << (m+1-t));      
        index_t s1   = -(i1 << (m+1-t));
        index_t js   = MOD(j-s,   m+1);
        index_t jMs1 = MOD(jM-s1, m+1);
        index_t u0   = BUILD(p, i,  js,   d, t, m+1);
        index_t u1   = BUILD(p, i1, jMs1, d, t, m+1);
        F t0 = in[u0];
        F t1 = in[u1];
        t0 = NEG(j,  s,  m+1) ? -t0 : t0;
        t1 = NEG(jM, s1, m+1) ? -t1 : t1;
        t1 = i1 > i ? -t1 : t1;     
        out[v] = t0 + t1;
    }
    /*
    double time = metric_time(m_ss_layers);
    double trans_bytes = 3*sizeof(F)*(1L << (d+t+m));
    metric_pop(m_ss_layers,
               "compress: "
               "d = %2ld, t = %2ld, m = %2ld         (%6.2lfGiB/s)", 
               d, t, m, trans_bytes/((double)(1L << 30))/(time/1000.0));
    */
}

/* Recursive (d,n),(d,n) -> (d,n) multiplication in F[x]/<x^N+1>. */

template <typename F>
void hostSSMul(index_t d, index_t n, const F* x, const F* y, F* z)
{
    assert(n >= 0);

    if(n <= 9) {
        bool haveBuf = false;
        F* out = z;
        if(x == z || y == z) {
            haveBuf = true;
            out = arrayAllocate<F>(1L << (d+n));
        }            
        if(n == 0)
            out[0] = x[0]*y[0];
        else
            ssMulBase(d, n-1, x, y, out);
        if(haveBuf) {
            arrayCopy(1L << (d+n), out, z);
            arrayDelete(out);
        }
        return;
    }

    index_t m = n/2;
    index_t t = n-m;
    assert(m <= t);

    F* masterBuf = arrayAllocate<F>(3*(1L << (d+t+m+1)));

    F* xe  = masterBuf;
    F* xes = masterBuf + (1L << (d+t+m+1));
    ssExpand(d, t, m, x, xe);
    F* xef = ssFFTForward(d, t, m, xe, xes);
    F* reclaimed;
    if(xef == xe)
        reclaimed = xes;
    else
        reclaimed = xe;

    F* ye  = masterBuf + 2*(1L << (d+t+m+1));
    F* yes = reclaimed;
    ssExpand(d, t, m, y, ye);
    F* yef = ssFFTForward(d, t, m, ye, yes);
    if(yef == ye)
        reclaimed = yes;
    else
        reclaimed = ye;

    F* zef = reclaimed;
    index_t calls = d+t;
    index_t block = 0;
    while(block + m + 1 < 26 && calls > 0) {
        block++;
        calls--;
    }

    F* xp = xef;
    F* yp = yef;
    F* zp = zef;
    for(index_t c = 0; c < (1L << calls); c++) {
        //ss_mul(block, m+1, xp, yp, zp);
        hostSSMul(block, m+1, xp, yp, zp);
        xp = xp + (1L << (block + m + 1));
        yp = yp + (1L << (block + m + 1));
        zp = zp + (1L << (block + m + 1));      
    }

    reclaimed = yef;
    // zef is also free at this point

    F* zes = reclaimed;
    F* ze = ssFFTInverse(d, t, m, zef, zes);

    ssCompress(d, t, m, ze, z);

    arrayDelete(masterBuf);

}



/* Quotient and remainder for polynomials ************************************/

/* Batch quotient and remainder for monic divisors. */

template <typename F>
void arrayQuoremMonic(index_t d, index_t n, const F* a,
                      index_t m, const F* b,
                      F* q, F* r)
{
    assert (n >= 1 && m >= 1);

    if (m > n) {
        for(index_t j = 0; j < (1L << d); j++)
            for(index_t i = 0; i < m; i++)
                r[j*m+i] = (i < n) ? a[j*n+i] : F::zero;
        // Caveat: should also set up a zero quotient.
        return;
    }

    /* Invariant: m <= n */
    index_t D = 1L << d;
    for(index_t i = 0; i < D; i++)
        assert(b[i*m+m-1] == F::one); // check that divisors are monic
    index_t k = 1+ceilLog2(n);
    index_t K = 1L << k;
    index_t l = 1+ceilLog2(n-m+1);
    index_t L = 1L << l;

    //metric_push(m_quorem);

    F* s = arrayAllocate<F>(1L << (d+k));
    F* t = arrayAllocate<F>(1L << (d+k));
    F* u = arrayAllocate<F>(1L << (d+k));

    // Reverse divisors

    //metric_push(m_quorem_layers);
    arrayRev(D, m, b, s); 
    /*
    metric_pop(m_quorem_layers,
               "reverse:   "
               "d = %2ld, k = %2ld, l = %2ld, "
               "K = %10ld, L = %10ld, "
               "n = %10ld, m = %10ld",
               d, k, l, K, L, n, m);*/
    
    // Compute truncated inverse of each reversed divisor (to u)

    //metric_push(m_quorem_layers);
    arrayInverseModXPowL(d, m, s, n-m+1, u);
    /*
    metric_pop(m_quorem_layers,
               "inverse:   "
               "d = %2ld, k = %2ld, l = %2ld, "
               "K = %10ld, L = %10ld, "
               "n = %10ld, m = %10ld",
               d, k, l, K, L, n, m);*/
    
    // Compute quotients (to q)

    //metric_push(m_quorem_layers);
    arrayScatter(D, n-m+1, L, u, s);
    arrayRev(D, n, a, t);
    arrayGather(D, n, n-m+1, t, u);
    arrayScatter(D, n-m+1, L, u, t);
    //ss_mul(d, l, s, t, u);
    hostSSMul(d, l, s, t, u);
    arrayGather(D, L, n-m+1, u, t);
    arrayRev(D, n-m+1, t, q);
    /*
    metric_pop(m_quorem_layers,
               "quotient:  "
               "d = %2ld, k = %2ld, l = %2ld, "
               "K = %10ld, L = %10ld, "
               "n = %10ld, m = %10ld",
               d, k, l, K, L, n, m);*/

    // Compute remainders (to r)

    //metric_push(m_quorem_layers);
    arrayScatter(D, n-m+1, K, q, t);
    arrayScatter(D, m, K, b, s);
    //ss_mul(d, k, s, t, s); // s = qb
    hostSSMul(d, k, s, t, s); // s = qb
    arrayGather(D, K, m-1, s, t);
    arrayGather(D, n, m-1, a, u);
    arraySub(D*(m-1), u, D*(m-1), t, D*(m-1), r);
    /*metric_pop(m_quorem_layers,
               "remainder: "
               "d = %2ld, k = %2ld, l = %2ld, "
               "K = %10ld, L = %10ld, "
               "n = %10ld, m = %10ld",
               d, k, l, K, L, n, m);

    metric_pop(m_quorem,
               "quorem:    "
               "d = %2ld, k = %2ld, l = %2ld, "
               "K = %10ld, L = %10ld, "
               "n = %10ld, m = %10ld",
               d, k, l, K, L, n, m);*/

    arrayDelete(u);
    arrayDelete(t);
    arrayDelete(s);
}

/* Non-monic quotient and remainder. */

template <typename F>
void arrayPolyQuorem(index_t n, const F* a, index_t m, const F* b, F* q, F* r)
{
    F t = b[m-1];
    F s = t.inv();

    F* bMonic = arrayAllocate<F>(m);
    arrayScalarMul<F>(m, b, s, bMonic);
    arrayQuoremMonic(0, n, a, m, bMonic, q, r);
    arrayScalarMul<F>(n-m+1, q, s, q);
    arrayDelete<F>(bMonic);
}


/* Miscellaneous functions for working with arrays representing polynomials **/

template <typename T>
index_t arrayPolyDeg(index_t n, const T* p)
{
    index_t deg;
    for (deg = n-1; deg >= 0; --deg) {
        if (p[deg] != T::zero)
            break;
    }
    return deg;  // NOTE: degree of an array of zeros is -1
}

template <typename T>
T arrayPolyEval(index_t n, const T* p, const T x)
{
    T y = p[n-1];
    for (index_t d = n-2; d >= 0; --d)
        y = y*x + p[d];
    return y;
}

template <typename T>
bool arrayPolyEq(index_t nf, const T* f,
                 index_t ng, const T* g)
{
    // if f is longer than g, swap them around
    if (nf > ng) {
        index_t ti = nf;
        nf = ng;
        ng = ti;
        const T* tp = f;
        f = g;
        g = tp;
    }
    // we can now assume g is longest, ng > nf
    for(index_t j = 0; j < ng; j++) {
        if (j >= nf) {
            if(g[j] != T::zero)
                return false;
        } else {
            if(f[j] != g[j])
                return false;
        }
    }
    return true;
}


/* Batch evaluation and interpolation with subproducts. **********************/

/* Batch evaluation at p points with the subproduct algorithm. */
template <typename F>
void arrayPolyBatchEval(index_t n, const F* f, 
                        index_t p, const F* u, 
                        F* evals) 
{
    assert(p >= 1);
    index_t k = ceilLog2(p); // Evaluate at 2^k points
    index_t K = 1L << k;
    
    /* Set up a subproduct tree with 2^{k+1}-1 nodes. */

    /* Level l = 0,1,...,k has 2^l nodes, each of degree 2^{k-l}.
     * Since the polynomials are monic, the tree will omit
     * the leading 1 from each polynomial. */

    F* a = arrayAllocate<F>(1L << (k+1));
    F* b = arrayAllocate<F>(1L << (k+1));
    F* s = arrayAllocate<F>(1L << (k+1));
    F* t = arrayAllocate<F>(1L << (k+1));
    F* q = arrayAllocate<F>(1L << (k+1));
    
    F* sub = arrayAllocate<F>((k+1)*K); 
      // (k+1)*2^k scalars, K = 2^k scalars at each level

    for(index_t i = 0; i < K; i++)
        sub[k*K+i] = i < p ? -u[i] : F::zero;

    /* Process the internal nodes one level at a time from bottom to top. */
    for(index_t l = k-1; l >= 0; l--) {
        
        /* For i = 0,...,1L << l, multiply nodes 2*i and 2*i+1 at level l+1. */

        arrayGather2 (l, k-l, k-(l+1), sub + (l+1)*K, q);
        arrayScatter2(l, k-(l+1), k-l, q, s);
        arrayGather2 (l, k-l, k-(l+1), sub + (l+1)*K, q, k-(l+1));
        arrayScatter2(l, k-(l+1), k-l, q, t);

        arrayMonic2(l, k-l, k-(l+1), s);
        arrayMonic2(l, k-l, k-(l+1), t);

        //ss_mul(l, k-l, s, t, sub + l*K); //  + i*(1L << (k-l))
        hostSSMul(l, k-l, s, t, sub + l*K); //  + i*(1L << (k-l))

        arrayAddOne2(l, k-l, sub + l*K); // cancel monic -1 wrap by adding 1
    }

    /* Descend down the subproduct tree and compute remainders. */

    /* Set up the root divisor & compute root remainder. */
    for(index_t j = 0; j < (1L << k); j++)
        b[j] = sub[0*K+j];
    b[(1L << k)] = F::one; // monic
    arrayQuoremMonic(0, n, f, (1L << k)+1, b, q, t); // note +1 here
    for(index_t i = n; i < (1L << k); i++)
        t[i] = F::zero;
       
    /* Now work down the levels. */
    F* r = t;
    F* w = s;
    for(index_t l = 1; l <= k; l++) {

        /* Invariant: r contains 2^{l-1} remainders of degree 2^{k-l+1}. */
        /* Compute next-level remainders to w. */

        arrayScatter(1L << l, 1L << (k-l), (1L << (k-l))+1, // note +1 here
                      sub + l*K, b);

        for(index_t j = 0; j < (1L << l); j++)
            b[((1L << (k-l))+1)*j + (1L << (k-l))] = F::one; // monic, note +1

        arrayInterleave(1L << (l-1), 1L << (k-l+1), r, r, a);

        arrayQuoremMonic(l, 1L << (k-l+1), a, (1L << (k-l))+1, b, q, w);

        /* Transpose r and w. */
        F* temp = w;
        w = r;
        r = temp;
    }
    assert(r != w);

    /* Copy result. */
    for(index_t i = 0; i < p; i++)
        evals[i] = r[i];

    arrayDelete(t);
    arrayDelete(s);
    arrayDelete(q);
    arrayDelete(b);
    arrayDelete(a);
    
    arrayDelete(sub);
}

/* Build a Lagrange polynomial for n points via subproducts. */
template <typename F>
void arrayLagrangeSubprod(index_t n, const F* u, const F* c, F* f)
{
    #if TRACK_TIME >= 4
    str_t nStr = std::to_string(n);
    logTimeStart("arrayLagrangeSubprod", nStr);
    #endif
    
    assert (n >= 1);
    index_t k = ceilLog2(n); // Evaluate at 2^k points
    index_t K = 1L << k;
    
    /* Ascend a subproduct tree with 2^{k+1}-1 nodes. */

    /* Level l = 0,1,...,k has 2^l nodes. */

    /* The full polynomials are monic of degree 2^{k-l},
     * the Lagrange polynomials are of degree 2^{k-l}-1. */

    F* t0 = arrayAllocate<F>(1L << (k+1));
    F* t1 = arrayAllocate<F>(1L << (k+1));
    F* t2 = arrayAllocate<F>(1L << (k+1));
    F* t3 = arrayAllocate<F>(1L << (k+1));
    F* x = arrayAllocate<F>(1L << (k+1));
    F* y = arrayAllocate<F>(1L << (k+1));
    F* z = arrayAllocate<F>(1L << (k+1));

    F* f0 = t0;
    F* f1 = t1;
    F* l0 = t2;
    F* l1 = t3;

    /* Prepare the leaf polynomials. */
    for(index_t i = 0; i < K; i++) {
        f0[i] = i < n ? -u[i] : F::zero;
           // the monic polynomial x-u[i]
        l0[i] = i < n ? c[i] : F::zero; 
           // coefficient for the Lagrange term that omits x-u[i]
    }

    /* Process the internal nodes one level at a time from bottom to top. */
    for(index_t l = k-1; l >= 0; l--) {

        // Multiply full nodes at level l+1 to get full node at level l. 

        arrayGather2 (l, k-l, k-(l+1), f0, z);
        arrayScatter2(l, k-(l+1), k-l, z, x);
        arrayGather2 (l, k-l, k-(l+1), f0, z, k-(l+1));
        arrayScatter2(l, k-(l+1), k-l, z, y);

        arrayMonic2(l, k-l, k-(l+1), x);
        arrayMonic2(l, k-l, k-(l+1), y);

        //ss_mul(l, k-l, x, y, f1);
        hostSSMul(l, k-l, x, y, f1);  

        arrayAddOne2(l, k-l, f1); // cancel monic -1 wrap by adding 1

        // Mix the full and Lagrange nodes at level l+1
        // to get Lagrange node at level l. 

        arrayGather2 (l, k-l, k-(l+1), f0, z);
        arrayScatter2(l, k-(l+1), k-l, z, x);
        arrayGather2 (l, k-l, k-(l+1), l0, z, k-(l+1));
        arrayScatter2(l, k-(l+1), k-l, z, y);

        arrayMonic2(l, k-l, k-(l+1), x);

        //ss_mul(l, k-l, x, y, l1);
        hostSSMul(l, k-l, x, y, l1);

        arrayGather2 (l, k-l, k-(l+1), l0, z);
        arrayScatter2(l, k-(l+1), k-l, z, x);
        arrayGather2 (l, k-l, k-(l+1), f0, z, k-(l+1));
        arrayScatter2(l, k-(l+1), k-l, z, y);

        arrayMonic2(l, k-l, k-(l+1), y);

        //ss_mul(l, k-l, x, y, z);
        hostSSMul(l, k-l, x, y, z);

        arrayAdd(1L << k, z, 1L << k, l1, 1L << k, l1);

        F* temp = f0;
        f0 = f1;
        f1 = temp;
        temp = l0;
        l0 = l1;
        l1 = temp;
    }
    for(index_t i = 0; i < n; i++)
        f[i] = l0[i+K-n];

    arrayDelete(z);
    arrayDelete(y);
    arrayDelete(x);
    arrayDelete(t3);
    arrayDelete(t2);
    arrayDelete(t1);
    arrayDelete(t0);
    
    #if TRACK_TIME >= 4
    logTimeEnd("arrayLagrangeSubprod", nStr);
    #endif
}

/* Interpolation with the subproduct algorithm. */
template <typename F>
void arrayPolyInterpolate(index_t n, const F* u, const F* y, F* f)
{
    #if TRACK_TIME >= 3
    str_t strN = std::to_string(n);
    logTimeStart("arrayPolyInterpolate", strN);
    #endif
    F* c = arrayAllocate<F>(n);
    for(index_t i = 0; i < n; i++)
        c[i] = F::one;
    arrayLagrangeSubprod(n, u, c, f);
    arrayPolyBatchEval(n, f, n, u, c);
    for(index_t i = 0; i < n; i++)
        c[i] = y[i]*c[i].inv();
    arrayLagrangeSubprod(n, u, c, f);   
    arrayDelete(c);
    #if TRACK_TIME >= 3
    logTimeEnd("arrayPolyInterpolate", strN);
    #endif
}



/* class for polynomials *****************************************************/

template <typename T>
class Poly {
private:
    bool init;
    index_t cap;
    std::shared_ptr<T> ptr;
    
    void resetToCap(index_t c) {
        assert ( c >= 0 );
        cap = c;
        ptr.reset(arrayAllocate<T>(c), ArrayDeleter<T>());
    }
    
    Poly(index_t c) : init {true}, cap {c}, ptr {nullptr} {
        resetToCap(c);
    }

public:
    
    Poly() : init {false}, cap {-1}, ptr{nullptr} {}
    
    Poly(const Poly& other) : init{other.init}, cap{-1}, ptr{nullptr} {
        if (init) {
            cap = other.cap;
            ptr = other.ptr;
        }
    }
    
    Poly(const T& val) : init {true}, cap {-1}, ptr {nullptr} {
        resetToCap(1);
        ptr.get()[0] = val;
    }
    
    Poly(const vec_t<T>& coeffs) : init {true}, cap {-1}, ptr {nullptr} {
        assert ( coeffs.size() < 1UL<<63 );
        index_t nCoeffs = (index_t)coeffs.size();
        resetToCap(nCoeffs);
        arrayCopy(nCoeffs, coeffs.data(), ptr.get());
    }
    
    bool initialized() const { return init; }
    
    index_t capacity() const {
        assert ( init );
        return cap;
    }
    
    index_t degree() const {
        assert ( init );
        return arrayPolyDeg(cap, ptr.get());
    }
    
    T& operator[] (index_t d) const {
        assert ( init && (d >= 0) && (d < cap ) );
        return ptr.get()[d];
    }
    
    T eval(const T x) const {
        return arrayPolyEval(cap, ptr.get(), x);
    }
    
    T operator() (const T x) const { return eval(x); }
    
    Poly& operator= (const Poly& other) {
        assert (other.init);
        init = true;
        cap = other.cap;
        ptr = other.ptr;
        return *this;
    }
    
    friend Poly operator+ (const Poly& a, const Poly& b) {
        assert (a.init && b.init);
        Poly r(a.cap > b.cap ? a.cap : b.cap);
        arrayAdd(a.cap, a.ptr.get(),
                 b.cap, b.ptr.get(),
                 r.cap, r.ptr.get());
        return r;
    }
    
    friend Poly operator- (const Poly& a) {
        assert (a.init);
        Poly b(a.cap);
        arrayNeg(a.cap, a.ptr.get(), b.ptr.get());
        return b;
    }
    
    friend Poly operator- (const Poly& a, const Poly& b) {
        assert (a.init && b.init);
        Poly r(a.cap > b.cap ? a.cap : b.cap);
        arraySub(a.cap, a.ptr.get(),
                 b.cap, b.ptr.get(),
                 r.cap, r.ptr.get());
        return r;
    }
    
    friend Poly operator* (const Poly& a, const Poly& b) {
        assert (a.init && b.init);
        index_t ad = a.degree();
        index_t bd = b.degree();
        index_t deg = ad > bd ? ad : bd;
        deg = deg <= 0 ? 1 : deg;
        index_t d = ceilLog2(deg + 1) + 1;
        index_t D = 1L << d;
        T* buf = arrayAllocate<T>(3*D);
        T* aa = buf;
        T* bb = buf + D;
        T* cc = buf + 2*D;
        arrayScatter(1, ad+1, D, a.ptr.get(), aa);
        arrayScatter(1, bd+1, D, b.ptr.get(), bb);
        hostSSMul(0, d, aa, bb, cc);
        deg = arrayPolyDeg(D, cc);
        if (deg < 0) { deg = 0; }
        Poly c(deg + 1);
        arrayCopy(deg + 1, cc, c.ptr.get());
        arrayDelete(buf);
        return c;
    }
    
    friend Poly operator* (const Poly& a, const T scalar) {
        assert (a.init);
        Poly b(a.cap);
        arrayScalarMul(a.cap, a.ptr.get(), scalar, b.ptr.get());
        return b;
    }
    
    friend Poly operator* (const T scalar, const Poly& a) {
        return a * scalar;
    }
    
    friend void quorem(const Poly& a, const Poly& b, Poly& q, Poly& r) {
        assert (a.init && b.init);
        index_t ad = a.degree();
        index_t bd = b.degree();
        assert (bd >= 0);  // cannot divide by zero
        index_t qd = ad - bd;
        if (qd < 0) {
            q = Poly(T::zero);
            r = a;
        }
        else {
            q = Poly(qd + 1);
            r = Poly(bd + 1);  // must be able to save zero remainder
            r.ptr.get()[bd] = T::zero;
            arrayPolyQuorem(ad + 1, a.ptr.get(),
                            bd + 1, b.ptr.get(),
                            q.ptr.get(),
                            r.ptr.get());
        }
    }
    
    friend Poly operator/ (const Poly& a, const Poly& b) {
        assert (a.init && b.init);
        Poly q, r;
        quorem(a, b, q, r);
        return q;
    }
    friend Poly operator% (const Poly& a, const Poly& b) {
        assert (a.init && b.init);
        Poly q, r;
        quorem(a, b, q, r);
        return r;
    }
    
    bool divides(const Poly& p) {
        return p % (*this) == Poly(T::zero);
    }
    
    friend bool operator== (const Poly& a, const Poly& b) {
        assert (a.init && b.init);
        return arrayPolyEq(a.cap, a.ptr.get(), b.cap, b.ptr.get());
    }
    
    friend bool operator!= (const Poly& a, const Poly& b) {
        return !(a == b);
    }
    
    static Poly x(index_t deg) {
        assert (deg >= 0);
        Poly<T> p(deg + 1);
        arrayZero(p.cap, p.ptr.get());
        p.ptr.get()[deg] = T::one;
        return p;
    }
    static Poly x() { return x(1); }
    
    static Poly interpolate(index_t n, const T* inputs, const T* targets) {
        assert (n >= 1);
        Poly<T> p(n);
        arrayPolyInterpolate(n, inputs, targets, p.ptr.get());
        return p;
    }
    
    void batchEval(index_t n, const T* inputs, T* evals) {
        assert (n >= 1);
        arrayPolyBatchEval(cap, ptr.get(), n, inputs, evals);
    }
    
    Poly take(index_t k) const {
        // Take k+1 largest-degree coefficients of f
        // and return them as a degree k polynomial;
        // return the zero polynomial if k is negative
        // (with degree capped at degree of polynomial)

        if (k < 0)
            return Poly(T::zero);
        
        index_t d = arrayPolyDeg(cap, ptr.get());
        if (k > d)
            k = d;
        Poly r(k+1);
        arrayCopy(k+1, ptr.get() + d - k, r.ptr.get());
        return r;
    }
    
    static Poly rand(index_t c) {
        assert (c >= 1);
        Poly a(c);
        arrayRand(a.cap, a.ptr.get());
        return a;
    }
    
    void copyCoeffsTo(T* dst) const {
        arrayCopy(cap, ptr.get(), dst);
    }
};


template <typename F>
std::ostream& operator<<(std::ostream& out, const Poly<F>& p)
{
    index_t deg = p.degree();
    if(deg < 0) {
        out << "0";
    } else {
        for(index_t i = deg; i >= 0; i--) {
            if(p[i] != F::zero) {
                if(i < deg)
                    out << " + ";
                if(p[i] != F::one || i == 0)
                    out << p[i];
                if(i > 0)
                    out << "x^{" << i << "}";
            }
        }
    }
    return out;
}



/********************************************* Extended Euclidean algorithm. */
/* Copied from six-subgraph */

template <typename F>
void gcdRecursive(index_t lvl,
                   const Poly<F> r_0,
                   const Poly<F> r_1,
                   index_t k,
                   index_t &h,
                   Poly<F> &s_h,
                   Poly<F> &t_h,
                   Poly<F> &s_hp1,
                   Poly<F> &t_hp1)
{
    // note: base case goes down to very low degree
    // should resort to a simpler algorithm earlier
    if(r_1.degree() < 0 || k < r_0.degree() - r_1.degree()) {
        h = 0;
        s_h = Poly<F>(F::one);
        t_h = Poly<F>(F::zero);
        s_hp1 = Poly<F>(F::zero);
        t_hp1 = Poly<F>(F::one);
        return;
    } else {
        if(k == 0 && k == r_0.degree() - r_1.degree()) {
            h = 1;
            s_h = Poly<F>(F::zero);
            t_h = Poly<F>(F::one);
            s_hp1 = Poly<F>(F::one);
            t_hp1 = Poly<F>(-r_0[r_0.degree()]*(r_1[r_1.degree()].inv()));
            return;
        }
    }
    index_t d = (k+1)/2;
    assert(d >= 1);

    Poly<F> f_s_h, f_t_h, f_s_hp1, f_t_hp1;
    index_t f_h;

    Poly<F> f_r_0 = r_0.take(2*d-2);
    Poly<F> f_r_1 = r_1.take(2*d-2-(r_0.degree()-r_1.degree()));
    
    gcdRecursive(lvl+1,
                  f_r_0,
                  f_r_1,
                  d-1,
                  f_h,
                  f_s_h,
                  f_t_h,
                  f_s_hp1,
                  f_t_hp1);

    index_t j = f_h + 1;
    index_t delta = f_t_hp1.degree();

    Poly<F> m_r_0 = r_0.take(2*k);
    Poly<F> m_r_1 = r_1.take(2*k-(r_0.degree()-r_1.degree()));

    // 4 x mul        
    Poly<F> r_jm1 = m_r_0*f_s_h   + m_r_1*f_t_h;
    Poly<F> r_j   = m_r_0*f_s_hp1 + m_r_1*f_t_hp1;

    if(r_j.degree() < 0 || 
       k < delta + r_jm1.degree() - r_j.degree()) {
        h = j-1;      
        s_h = f_s_h;
        t_h = f_t_h;
        s_hp1 = f_s_hp1;
        t_hp1 = f_t_hp1;
        return;
    }

    // quorem
    Poly<F> q_j, r_jp1;
    quorem(r_jm1, r_j, q_j, r_jp1); 
    
    index_t dstar = k - delta - (r_jm1.degree() - r_j.degree());
    assert(dstar >= 0);

    Poly<F> l_r_0 = r_j.take(2*dstar);
    Poly<F> l_r_1 = r_jp1.take(2*dstar-(r_j.degree()-r_jp1.degree()));
    Poly<F> l_s_h, l_t_h, l_s_hp1, l_t_hp1;
    index_t l_h;
    gcdRecursive(lvl+1,
                  l_r_0,
                  l_r_1,
                  dstar,
                  l_h,
                  l_s_h,
                  l_t_h,
                  l_s_hp1,
                  l_t_hp1);

    h = l_h + j;

    // 10 x mul 
    // (could use Strassen's here to get rid of at least one)
    Poly<F> u11 = f_s_hp1;
    Poly<F> u12 = f_t_hp1;
    Poly<F> u21 = f_s_h-q_j*f_s_hp1;
    Poly<F> u22 = f_t_h-q_j*f_t_hp1;
    s_h   = l_s_h*u11 + l_t_h*u21;
    t_h   = l_s_h*u12 + l_t_h*u22;
    s_hp1 = l_s_hp1*u11 + l_t_hp1*u21;
    t_hp1 = l_s_hp1*u12 + l_t_hp1*u22;   
}

template <typename F>
void gcd(const Poly<F> &r_0,
         const Poly<F> &r_1,
         index_t k,
         index_t &h,
         Poly<F> &r_h,
         Poly<F> &s_h,
         Poly<F> &t_h,
         Poly<F> &r_hp1,
         Poly<F> &s_hp1,
         Poly<F> &t_hp1)
{
    gcdRecursive(0, r_0, r_1, k, h, s_h, t_h, s_hp1, t_hp1);
    r_h   = r_0*s_h   + r_1*t_h;
    r_hp1 = r_0*s_hp1 + r_1*t_hp1;  

    assert(r_h   == r_0*s_h   + r_1*t_h);
    assert(r_hp1 == r_0*s_hp1 + r_1*t_hp1);
    if(r_hp1.degree() < 0) {
        assert(r_h.divides(r_0));
        assert(r_h.divides(r_1));
    } else {
        assert(r_0.degree() - r_h.degree() <= k);
        assert(r_0.degree() - r_hp1.degree() > k);
    }
}


/****************************************** Reed--Solomon encoding/decoding. */
/* Copied from six-subgraph */

template <typename F>
void rsEncode(index_t d, index_t e, F *p, F *q)
{
    assert(d >= 0 && d+1 <= e);
    Poly<F> f = Poly<F>::x(d);
    for(index_t i = 0; i <= d; i++)
        f[i] = p[i];
    F *u = arrayAllocate<F>(e);
    for(index_t i = 0; i < e; i++)
        u[i] = F(i);
    f.batchEval(e, u, q);
    arrayDelete(u);
}

template <typename F>
void corrupt(index_t d, index_t e, F *f)
{
    assert(d >= 0 && d+1 <= e);
    index_t n = (e-d-1)/2;
    index_t *q = arrayAllocate<index_t>(e);
    randperm(e, q);
    for(index_t i = 0; i < n; i++)
        f[q[i]] = F::rand();
    arrayDelete(q);
}

template <typename F>
bool rsDecodeXY(index_t d, index_t e, F *x, F *y, Poly<F> &p)
{
    #if TRACK_TIME
    logTimeStart("rsDecodeXY", std::to_string(d), true);
    #endif
    
    // Gao's decoder
    assert(d >= 0 && d+1 <= e);

    F *u = arrayAllocate<F>(e+1);
    F *c = arrayAllocate<F>(e+1);
    F *f = arrayAllocate<F>(e+1);
    Poly<F> g0, g1, r0, s0, t0, r1, s1, t1;

    for(index_t i = 0; i < e+1; i++) {
        u[i] = (i < e) ? x[i] : F::zero;
        c[i] = (i == e) ? F::one : F::zero;
    }
    arrayLagrangeSubprod(e+1, u, c, f);
    g0 = Poly<F>::x(e);
    for(index_t i = 0; i < e; i++)
        g0[i] = f[i];

    arrayDelete(f);
    arrayDelete(u);
    arrayDelete(c);

    g1 = Poly<F>::interpolate(e, x, y);

    index_t k = (e-d-1)/2;  
    index_t h;

    gcd(g0, g1, k, h, r0, s0, t0, r1, s1, t1);

    bool success = r1.degree() < (e+d+1)/2;
    if(success) {
        success = (t1.degree() <= r1.degree()) && 
                  (r1.degree() - t1.degree() <= d);
        if(success) {
            quorem(r1, t1, r0, s0);
            success = s0.degree() < 0;
            if(success)
                p = r0;
        }
    }
    
    #if TRACK_TIME
    logTimeEnd("rsDecodeXY", std::to_string(d), true);
    #endif
    return success;
}

#endif  // __IDIGM_POLY_H__


