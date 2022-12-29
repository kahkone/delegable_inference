#ifndef __IDIGM_CRT_UTIL_H__
#define __IDIGM_CRT_UTIL_H__

#include <tuple>
#include <gmp.h>
#include <cmath>
#include <limits>
#include "montgomery.cpp"
#include "util.cpp"
#include "factor.cpp"
#include "factorgraph.cpp"
#include "primes_list.h"


class MPInt {
private:
    mpz_t* ptr;
    bool init;
    
    static MPInt performBinaryOp(void (*op)(mpz_t, const mpz_t, const mpz_t),
                                 const MPInt& l, const MPInt& r)
    {
        assert (l.init && r.init);
        mpz_t result; mpz_init(result);
        op(result, *l.ptr, *r.ptr);
        MPInt mpIntResult(result);
        mpz_clear(result);
        return mpIntResult;
    }

public:
    
    static const MPInt zero, one;
    
    MPInt(): ptr{nullptr}, init{false} {}
    MPInt(mpz_t w): ptr{new mpz_t[1]}, init{true} {
        mpz_init(*ptr);
        mpz_set(*ptr, w);
    };
    MPInt(long int x): ptr{new mpz_t[1]}, init{true} {
        mpz_init(*ptr);
        mpz_set_si(*ptr, x);
    };
    MPInt(const MPInt& other): ptr{new mpz_t[1]}, init{true} {
        assert (other.init);
        mpz_init(*ptr);
        mpz_set(*ptr, *other.ptr);
    };
    ~MPInt() {
        if (init) {
            mpz_clear(*ptr);
            if (ptr != nullptr) delete[] ptr;
        }
    }
    
    MPInt& operator=(const MPInt& other) {
        assert (other.init);
        if (!init) {
            ptr = new mpz_t[1];
            mpz_init(*ptr);
            init = true;
        }
        mpz_set(*ptr, *other.ptr);
        return *this;
    }
    
    MPInt operator-() const {
        assert (init);
        mpz_t result; mpz_init(result);
        mpz_neg(result, *ptr);
        MPInt mpIntResult(result);
        mpz_clear(result);
        return mpIntResult;
    }
    
    MPInt& operator+=(const MPInt& other) {
        assert (init && other.init);
        mpz_add(*ptr, *ptr, *other.ptr);
        return *this;
    }
    MPInt& operator-=(const MPInt& other) {
        assert (init && other.init);
        mpz_sub(*ptr, *ptr, *other.ptr);
        return *this;
    }
    MPInt& operator*=(const MPInt& other) {
        assert (init && other.init);
        mpz_mul(*ptr, *ptr, *other.ptr);
        return *this;
    }
    
    friend MPInt operator+(const MPInt& l, const MPInt& r) {
        return MPInt::performBinaryOp(&mpz_add, *l.ptr, *r.ptr);
    }
    friend MPInt operator-(const MPInt& l, const MPInt& r) {
        return MPInt::performBinaryOp(&mpz_sub, *l.ptr, *r.ptr);
    }
    friend MPInt operator*(const MPInt& l, const MPInt& r) {
        return MPInt::performBinaryOp(&mpz_mul, *l.ptr, *r.ptr);
    }
    friend MPInt operator/(const MPInt& l, const MPInt& r) {
        assert ( mpz_cmp_si(*r.ptr, 0) != 0 );
        return MPInt::performBinaryOp(&mpz_tdiv_q, *l.ptr, *r.ptr);
    }
    friend MPInt operator%(const MPInt& l, const MPInt& r) {
        assert ( mpz_cmp_si(*r.ptr, 0) != 0 );
        return MPInt::performBinaryOp(&mpz_tdiv_r, *l.ptr, *r.ptr);
    }
    
    friend MPInt mod(const MPInt& l, const MPInt& r) {
        assert ( mpz_cmp_si(*r.ptr, 0) != 0 );
        return MPInt::performBinaryOp(&mpz_mod, *l.ptr, *r.ptr);
    }
    
    // modular inverse of x mod *ptr
    MPInt modInverse(const MPInt& x) const {
        assert ( x.init );
        assert ( mpz_cmp_si(*x.ptr, 0) != 0 );
        assert ( mpz_cmp_si(*ptr, 0) != 0 );
        mpz_t result; mpz_init(result);
        int inverseExists = mpz_invert(result, *x.ptr, *ptr);
        assert ( inverseExists );
        MPInt mpIntResult(result);
        mpz_clear(result);
        return mpIntResult;
    }
    
    friend bool operator==(const MPInt& l, const MPInt& r) {
        assert ( l.init && r.init );
        return (mpz_cmp(*l.ptr, *r.ptr) == 0);
    }
    friend bool operator>=(const MPInt& l, const MPInt& r) {
        assert ( l.init && r.init );
        return (mpz_cmp(*l.ptr, *r.ptr) >= 0);
    }
    friend bool operator<=(const MPInt& l, const MPInt& r) {
        assert ( l.init && r.init );
        return (mpz_cmp(*l.ptr, *r.ptr) <= 0);
    }
    friend bool operator<(const MPInt& l, const MPInt& r) {
        assert ( l.init && r.init );
        return (mpz_cmp(*l.ptr, *r.ptr) < 0);
    }
    friend bool operator>(const MPInt& l, const MPInt& r) {
        assert ( l.init && r.init );
        return (mpz_cmp(*l.ptr, *r.ptr) > 0);
    }
    
    void copyToMPZ(mpz_t& m) const {
        mpz_init(m);
        mpz_set(m, *ptr);
    }
    
    long int toLongInt() const {
        // abort if *ptr does not fit into a long int
        assert (mpz_cmpabs_ui(*ptr, 0xFFFFFFFFFFFFFFFF) < 0);
        return mpz_get_si(*ptr);
    }
    int toInt() const {
        // abort if *ptr does not fit into an int
        assert (mpz_cmpabs_ui(*ptr, 0xFFFFFFFF) < 0);
        return (int)mpz_get_si(*ptr);
    }
    
    static MPInt rand() { return MPInt((uint_t)std::rand()); }
    static MPInt randInRange0To(MPInt upper) {
        MPInt result = MPInt::rand();
        while (result < upper)
            result *= MPInt::rand();
        return result % (upper + 1);
    }
    static MPInt randInRange(MPInt lower, MPInt upper) {
        return randInRange0To(upper - lower) + lower;
    }
    
    static str_t to_string(const MPInt& x) {
        assert (x.init);
        char* arr = new char[mpz_sizeinbase(*x.ptr, 10)+2];
        mpz_get_str(arr, 10, *x.ptr);
        str_t result(arr);
        delete[] arr;
        return result;
    }
};

const MPInt MPInt::zero = MPInt(0L);
const MPInt MPInt::one = MPInt(1L);

std::ostream& operator<<(std::ostream& out, const MPInt& x)
{
    mpz_t m; mpz_init(m);
    x.copyToMPZ(m);
    out << m;
    mpz_clear(m);
    return out;
}


int getFloatSign(float f) {
    int x = *((int*)(&f));
    return (x & (1<<31)) != 0;
}

// Return the 8-bit exponent field of f, uninterpreted, as an integer
int getFloatExpField(float f)
{
    int x = *((int*)(&f));
    return (x >> 23) & 0xFF;
}

// Return exponent of non-NaN, non-Inf float f, interpreted as per IEEE 754
int getFloatExp(float f)
{
    constexpr int expWidth {8};
    int rawExp = getFloatExpField(f);
    assert (rawExp < (1 << expWidth) - 1);  // f must be neither NaN nor Inf
    constexpr int expBias {127};
    int exp = rawExp >= 1 ? rawExp - expBias : -126;
    return exp;
}

int getFloatSignificand(float f)
{
    int e {getFloatExpField(f)};
    int x = *((int*)(&f));
    if (e >= 1)
        return (1<<23) + (x & 0x7FFFFF);
    return 0 + (x & 0x7FFFFF);
}

int getLSBExp(int x)
{
    if (x == 0)
        return 0;
    int lsbPos;
    for (lsbPos = 0; lsbPos < 32; ++lsbPos) {
        if (x & (1 << lsbPos))
            break;
    }
    return lsbPos;
}

void floatToMPZ(float f, mpz_t& result, int expOf2=149)
{
    assert ( std::isfinite(f) );
    mpz_init(result);
    if (f == 0.0) {
        mpz_set_si(result, 0);
        return;
    }
    /* IEEE 754-2019 specifies that
     *    emax = 127.
     *    emin = 1 - emax = -126.
     *    p = 24.
     *    The exponent bias b = 127.
     *    For T != 0, the 23-bit trailing significand T is interpreted
     *    as being multiplied by 2^(1-p) = 2^(-23).
     * By IEEE 754-2019 section 3.4, the smallest possible unsigned float is
     *    2^(-126)*(0+2^(-23)) = 2^(-149).
     * Thus: Multiplying any float by 2^(149) gives an integer.
     * For floats with larger least significant bit, expOf2 may be
     * chosen to be smaller than 149.
     */
    assert (expOf2 < 150);
    int exp = getFloatExp(f);
    // getFloatSignificand already returns an integer, i.e. multiplies the
    // significand by 2^(23)
    int scalingExp = exp + expOf2 - 23;
    int significand = getFloatSignificand(f);
    int finalLSBExp = exp + expOf2 - (23 - getLSBExp(significand));
    assert (finalLSBExp >= 0);
    
    mpz_set_si(result, significand);
    if (scalingExp >= 0)
        mpz_mul_2exp(result, result, (uint_t)scalingExp);
    else
        mpz_cdiv_q_2exp(result, result, (uint_t)(-scalingExp));
    
    if (getFloatSign(f))
        mpz_neg(result, result);
}

MPInt floatToMPInt(float f, int expOf2=149)
{
    mpz_t rawResult;
    floatToMPZ(f, rawResult, expOf2);
    MPInt result(rawResult);
    mpz_clear(rawResult);
    return result;
}

Factor<MPInt> floatToMPIntFactor(const Factor<float>& f, int expOf2=149)
{
    long_uint_t vol = f.volume;
    MPInt* intData = new MPInt[vol];
    for (long_uint_t k = 0; k < vol; ++k) {
        intData[k] = floatToMPInt(f.valueAt(k), expOf2);
    }
    Factor<MPInt> fZ(f.label, f.vars, intData);
    delete[] intData;
    return fZ;
}

// Return exponent corresponding to least significant bit in IEEE-754 float f.
// If (f == 0), return ifZeroReturn.
int getLSBExp(float f, int ifZeroReturn=0)
{
    int exp = getFloatExp(f);
    int significand = getFloatSignificand(f);
    int lsbPos = getLSBExp(significand);
    assert (lsbPos < 24);
    int rawExp = getFloatExpField(f); assert (rawExp >= 0);
    if (significand == 0 && rawExp == 0)
        return ifZeroReturn;
    return exp + (lsbPos - 23);
}

/* Return the smallest possible exponent k such that
 * for all values v in f, (v * 2^k) is an integer. */
int intConversionExp(const Factor<float>& f,
                     bool nonnegative=true, int ifZeroReturn=0)
{
    int startExp {999};
    int minLSBExp {startExp};
    // the smallest possible exponent of any bit in a IEEE 754 float is -149
    int zeroExp {-150};
    for (long_uint_t p = 0; p < f.volume; ++p) {
        int lsbExp = getLSBExp(f.valueAt(p), zeroExp);
        if (lsbExp == zeroExp)
            continue;
        if (lsbExp < minLSBExp)
            minLSBExp = lsbExp;
    }
    if (minLSBExp == startExp)  // all values in f were 0.0
        return ifZeroReturn;
    int conversionExp = -minLSBExp;
    if (nonnegative)
        return conversionExp >= 0 ? conversionExp : 0;
    return conversionExp;
}

/* Return the smallest possible exponent k such that
 * for all values v in f, v <= 2^k. */
int upperBoundExp(const Factor<float>& f)
{
    int ubExp = -150;
    for (long_uint_t k = 0; k < f.getVolume(); ++k) {
        int exp = getFloatExp(f.valueAt(k)) + 1;
        if (exp > ubExp)
            ubExp = exp;
    }
    return ubExp;
}

/* conversionExps: map in which to store the exponents of 2 used in converting
 * each Factor from float to MPInt. */
FactorGraph<MPInt> floatToMPIntFG(const FactorGraph<float>& fg,
                                  map_t<label_t,int>& conversionExps,
                                  map_t<label_t,int>& upperBoundExps,
                                  bool nonnegative=true)
{
    assert ( conversionExps.size() == 0 );
    assert ( upperBoundExps.size() == 0 );
    FactorGraph<MPInt> fgZ;
    int zeroExp {0};
    for (Factor<float> f : fg.getFactors()) {
        int cExp = intConversionExp(f, nonnegative, zeroExp);
        conversionExps[f.label] = cExp;
        fgZ.insertFactor(floatToMPIntFactor(f, cExp));
        upperBoundExps[f.label] = upperBoundExp(f);
    }
    fgZ.setBoundary(fg.getBoundary());
    return fgZ;
}

float mpIntToFloat(MPInt x, int conversionExp, bool allowSubnormOrInf=false)
{
    mpf_t mpfX; mpf_init(mpfX);
    mpz_t mpzX;
    x.copyToMPZ(mpzX);
    mpf_set_z(mpfX, mpzX);
    if (conversionExp > 0)
        mpf_mul_2exp(mpfX, mpfX, (uint_t)conversionExp);
    else if (conversionExp < 0)
        mpf_div_2exp(mpfX, mpfX, (uint_t)(-conversionExp));

    mpf_t maxFloat; mpf_init(maxFloat);
    mpf_set_d(maxFloat, std::numeric_limits<float>::max());
    mpf_t minFloat; mpf_init(minFloat);
    mpf_set_d(minFloat, std::numeric_limits<float>::denorm_min());
    mpf_t lowestFloat; mpf_init(lowestFloat);
    mpf_set_d(lowestFloat, std::numeric_limits<float>::lowest());
    
    mpf_t mpfXAbs; mpf_init(mpfXAbs);
    mpf_abs(mpfXAbs, mpfX);
    bool underflow = (mpf_cmp_d(mpfX, 0.0) != 0)
                     && (mpf_cmp(mpfXAbs, minFloat) < 0);
    bool overflow = (mpf_cmp(mpfX, lowestFloat) < 0)
                    || (mpf_cmp(mpfX, maxFloat) > 0);
    if ( ! allowSubnormOrInf) {
        assert ( ! underflow );
        assert ( ! overflow );
    }
    // TODO: implement logging, log warnings if (underflow || overflow)
    
    float result = (float)mpf_get_d(mpfX);
    mpz_clear(mpzX);
    mpf_clears(mpfX, maxFloat, minFloat, lowestFloat, mpfXAbs, NULL);
    return result;
}


// Convert an MPInt to Montgomery form.
template <typename P>
mont_t mpIntToMontgomery(const MPInt x)
{
    MPInt xModN = x % MPInt(P::modulus);
    if (xModN < 0L)
        xModN += P::modulus;
    #if DEBUG_ASSERTS
    // xModN can safely be cast to uint_t
    assert ( xModN >= 0L );
    assert ( xModN < 1L<<32 );
    #endif
    return to_montgomery<P>((uint_t)xModN.toLongInt());
}

template <typename P>
Zp<P> mpIntToMod(const MPInt x)
{
    return Zp<P>(mpIntToMontgomery<P>(x), true);
}

template <typename P>
Factor<Zp<P>> mpIntToModFactor(const Factor<MPInt>& fIn)
{
    long_uint_t vol = fIn.volume;
    Zp<P>* dataOut = new Zp<P>[vol];
    for (long_uint_t k = 0; k < vol; ++k)
        dataOut[k] = mpIntToMod<P>(fIn.valueAt(k));
    Factor<Zp<P>> fOut(fIn.label, fIn.vars, dataOut);
    delete[] dataOut;
    return fOut;
}

Factor<float> mpIntToFloatFactor(
    const Factor<MPInt>& fIn,
    int conversionExp,
    bool allowSubnormOrInf=false)
{
    long_uint_t vol = fIn.volume;
    float* dataOut = new float[vol];
    for (long_uint_t k = 0; k < vol; ++k)
        dataOut[k] = mpIntToFloat(
                        fIn.valueAt(k), conversionExp, allowSubnormOrInf);
    Factor<float> fOut(fIn.label, fIn.vars, dataOut);
    delete[] dataOut;
    return fOut;
}

template <typename F>
Factor<uint_t> modToUintFactor(const Factor<F>& fIn)
{
    long_uint_t vol = fIn.volume;
    uint_t* dataOut = new uint_t[vol];
    for (long_uint_t k = 0; k < vol; ++k)
        dataOut[k] = (fIn.valueAt(k)).value();
    Factor<uint_t> fOut(fIn.label, fIn.vars, dataOut);
    delete[] dataOut;
    return fOut;
}

template <typename F>
Factor<F> uintToModFactor(const Factor<uint_t>& fIn)
{
    long_uint_t vol = fIn.volume;
    F* dataOut = new F[vol];
    for (long_uint_t k = 0; k < vol; ++k)
        dataOut[k] = F(fIn.valueAt(k));
    Factor<F> fOut(fIn.label, fIn.vars, dataOut);
    delete[] dataOut;
    return fOut;
}

Factor<MPInt> uintToMPIntFactor(const Factor<uint_t>& fIn)
{
    long_uint_t vol = fIn.volume;
    MPInt* dataOut = new MPInt[vol];
    for (long_uint_t k = 0; k < vol; ++k)
        dataOut[k] = MPInt(fIn.valueAt(k));
    Factor<MPInt> fOut(fIn.label, fIn.vars, dataOut);
    delete[] dataOut;
    return fOut;
}

template <typename P>
FactorGraph<Zp<P>> mpIntToModFG(const FactorGraph<MPInt>& fgIn)
{
    FactorGraph<Zp<P>> fgOut;
    for (Factor<MPInt> f : fgIn.getFactors())
        fgOut.insertFactor(mpIntToModFactor<P>(f));
    fgOut.setBoundary(fgIn.getBoundary());
    return fgOut;
}

/* For all i, gs[i] is expected to take values in {0, ..., uintModuli[i]-1}.
 * Let N be the product of all uintModuli.
 * crtReconstruct first computes the Factor f with values f::v in {0, ..., N-1}
 * such that for all i, f::v % moduli[i] == gs[i]::v.
 * If recoverNegatives is true, then
 *     reconstructed values v satisfying 0 <= v < N/2 will remain positive,
 *     and values v satisfying N/2 <= v < N will be reinterpreted as (v-N).
 */
Factor<MPInt> crtReconstruct(vec_t<Factor<uint_t>> gs,
                             vec_t<uint_t> uintModuli,
                             label_t resultId=0,
                             bool recoverNegatives=false)
{
    uint_t k = uintModuli.size();
    #if DEBUG_ASSERTS
    assert ( k > 0 );
    assert ( gs.size() == k );
    #endif
    
    vec_t<MPInt> moduli(k);
    for (uint_t i = 0; i < k; ++i)
        moduli[i] = MPInt(uintModuli[i]);
    
    // moduliProds[i] = \prod_{j!=i} moduli[j]
    vec_t<MPInt> moduliProds(k);
    MPInt allModuliProd = vectorProduct(moduli);
    for (uint_t i = 0; i < k; ++i)
        moduliProds[i] = allModuliProd / moduli[i];

    // invModuliProds[i] = \prod_{j!=i} (inverse of moduli[j] mod moduli[i])
    vec_t<MPInt> invModuliProds(k);
    for (uint_t i = 0; i < k; ++i) {
        invModuliProds[i] = MPInt(1);
        for (uint_t j = 0; j < k; ++j) {
            if (j == i)
                continue;
            invModuliProds[i] *= moduli[i].modInverse(moduli[j]);
        }
    }
    
    // permute Factors in gs to be in same order
    vec_t<Var> orderedVars = gs[0].vars;
    vec_t<Factor<MPInt>> hs;
    for (Factor<uint_t> g : gs)
        hs.push_back(uintToMPIntFactor(permuteFactor<uint_t>(orderedVars, g)));
    
    long_uint_t vol = hs[0].getVolume();
    #if DEBUG_ASSERTS
    for (uint_t i = 0; i < k; ++i)
        assert ( hs[i].getVolume() == vol );
    #endif
    
    MPInt* reconstructedData = new MPInt[vol];
    for (long_scalar_t p = 0; p < vol; ++p) {
        MPInt prodSum(0L);
        for (uint_t i = 0; i < k; ++i) {
            prodSum += hs[i].valueAt(p) * invModuliProds[i] * moduliProds[i];
        }
        prodSum = prodSum % allModuliProd;
        if (recoverNegatives) {
            if (prodSum >= allModuliProd / 2)
                prodSum -= allModuliProd;
        }
        reconstructedData[p] = prodSum;
    }
    
    Factor<MPInt> result(resultId, orderedVars, reconstructedData);
    delete[] reconstructedData;
    return result;
}

//NOTE: IDIGM_PRIMES defined in primes_list.h (generated by batch-gen.py)
vec_t<uint_t> primesToReach(MPInt x, vec_t<uint_t> primesList=IDIGM_PRIMES)
{
    MPInt y(1);
    vec_t<uint_t> primes;
    for (uint_t p : primesList) {
        if (y >= x)
            break;
        y *= MPInt(p);
        primes.push_back(p);
    }
    assert ( y >= x );
    return primes;
}



#endif  // __IDIGM_CRT_UTIL_H__

