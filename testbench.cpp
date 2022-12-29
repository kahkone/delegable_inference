#ifndef __IDIGM_TESTBENCH_H__
#define __IDIGM_TESTBENCH_H__

#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "util.cpp"
#include "factorgraph.cpp"
#include "eval_algo.cpp"


template <typename T>
class Array2D {
private:
    long_uint_t nRows;
    long_uint_t nCols;
    T* data;
public:
    Array2D(): nRows {0}, nCols {0}, data {nullptr} {};
    Array2D(long_uint_t r, long_uint_t c, const T* dataSrc):
        nRows {r}, nCols {c}, data {new T[r*c]}
    {
        long_uint_t size = r*c;
        assert ( size < 1UL<<63 );
        arrayCopy((index_t)size, dataSrc, data);
    }
    Array2D(const Array2D& other):
        nRows {other.nRows}, nCols {other.nCols}, data {new T[nRows*nCols]}
    {
        long_uint_t size = nRows*nCols;
        assert ( size < 1UL<<63 );
        arrayCopy((index_t)size, other.data, data);
    }
    ~Array2D()
    {
        if (data != nullptr)
            delete[] data;
    }
    
    Array2D& operator= (const Array2D& other)
    {
        if (data != nullptr)
            delete[] data;
        nRows = other.nRows;
        nCols = other.nCols;
        long_uint_t size = nRows * nCols;
        assert ( size < 1UL<<63 );
        data = new T[size];
        arrayCopy((index_t)size, other.data, data);
        return *this;
    }
    
    bool operator== (const Array2D& other) const
    {
        if (nRows != other.nRows || nCols != other.nCols)
            return false;
        for (uint_t r = 0; r < nRows; ++r) {
            for (uint_t c = 0; c < nCols; ++c) {
                if (valueAt(r,c) != other.valueAt(r,c))
                    return false;
            }
        }
        return true;
    }
    bool operator!= (const Array2D& other) const
    {
        return ( ! (*this == other));
    }
    
    // read/write access to single element
    T& operator() (long_uint_t r, long_uint_t c)
    {
        assert ( r < nRows );
        assert ( c < nCols );
        return data[r*nCols + c];
    }
    // read-only access to single element
    T valueAt(long_uint_t r, long_uint_t c) const
    {
        assert ( r < nRows );
        assert ( c < nCols );
        return data[r*nCols + c];
    }
    // read-only access to all data
    const T* getData() const { return data; }
    
    long_uint_t getNRows() const { return nRows; }
    long_uint_t getNCols() const { return nCols; }
    
    static Array2D rand(long_uint_t r, long_uint_t c, T min=0, T max=0)
    {
        assert (min <= max);
        long_uint_t size {r*c};
        vec_t<T> v(size);
        if (min != max) {
            for (long_uint_t i = 0; i < size; ++i)
                v[i] = genRandScalar(min, max);
        }
        else {
            for (long_uint_t i = 0; i < size; ++i)
                v[i] = genRandScalar<T>();
        }
        return Array2D(r, c, v.data());
    }
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const Array2D<T>& a)
{
    out << "[\n";
    for (long_uint_t r = 0; r < a.getNRows(); ++r) {
        out << "  [";
        for (long_uint_t c = 0; c < a.getNCols(); ++c) {
            out << a.valueAt(r,c);
            if (c != a.getNCols() - 1)
                out << ", ";
        }
        out << "]\n";
    }
    out << "]\n";
    return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const vec_t<Array2D<T>>& v)
{
    out << "[\n";
    long_uint_t s = v.size();
    for (long_uint_t n = 0; n < s; ++n) {
        out << v[n];
        if (n < s-1)
            out << ",\n";
    }
    out << "]\n";
    return out;
}


template <typename T>
T randInRange(T min, T max)
{
    T result = min;
    result += (T)rand() % (1 + max - min);
    return result;
}


template <typename I, typename Params>
class InstanceGenerator {
private:
    Params params;
    I (*generator)(Params);
public:
    InstanceGenerator(Params p, I (*g)(Params)): params {p}, generator {g} {}
    I operator()() const { return generator(params); }
    Params getParameters() const { return params; }
};


template <typename T>
struct ParamsMatMul {
    uint_t minNOperands;
    uint_t maxNOperands;
    uint_t minSize;
    uint_t maxSize;
    T minValue;
    T maxValue;
    
    ParamsMatMul(uint_t minNOps, uint_t maxNOps,
                 uint_t minS, uint_t maxS,
                 T minV, T maxV):
        minNOperands {minNOps}, maxNOperands {maxNOps},
        minSize {minS}, maxSize {maxS},
        minValue {minV}, maxValue {maxV}
    {
        assert (minNOps > 0);
        assert (minNOps <= maxNOps);
        assert (minS > 0);
        assert (minS <= maxS);
        assert (minV <= maxV);
    }
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const ParamsMatMul<T>& p)
{
    out << "ParamsMatMul {"
        << p.minNOperands << ", " << p.maxNOperands << ", "
        << p.minSize << ", " << p.maxSize << ", "
        << p.minValue << ", " << p.maxValue << "}";
    return out;
}


template <typename T>
vec_t<Array2D<T>> instanceGenMatMul(ParamsMatMul<T> params)
{
    uint_t nOps = randInRange(params.minNOperands, params.maxNOperands);
    vec_t<Array2D<T>> operands;
    uint_t r = randInRange(params.minSize, params.maxSize);
    uint_t c = randInRange(params.minSize, params.maxSize);
    for (uint_t i = 0; i < nOps; ++i) {
        operands.push_back(
            Array2D<T>::rand(r, c, params.minValue, params.maxValue));
        r = c;
        c = randInRange(params.minSize, params.maxSize);
    }
    return operands;
}

// NOTE: No checks for numeric overflow or underflow
template <typename T>
Array2D<T> solveMatMulBase(vec_t<Array2D<T>> operands)
{
    #if TRACK_TIME
    logTimeStart("solveMatMulBase", "", true);
    #endif
    
    uint_t n = operands.size();
    assert ( n > 0 );
    if (n == 1)
        return operands[0];
    
    Array2D<T> result = operands[0];
    for (uint_t i = 1; i < n; ++i) {
        uint_t ra = result.getNRows();
        uint_t ca = result.getNCols();
        const T* dataA = result.getData();
        
        uint_t rb = operands[i].getNRows();
        uint_t cb = operands[i].getNCols();
        const T* dataB = operands[i].getData();
        assert ( rb == ca );
        T* permutedDataB = new T[cb*rb];
        permuteArray<T>(2, {1,0}, {rb,cb}, dataB, permutedDataB);
        
        T* resultData = new T[ra*cb];
        naiveMatrixMul(ra, ca, cb, dataA, permutedDataB, resultData);
        result = Array2D(ra, cb, resultData);
        
        delete[] permutedDataB;
        delete[] resultData;
    }
    
    #if TRACK_TIME
    logTimeEnd("solveMatMulBase", "", true);
    logTimeTrackingSummaries();
    #endif
    
    return result;
}

// NOTE: Assumption: F is Zp<P> for some Montgomery params P
template <typename F>
Array2D<F> solveMatMulDFGC(vec_t<Array2D<F>> operands)
{
    #if TRACK_TIME
    logTimeStart("solveMatMulDFGC", "", true);
    logTimeStart("solveMatMulDFGC:conversionToFG");
    #endif
    
    uint_t n = operands.size();
    assert ( n >= 2 );
    
    // Translate operands (matrices) into Factors; construct a FactorGraph
    FactorGraph<F> fg;
    vec_t<label_t> factorIds;
    for (uint_t i = 0; i < n; ++i) {
        Var r(i, operands[i].getNRows());
        Var c(i+1, operands[i].getNCols());
        Factor<F> f(fg.getFreeFactorId(), {r,c}, operands[i].getData());
        fg.insertFactor(f);
        factorIds.push_back(f.getId());
    }
    fg.setBoundary({
        Var(0, operands[0].getNRows()),
        Var(n, operands[n-1].getNCols())
    });
    vec_t<Var> cutset = setToVector( fg.getBoundary() );
    vec_t<IdTriple> contractionOrder { 
        {factorIds[0], factorIds[1], fg.getFreeFactorId()}
    };
    for (uint_t i = 1; i < n-1; ++i) {
        contractionOrder.push_back({
            contractionOrder[i-1].newLabel,
            factorIds[i+1],
            fg.getFreeFactorId()
        });
    }
    #if TRACK_TIME
    logTimeStart("solveMatMulDFGC:etasTauElls");
    #endif
    map_t<Var, map_t<uint_t, F>> etas = constructEtas<F>(cutset);
    map_t<vec_t<scalar_t>, F> tau = constructTau<F>(cutset);
    map_t<Var, Poly<F>> ells = interpolateElls(cutset, etas, tau);
    #if TRACK_TIME
    logTimeEnd("solveMatMulDFGC:etasTauElls");
    logTimeEnd("solveMatMulDFGC:conversionToFG");
    #endif
    
    // Solve the inference problem over the FactorGraph
    Poly<F> proofPoly;
    Factor<F> factorResult = contractViaProof(
        fg, cutset, etas, tau, ells, proofPoly, contractionOrder);
    
    // Translate back to Array2D
    vec_t<Var> vars = factorResult.getVars();
    assert ( vars.size() == 2 );
    if (vars[0].label != 0) {
        factorResult = permuteFactor({vars[1],vars[0]}, factorResult);
        vars = factorResult.getVars();
    }
    uint_t nr = vars[0].size;
    uint_t nc = vars[1].size;
    #if DEBUG_ASSERTS
    assert ( nr == operands[0].getNRows() );
    assert ( nc == operands[n-1].getNCols() );
    #endif
    Array2D<F> result(nr, nc, factorResult.dataAt(0));
    
    #if TRACK_TIME
    logTimeEnd("solveMatMulDFGC", "", true);
    logTimeTrackingSummaries();
    clearTimeTrackingData();
    #endif
    
    // Verify the proof polynomial
    assert ( verifyProofPolynomial(proofPoly, fg, cutset, etas, tau) );
    #if TRACK_TIME
    logTimeTrackingSummaries();
    #endif
    
    return result;
}


struct ParamsPermanent {
    uint_t minN;
    uint_t maxN;
    
    ParamsPermanent(uint_t minSize, uint_t maxSize):
        minN {minSize}, maxN {maxSize}
    {
        assert (minN > 0);
        assert (minN <= maxN);
    }
};

std::ostream& operator<<(std::ostream& out, const ParamsPermanent& p)
{
    out << "ParamsPermanent {" << p.minN << ", " << p.maxN << "}";
    return out;
}


template <typename T>
Array2D<T> instanceGenPermanent(ParamsPermanent params)
{
    uint_t n = randInRange(params.minN, params.maxN);
    return Array2D<T>::rand(n,n);
}


// Compute the permanent of an n-by-n matrix using Ryser's formula
template <typename T>
T solvePermanentBase(Array2D<T> matrix)
{
    #if TRACK_TIME
    logTimeStart("solvePermanentBase", "", true);
    #endif
    uint_t n = matrix.getNCols(); 
    assert ( n == matrix.getNRows() );
    assert ( n < 31 );
    vec_t<uint_t> elems;
    for (uint_t e = 0; e < n; ++e) { elems.push_back(e); }
    index_t nSubsets = 1L << n;
    T permanent {0};
    for (index_t s = 0; s < nSubsets; ++s) {
        vec_t<uint_t> set = subset(elems, s);
        T subproduct {1};
        for (uint_t i = 0; i < n; ++i) {
            T subsum {0};
            for (uint_t j : set)
                subsum += matrix.valueAt(i, j);
            subproduct *= subsum;
        }
        long sign = pow(-1L, set.size());
        // NOTE: writing  "permanent += subproduct * pow(-1L, set.size());" 
        // would lead to incorrect results when T=Zp<_>; hence the if-else
        if (sign < 0)
            permanent -= subproduct;
        else
            permanent += subproduct;
    }
    long sign = pow(-1, n);
    if (sign < 0)
        permanent = -permanent;
    
    #if TRACK_TIME
    logTimeEnd("solvePermanentBase", "", true);
    logTimeTrackingSummaries();
    #endif
    return permanent;
}

// NOTE: exponential complexity
template <typename T>
T solvePermanentDFGC(Array2D<T> matrix)
{
    #if TRACK_TIME
    logTimeStart("solvePermanentDFGC", "", true);
    logTimeStart("solvePermanentDFGC:conversionToFG");
    #endif
    uint_t n = matrix.getNCols(); 
    assert ( n == matrix.getNRows() );
    assert ( n > 0 );
    assert ( n < 26 );  // exponential memory complexity
    
    // Translate the matrix into a FactorGraph
    FactorGraph<T> fg;
    vec_t<Var> cutset = {};
    vec_t<Var> vars;
    for (uint_t i = 0; i < n; ++i) {
        Var v(i,2);
        vars.push_back(v);
        if (i % 2 == 0)
            cutset.push_back(v);
    }
    vec_t<Factor<T>> factors;
    for (uint_t k = 0; k < n; ++k) {
        long nSubsets = pow(2,n);
        assert ( nSubsets < 1L<<62 );
        vec_t<T> data((long_uint_t)nSubsets);
        for (long w = 0; w < nSubsets; ++w) {
            T sum = T(0);
            for (uint_t pos = 0; pos < n; ++pos) {
                if ( (w >> pos) & 1 )
                    sum += matrix.valueAt(k,pos);
            }
            if ( ! ((w >> k) & 1) )
                sum = -sum;
            data[(long_uint_t)w] = sum;
        }
        Factor<T> f(fg.getFreeFactorId(), vars, data.data());
        fg.insertFactor(f);
    }
    
    #if TRACK_TIME
    logTimeStart("solvePermanentDFGC:etasTauElls");
    #endif
    map_t<Var, map_t<uint_t, T>> etas = constructEtas<T>(cutset);
    map_t<vec_t<scalar_t>, T> tau = constructTau<T>(cutset);
    map_t<Var, Poly<T>> ells = interpolateElls(cutset, etas, tau);
    #if TRACK_TIME
    logTimeEnd("solvePermanentDFGC:etasTauElls");
    logTimeEnd("solvePermanentDFGC:conversionToFG");
    #endif
    
    // Solve the inference problem over the FactorGraph
    Poly<T> proofPoly;
    Factor<T> factorResult = contractViaProof(
        fg, cutset, etas, tau, ells, proofPoly);
    
    // Recover the scalar result
    assert ( factorResult.isScalar() );
    T result = factorResult.scalarValue();
    
    #if TRACK_TIME
    logTimeEnd("solvePermanentDFGC", "", true);
    logTimeTrackingSummaries();
    clearTimeTrackingData();
    #endif
    
    // Verify the proof polynomial
    assert ( verifyProofPolynomial(proofPoly, fg, cutset, etas, tau) );
    #if TRACK_TIME
    logTimeTrackingSummaries();
    #endif
    
    return result;
}

// NOTE: O(n * n!) complexity
template <typename T>
T solvePermanentDFGC_slow(Array2D<T> matrix)
{
    #if TRACK_TIME
    logTimeStart("solvePermanentDFGC_slow");
    #endif
    
    
    uint_t n = matrix.getNCols(); 
    assert ( n == matrix.getNRows() );
    assert ( n > 0 );
    assert ( n < 13 );  // Assuming we have less than 81 Gb of memory
    
    // Translate the matrix into a FactorGraph
    FactorGraph<T> fg;
    uint_t nPermutations = factorial(n);
    Var v(0, nPermutations);
    vec_t<uint_t> perm;
    for (uint_t i = 0; i < n; ++i)
        perm.push_back(i);
    for (uint_t k = 0; k < n; ++k) {
        std::sort(perm.begin(), perm.end());
        vec_t<T> kProjections;
        do {
            kProjections.push_back(matrix.valueAt(k, perm[k]));
        } while (std::next_permutation(perm.begin(), perm.end()));
        Factor<T> f(fg.getFreeFactorId(), {v}, kProjections.data());
        fg.insertFactor(f);
    }
    vec_t<Var> cutset = {};
    
    #if TRACK_TIME
    logTimeStart("solvePermanentDFGC_slow:etasTauElls");
    #endif
    map_t<Var, map_t<uint_t, T>> etas = constructEtas<T>(cutset);
    map_t<vec_t<scalar_t>, T> tau = constructTau<T>(cutset);
    map_t<Var, Poly<T>> ells = interpolateElls(cutset, etas, tau);
    #if TRACK_TIME
    logTimeEnd("solvePermanentDFGC_slow:etasTauElls");
    logTimeEnd("solvePermanentDFGC_slow:conversionToFG");
    #endif
    
    // Solve the inference problem over the FactorGraph
    Poly<T> proofPoly;
    Factor<T> factorResult = contractViaProof(
        fg, cutset, etas, tau, ells, proofPoly);
    
    // Recover the scalar result
    assert ( factorResult.isScalar() );
    T result = factorResult.scalarValue();
    
    #if TRACK_TIME
    logTimeEnd("solvePermanentDFGC_slow");
    logTimeTrackingSummaries();
    clearTimeTrackingData();
    #endif
    
    // Verify the proof polynomial
    assert ( verifyProofPolynomial(proofPoly, fg, cutset, etas, tau) );
    #if TRACK_TIME
    logTimeTrackingSummaries();
    #endif
    
    return result;
}


struct TBOptions {
    uint_t numberOfRuns;
    int verbosity;
};

template <typename I, typename P, typename S>
int runExperiments(
    vec_t<InstanceGenerator<I, P>> instanceGenerators,
    vec_t<S (*)(I)> solutionAlgos,
    TBOptions opts)
{
    uint_t nAlgos = solutionAlgos.size();
    assert ( nAlgos > 0 );
    uint_t nRuns = opts.numberOfRuns;
    assert ( nRuns > 0 );
    uint_t nInstances = instanceGenerators.size();

    int returnStatus = 0;
    
    // Generate instances, compute and compare solutions, record metrics
    vec_t<vec_t<S>> solutions(nAlgos);
    for (uint_t i = 0; i < nInstances; ++i) {
        auto generator = instanceGenerators[i];
        if (opts.verbosity > 0)
            std::cout << "Beginning runs with InstanceParams "
                      << generator.getParameters() << "\n";
        for (uint_t r = 0; r < nRuns; ++r) {
            if (opts.verbosity > 0)
                std::cout << "Beginning run " << r << "\n";
            I instance = generator();
            if (opts.verbosity > 1)
                std::cout << "Instance: " << instance << "\n";
            for (uint_t a = 0; a < nAlgos; ++a) {
                if (opts.verbosity > 0)
                    std::cout << "Executing algorithm " << a << "\n";
                solutions[a].push_back(solutionAlgos[a](instance));
            }
            if (opts.verbosity > 0)
                std::cout << "Checking equality of algorithms' outputs\n";
            // NOTE: Assuming operator== is transitive on type S
            S solution0 = solutions[0][r];
            for (uint_t a = 1; a < nAlgos; ++a) {
                if (solutions[a][r] != solution0) {
                    std::cout << "\nFAILURE: algorithms 0 and " << a
                        << " disagree on input " << instance <<"\n";
                    std::cout << "Algorithm 0 produced output "
                        << solution0 << "\nbut algorithm " << a 
                        << "produced output " << solutions[a][r] << "\n";
                    returnStatus = 1;
                }
            }
            if (opts.verbosity > 0)
                std::cout << "Run " << r << " complete.\n";
        }
    }
    
    return returnStatus;
}


#endif  // __IDIGM_TESTBENCH_H__
