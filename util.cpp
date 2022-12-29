#ifndef __IDIGM_UTIL_H__
#define __IDIGM_UTIL_H__

#include <set>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <exception>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <sys/time.h>
#include <ctime>
#include <chrono>
#include "montgomery.cpp"


typedef unsigned int scalar_t;
typedef unsigned int uint_t;
typedef long unsigned int long_scalar_t;
typedef long unsigned int long_uint_t;
typedef long int index_t;

template <typename T>
using vec_t = std::vector<T>;

template <typename T>
using set_t = std::set<T>;

template <typename K, typename V>
using map_t = std::map<K,V>;

using str_t = std::string;

str_t timestamp() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    long int microsecs = 1000000 * tv.tv_sec + tv.tv_usec;
    return std::to_string(microsecs);
}



/* Functions for working with set_t ******************************************/

template <typename T>
set_t<T> vectorToSet(vec_t<T> v)
{
    set_t<T> s;
    s.insert(v.begin(), v.end());
    return s;
}

template <typename T, typename Predicate>
set_t<T> filterSet(Predicate p, const set_t<T>& inputSet)
{
    set_t<T> filteredSet;
    for (T x : inputSet)
        if (p(x)) filteredSet.insert(x);
    return filteredSet;
}

template <typename T>
bool isIn(set_t<T> s, T x)
{
    auto it = s.find(x);
    return ( it != s.end() );
}

template <typename T>
set_t<T> complementSet(const set_t<T>& s1, const set_t<T>& s2)
{
    set_t<T> r;
    for (auto it = s1.begin(); it != s1.end(); ++it) {
        if ( !(isIn(s2, *it)) )
            r.insert(*it);
    }
    return r;
}

template <typename T>
set_t<T> intersectSets(const set_t<T>& s1, const set_t<T>& s2)
{
    set_t<T> r;
    for (auto it = s1.begin(); it != s1.end(); ++it) {
        if (isIn(s2, *it)) r.insert(*it);
    }
    return r;
}

template <typename T>
set_t<T> unionSets(set_t<T> s1, const set_t<T>& s2)
{
    for (auto it = s2.begin(); it != s2.end(); ++it) {
        if (!(isIn(s1, *it))) s1.insert(*it);
    }
    return s1;
}

template <typename T>
void printSet(const set_t<T>& s)
{
    std::cout << "{";
    for (T e : s)
        std::cout << e << ", ";
    std::cout << "}\n";
}

template <typename T>
vec_t<T> setToVector(const set_t<T>& s)
{
    vec_t<T> v;
    v.insert(v.begin(), s.begin(), s.end());
    return v;
}

template <typename T>
bool isSuperset(const set_t<T>& super, const set_t<T>& sub)
{
    for (T x : sub) {
        bool superContainsX {false};
        for (T y : super) {
            if (x == y) {
                superContainsX = true;
                break;
            }
        }
        if (!superContainsX)
            return false;
    }
    return true;
}


template <typename T>
vec_t<T> randSubset(vec_t<T>, uint_t);

template <typename T>
bool isIn(const vec_t<T>&, T);

template <typename T>
set_t<T> randSubset(const set_t<T>& s, uint_t n)
{
    assert ( n <= s.size() );
    vec_t<uint_t> allInds(n);
    for (uint_t i = 0; i < n; ++i)
        allInds[i] = i;
    vec_t<uint_t> randInds = randSubset<uint_t>(allInds, n);
    set_t<T> subset;
    uint_t ind = 0;
    for (auto it = s.begin(); it != s.end(); ++it) {
        if (isIn(randInds, ind))
            subset.insert(*it);
        ++ind;
    }
    return subset;
}

template <typename T>
bool equal(const set_t<T>& s1, const set_t<T>& s2)
{
    return isSuperset(s1, s2) && isSuperset(s2, s1);
}


/* Functions for working with vec_t ************************************/


// Fold a vector with (*)
template <typename T>
T vectorProduct(const vec_t<T> vec)
{
    // NOTE: Returning 1 for empty vectors is expected
    T p {1};
    for (T x : vec) p *= x;
    return p;
}

template <typename T>
bool isIn(const vec_t<T>& v, T x)
{
    bool found {false};
    for (T e : v) {
        if (e == x) {
            found = true;
            break;
        }
    }
    return found;
}

// Return v1 - v2
template <typename T>
vec_t<T> complement(vec_t<T> v1, const vec_t<T>& v2)
{
    typename vec_t<T>::iterator it;
    for (it = v1.begin(); it != v1.end(); ) {
        if (isIn(v2, *it))
            it = v1.erase(it);
        else
            ++it;
    }
    return v1;
}

template <typename T>
vec_t<T> intersect(const vec_t<T>& v1, const vec_t<T>& v2)
{
    vec_t<T> r;
    for (auto it = v1.begin(); it != v1.end(); ++it) {
        if (isIn(v2, *it)) r.push_back(*it);
    }
    return r;
}

template <typename T>
vec_t<T> unionVectors(vec_t<T> v1, const vec_t<T>& v2)
{
    for (auto it = v2.begin(); it != v2.end(); ++it) {
        if (!(isIn(v1, *it))) v1.push_back(*it);
    }
    return v1;
}

template <typename T>
str_t vecToStr(const vec_t<T>& v)
{
    std::stringstream s; s << "[";
    for (uint_t i = 0; i < v.size(); ++i) {
        s << v[i];
        if (i < v.size() - 1)
            s << ", ";
    }
    s << "]";
    return s.str();
}

template <typename T>
void printVector(const vec_t<T>& v)
{
    std::cout << vecToStr(v) << "\n";
}

template <typename T>
index_t indexOf(T x, const vec_t<T>& vec)
{
    long unsigned int n = vec.size();
    if (n > 1UL<<62) {
        throw std::invalid_argument(
            "Expected vector of size less than 1UL<<62.\n");
    }
    for (long unsigned int j = 0; j < n; ++j)
        if (x == vec[j]) return (index_t)j;
    return -1;
}

template <typename T>
bool isSuperset(vec_t<T> super, vec_t<T> sub)
{
    for (T x : sub) {
        bool superContainsX {false};
        for (T y : super) {
            if (x == y) {
                superContainsX = true;
                break;
            }
        }
        if (!superContainsX)
            return false;
    }
    return true;
}

template <typename T>
vec_t<T> moveToEnd(const vec_t<T>& host,
                         const vec_t<T>& tail)
{
    #if DEBUG_ASSERTS
    assert ( isSuperset(host, tail) );
    #endif
    
    // NOTE: duplicates are deleted from host
    vec_t<T> r = complement(host, tail);
    r.insert(r.end(), tail.begin(), tail.end());
    return r;
}

template <typename T>
vec_t<T> moveToFront(const vec_t<T>& host,
                           const vec_t<T>& head)
{
    #if DEBUG_ASSERTS
    assert ( isSuperset(host, head) );
    #endif
    
    // NOTE: duplicates are deleted from host
    vec_t<T> r = complement(host, head);
    r.insert(r.begin(), head.begin(), head.end());
    return r;
}

template <typename T>
vec_t<scalar_t> getPermutation(const vec_t<T>& v1,
                                     const vec_t<T>& v2)
{
    #if DEBUG_ASSERTS
    assert ( isSuperset(v1, v2) );
    assert ( isSuperset(v2, v1) );
    assert ( v1.size() == v2.size() );
    #endif
    
    vec_t<scalar_t> perm (v1.size(), 0);
    for (scalar_t i = 0; i < v2.size(); ++i) {
        perm[i] = (scalar_t) indexOf(v2[i], v1);
    }
    return perm;
}

template <typename T>
vec_t<T> randSubset(vec_t<T> vec, scalar_t n)
{
    assert ( n <= vec.size() );
    vec_t<T> subset;
    for (scalar_t i = 0; i < n; ++i) {
        scalar_t j = (scalar_t)rand() % vec.size();
        subset.push_back(vec[j]);
        vec.erase(vec.begin() + j);
    }
    return subset;
}

/* Interpret i as a binary string, each bit corresponding to an element of s;
 * return the subset of s encoded by that binary string.
 * NOTE: Only works for sets with at most 64 elements. */
template <typename T>
vec_t<T> subset(const vec_t<T>& s, index_t i)
{
    assert ( s.size() <= 64 );
    vec_t<T> subset;
    for (T x : s) {
        if (i & 1)
            subset.push_back(x);
        i = i >> 1;
    }
    return subset;
}

template <typename T>
T foldl(T (*func)(T,T), const T& t0, const vec_t<T>& ts)
{
    T result = t0;
    for (T t : ts)
        result = func(result, t);
    return result;
}


/* Functions for working with map_t ******************************************/

template <typename K, typename T>
inline
bool isIn(const map_t<K,T>& m, K key)
{
    return (m.count(key) > 0);
}

template <typename K, typename T>
T sumValues(const map_t<K,T>& m)
{
    T sum {0};
    for (auto it = m.begin(); it != m.end(); ++it)
        sum += it -> second;
    return sum;
}

template <typename K, typename T>
void printMap(const map_t<K,T>& m)
{
    std::cout << "{";
    for (auto pair : m)
        std::cout << "{" << pair.first << ": " << pair.second << "} ";
    std::cout << "}\n";
}

template <typename K, typename V>
vec_t<K> getKeys(const map_t<K,V>& m)
{
    vec_t<K> keys;
    for (auto p : m) {
        keys.push_back(p.first);
    }
    return keys;
}



/* Time tracking *************************************************************/

using wclock = std::chrono::steady_clock;
using timepoint_t = std::chrono::time_point<wclock>;
using duration_t = std::chrono::duration<double>;

double durationToDouble(duration_t d)
{
    return d.count() * duration_t::period::num / duration_t::period::den;
}

#if TRACK_TIME
constexpr int TIME_STACK_CAPACITY = 4096;
int TIME_STACK_PTR = -1;

#if TRACK_TIME_WALL
timepoint_t TIME_STACK[TIME_STACK_CAPACITY];  // wall time
#else
std::clock_t TIME_STACK[TIME_STACK_CAPACITY]; // processor time
#endif

str_t TIME_LOG_FILENAME = "time-log_" + timestamp();
std::ofstream TIME_LOG(TIME_LOG_FILENAME);

map_t<str_t, uint_t> CALL_COUNTS;
map_t<str_t, double> TIME_SUMS;
map_t<str_t, map_t<str_t,uint_t>> MSG_COUNTS;
map_t<str_t, map_t<str_t,double>> TIMES_BY_MSG;

void setTimeLogPrecision(int p=9)
{
    TIME_LOG << std::fixed << std::setprecision(p);
}

void setTimeLogFile(str_t filename, bool deleteOld=true)
{
    TIME_LOG = std::ofstream(filename);
    assert ( TIME_LOG );
    if (deleteOld)
        std::remove(TIME_LOG_FILENAME.c_str());
    TIME_LOG_FILENAME = filename;
}

void clearTimeTrackingData()
{
    TIME_STACK_PTR = -1;
    CALL_COUNTS.clear();
    TIME_SUMS.clear();
    MSG_COUNTS.clear();
    TIMES_BY_MSG.clear();
}

void pushTime()
{
    assert ( TIME_STACK_PTR + 1 < TIME_STACK_CAPACITY );
    TIME_STACK_PTR++;
#if TRACK_TIME_WALL
    TIME_STACK[TIME_STACK_PTR] = wclock::now();
#else
    TIME_STACK[TIME_STACK_PTR] = std::clock();
#endif
}

double popTime()
{
    assert ( TIME_STACK_PTR >= 0 );
#if TRACK_TIME_WALL
    duration_t deltaW = wclock::now() - TIME_STACK[TIME_STACK_PTR];
    TIME_STACK_PTR--;
    return durationToDouble(deltaW);
#else
    std::clock_t deltaP = std::clock() - TIME_STACK[TIME_STACK_PTR];
    TIME_STACK_PTR--;
    return (double)deltaP / CLOCKS_PER_SEC;
#endif
}

inline
void incrementOrInit(map_t<str_t, uint_t>& m, str_t key, uint_t init=1)
{
    m[key] = isIn(m, key) ? m[key] + 1 : init;
}

void logTimeStart(str_t name="", str_t msg="", bool writeToLog=false)
{
    if (writeToLog) {
        TIME_LOG << "start [" << TIME_STACK_PTR+1 << "] "
                 << name << " " << msg << "\n";
    }
    incrementOrInit(CALL_COUNTS, name);
    if (isIn(MSG_COUNTS, name))
        incrementOrInit(MSG_COUNTS[name], msg);
    else
        MSG_COUNTS[name] = {{msg, 1}};
    pushTime();
}

void logTimeEnd(str_t name="", str_t msg="", bool writeToLog=false)
{
    double delta = popTime();
    
    if (isIn(TIME_SUMS, name)) {
        TIME_SUMS[name] += delta;
    } else {
        TIME_SUMS[name] = delta;
    }
    
    if ( ! isIn(TIMES_BY_MSG[name], msg))
        TIMES_BY_MSG[name][msg] = 0.0;
    TIMES_BY_MSG[name][msg] += delta;
    
    if (writeToLog) {
        TIME_LOG << "end   [" << TIME_STACK_PTR+1 << "] " << delta
                 << " " << name << "\n";
    }
}

void logTimeTrackingSummaries(bool flush=false)
{
    TIME_LOG << "\nSummary:\n{\n";
    for (auto c = CALL_COUNTS.begin(); c != CALL_COUNTS.end(); ) {
        auto entry = *c;
        str_t name = entry.first;
        TIME_LOG << "  \"" << name << "\": {\n"
                 << "    \"calls\": " << entry.second << ",\n";
        if (isIn(TIME_SUMS, name)) {
            TIME_LOG << "    \"time\": "<< TIME_SUMS[name] <<",\n";
        }
        if (isIn(MSG_COUNTS, name)) {
            TIME_LOG << "    \"msg_counts\": {\n";
            for (auto m = MSG_COUNTS[name].begin();
                 m != MSG_COUNTS[name].end(); )
            {
                str_t msg = m -> first;
                uint_t count = m -> second;
                TIME_LOG << "      \"" << msg << "\" : " << count;
                if ( ++m != MSG_COUNTS[name].end() )
                    TIME_LOG << ",";
                TIME_LOG << "\n";
            }
            TIME_LOG << "    },\n";
        }
        if (isIn(TIMES_BY_MSG, name)) {
            TIME_LOG << "    \"times_by_msg\": {\n";
            for (auto m = TIMES_BY_MSG[name].begin();
                 m != TIMES_BY_MSG[name].end(); )
            {
                str_t msg = m -> first;
                double time = m -> second;
                TIME_LOG << "      \"" << msg << "\" : " << time;
                if ( ++m != TIMES_BY_MSG[name].end() )
                    TIME_LOG << ",";
                TIME_LOG << "\n";
            }
            TIME_LOG << "    }\n";
        }
        TIME_LOG << "  }";
        if ( ++c != CALL_COUNTS.end() )
            TIME_LOG << ",";
        TIME_LOG << "\n";
    }
    TIME_LOG << "\n}\n";
    
    if (flush)
        clearTimeTrackingData();
}
#endif  // TRACK_TIME


/* Memory management *********************************************************/

// TODO: Implement memory tracking data structures

template <typename T>
T* arrayAllocate(scalar_t n)
{
    // TODO: Implement memory tracking
    T* p = new T[n];
    return p;
}

template <typename T>
void arrayDelete(const T* p)
{
    // TODO: Implement memory tracking
    delete[] p;
}

template <typename T>
struct ArrayDeleter {
    void operator() (const T* p) {
        arrayDelete(p);
    }
};

template <typename T>
void arrayCopy(index_t n, const T* f, T* g)
{
    for(index_t i = 0; i < n; ++i)
        g[i] = f[i];
}


/* Miscellaneous functions ***************************************************/

/* Generate a float uniformly at random in range [min, max].
 * if min <= max, generate a float uniformly at random in range
 * [std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max()].
 */
float genRandFloatInRange(float min=0.0, float max=0.0)
{
    assert ( std::isfinite(max) ); assert ( std::isfinite(min) );
    float f;
    double maxFloat = std::numeric_limits<float>::max();
    double minFloat = std::numeric_limits<float>::lowest();
    if (max > min) {
        maxFloat = max;
        minFloat = min;
    }
    uint_t maxUint = std::numeric_limits<uint_t>::max();
    uint_t r = (uint_t)(rand() << 1);
    double relativePos = (double)r / ((double)maxUint);
    f = (float)(relativePos * (maxFloat - minFloat) + minFloat);
    assert ( std::isfinite(f) );
    return f;
}

/* Generate 32- or 64-bit sequences uniformly at random, cast to type T.
 * NOTE: Not necessarily the same as generating numbers 
 *       of type T uniformly at random! */
template <typename T>
T genFromRandBits()
{
    size_t nBytes = sizeof(T);
    assert (nBytes == 4 || nBytes == 8);
    
    long int r = rand();
    int sign = rand() % 2;
    if (sign)
        r *= -1;
    if (nBytes == 8) {
        r = r << 32;
        r += rand();
    }
    T s = *((T*)(&r));
    return s;
}


template <typename T>
T genRandScalar() {
    if constexpr (std::is_same_v<T, float>)
        return genRandFloatInRange(-1.0, 1.0);
    if constexpr (std::is_same_v<T, double>)
        return (double)genRandFloatInRange(-1.0, 1.0);
    else if constexpr (std::is_same_v<T, uint_t>)
        return (((uint_t)rand())<<1) + ((uint_t)rand()%2);
    else if constexpr (std::is_same_v<T, int>)
        return rand();
    else
        return T((uint_t)rand());
}


template <typename P>
Zp<P> genRandScalar(Zp<P> min, Zp<P> max)
{
    #if DEBUG_ASSERTS
    assert ( min <= max );
    #endif
    uint_t diff = (max - min).value();
    return min + Zp<P>((uint_t)rand() % (1 + diff));
}

float genRandScalar(float min, float max)
{ return genRandFloatInRange(min, max); }

uint_t genRandScalar(uint_t min, uint_t max)
{
    #if DEBUG_ASSERTS
    assert ( min <= max );
    #endif
    uint_t r = (((uint_t)rand())<<1) + ((uint_t)rand()%2);
    return min + (r % (1 + max - min));
}

int genRandScalar(int min, int max)
{
    #if DEBUG_ASSERTS
    assert ( min <= max );
    #endif
    return min + (rand() % (1 + max - min));
}


long pow(long base, uint_t exp)
{
    assert ( ! (base == 0 && exp == 0) );
    long result = 1;
    for (uint i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

long_uint_t factorial(long_uint_t x)
{
    if (x <= 1)
        return 1;
    assert ( x <= 20 );  // else the result does not fit into long_uint_t
    return x * factorial(x-1);
}



/* Functions for working with arrays *****************************************/

template <typename T>
T arrayProduct(size_t n, const T* elements)
{
    T p {1};
    for (size_t i = 0; i < n; ++i) {
        p *= elements[i];
    }
    return p;
}

template <typename T>
void printArray(size_t n, const T* arr)
{
    std::cout << "[";
    for (size_t i = 0; i < n; ++i) {
        if (i == n-1)
            std::cout << arr[i] << "]";
        else
            std::cout << arr[i] << ", ";
    }
    std::cout << "\n";
}

template <typename T>
void printNDArray(size_t nDims, const scalar_t* dimSizes, const T* arr,
                  size_t level=0, long_scalar_t maxVol=2000)
{
    assert ( arr != NULL );
    assert ( (nDims >= 0) && (level >= 0) );
    if (nDims == 0) {
        std::cout << "Scalar(" << arr[0] << ")\n";
        return;
    }
    long_scalar_t vol = arrayProduct(nDims - 1, dimSizes + 1);
    if (level == 0) {
        if (vol * dimSizes[0] > maxVol) {
            std::cout << "printNDArray: array volume exceeds maxVol ("
                << maxVol <<"); aborting.\n";
            return;
        }
    }
    
    for (size_t l = 0; l < level; ++l) std::cout << "  ";
    if (nDims == 1) {
        printArray(dimSizes[0], arr);
        return;
    }
    std::cout << "[\n";
    for (size_t a = 0; a < dimSizes[0]; ++a) {
        printNDArray(nDims - 1, dimSizes + 1, arr + a*vol, level + 1);
    }
    for (size_t l = 0; l < level; ++l) std::cout << "  ";
    std::cout << "]\n";
}

template <typename T>
void naiveMatrixMul(uint_t na, uint_t nc, uint_t nb,
                    const T* matA, const T* matB, T* matC)
{
    #if TRACK_TIME >= 4
    str_t msg = std::to_string(na) 
                + "," + std::to_string(nc)
                + "," + std::to_string(nb);
    logTimeStart("naiveMatrixMul", msg);
    #endif
    // matA : na x nc
    // matB : nb x nc (transposed)
    // matC : na x nb
    uint_t rnc {0}, cnc {0};
    for (uint_t r = 0; r < na; ++r) {
        rnc = r * nc;
        for (uint_t c = 0; c < nb; ++c) {
            cnc = c * nc;
            T acc = 0L;
            for (uint_t i = 0; i < nc; ++i) {
                acc += matA[rnc + i] * matB[cnc + i];
            }
            matC[r*nb + c] = acc;
        }
    }
    #if TRACK_TIME >= 4
    logTimeEnd("naiveMatrixMul", msg);
    #endif
}


/* Translate single linear coordinate into nDims-tuple of coordinates */
void getCoords(long_scalar_t k, scalar_t nDims, const scalar_t* dimSizes,
               scalar_t* coords)
{
    long_scalar_t rem = k;
    long_scalar_t vol = arrayProduct<scalar_t>(nDims, dimSizes);
    #if DEBUG_ASSERTS
    assert ( k < vol );
    #endif
    scalar_t q;
    for (scalar_t d = 0; d < nDims; ++d) {
        vol /= dimSizes[d];
        q = rem / vol;
        coords[d] = q;
        rem -= q * vol;
    }
}

long_scalar_t linearizeCoords(scalar_t nDims, const scalar_t* dimSizes,
                              const scalar_t* coords)
{
    #if TRACK_TIME >= 4
    str_t volStr = std::to_string(arrayProduct<scalar_t>(nDims,dimSizes));
    logTimeStart("linearizeCoords", volStr);
    #endif
    long_scalar_t linearCoord {0};
    long_scalar_t vol = arrayProduct<scalar_t>(nDims, dimSizes);
    for (scalar_t d = 0; d < nDims; ++d) {
        vol /= dimSizes[d];
        linearCoord += vol * coords[d];
    }
    #if TRACK_TIME >= 4
    logTimeEnd("linearizeCoords", volStr);
    #endif
    return linearCoord;
}

/* Construct an array that is a copy but with one element removed */
template <typename T>
void elideElementAt(scalar_t ind, scalar_t n, const T* arr, T* out)
{
    scalar_t j = 0;
    for (scalar_t i = 0; i < n; ++i) {
        if (i == ind)
            continue;
        else {
            out[j] = arr[i];
            ++j;
        }
    }
}

/* NOTE: Inefficient */
template <typename T>
void arraySliceAtDim(scalar_t sliceDim, scalar_t sliceInd,
                     scalar_t nDims, const scalar_t* dimSizes,
                     const T* arr, T* output)
{
    scalar_t* outDimSizes = new scalar_t[nDims-1];
    elideElementAt(sliceDim, nDims, dimSizes, outDimSizes);
    long_scalar_t outVol = arrayProduct(nDims-1, outDimSizes);
    assert (outVol == arrayProduct(nDims, dimSizes) / dimSizes[sliceDim]);
    scalar_t* outCoords = new scalar_t[nDims-1];
    scalar_t* arrCoords = new scalar_t[nDims];
    
    for (long_scalar_t outPos = 0; outPos < outVol; ++outPos) {
        getCoords(outPos, nDims-1, outDimSizes, outCoords);
        arrCoords[sliceDim] = sliceInd;
        scalar_t j = 0;
        for (scalar_t i = 0; i < nDims-1; ++i) {
            if (i == sliceDim)
                ++j;
            arrCoords[j] = outCoords[i];
            ++j;
        }
        long_scalar_t arrPos = linearizeCoords(nDims, dimSizes, arrCoords);
        output[outPos] = arr[arrPos];
    }
    delete[] outDimSizes;
    delete[] outCoords;
    delete[] arrCoords;
}

/* NOTE: Inefficient */
/* dimsToInds: map from indices of slice dimensions to indices within
 * those dimensions. E.g. {{1,3}, {7,99}} would slice dimension 1 along
 * position 3 and dimension 7 along positions 99.
 */
template <typename T>
void arraySlice(map_t<uint_t, uint_t> dimsToInds,
                vec_t<scalar_t> dimSizes, const T* arr,
                T* output)
{
    assert (dimsToInds.size() <= dimSizes.size());
    
    if (dimsToInds.size() == 0) {
        index_t vol = (index_t)vectorProduct(dimSizes);
        arrayCopy(vol, arr, output);
        return;
    }
    
    if (dimsToInds.size() == 1) {
        uint_t sliceDim = dimsToInds.begin() -> first;
        uint_t sliceInd = dimsToInds.begin() -> second;
        arraySliceAtDim(sliceDim, sliceInd,
                        dimSizes.size(), dimSizes.data(),
                        arr, output);
        return;
    }
    
    // Recursive call: Slice along all but first slice dimension
    
    vec_t<uint_t> sliceDims;
    for (auto slicePair : dimsToInds)
        sliceDims.push_back(slicePair.first);
    std::sort(sliceDims.begin(), sliceDims.end());
    
    uint_t curDim = sliceDims[0];
    
    long_scalar_t recOutputVol = vectorProduct(dimSizes);
    for (auto slicePair : dimsToInds)
        recOutputVol /= dimSizes[slicePair.first];
    recOutputVol *= dimSizes[curDim];
    T* recOutput = new T[recOutputVol];
    
    map_t<uint_t, uint_t> recDimsToInds(dimsToInds);
    recDimsToInds.erase(curDim);
    
    arraySlice(recDimsToInds, dimSizes, arr, recOutput);
    
    // Slice output of recursive call along remaining (first) slice dimension
    
    uint_t d {0};
    for (auto it = dimSizes.begin(); it != dimSizes.end(); ) {
        if (isIn(sliceDims, d) && d != curDim)
            it = dimSizes.erase(it);
        else
            ++it;
        ++d;
    }
    arraySliceAtDim<T>(curDim, dimsToInds.at(curDim),
                       dimSizes.size(), dimSizes.data(),
                       recOutput, output);
    
    delete[] recOutput;
}



/* Permutation for vec_t and arrays ******************************************/

void generateRandPerm(uint_t n, uint_t* perm)
{
    for (uint_t i = 0; i < n; ++i) {
        perm[i] = i;
    }
    for (uint_t i = 0; i < n; ++i) {
        uint_t swapInd = i + ((uint_t)rand() % (n-i));
        uint_t swapVal = perm[swapInd];
        perm[swapInd] = perm[i];
        perm[i] = swapVal;
    }
}
    
void inversePermutation(uint_t n, const uint_t* perm, uint_t* inverse)
{
    for (uint_t i = 0; i < n; ++i) {
        inverse[perm[i]] = i;
    }
}


template <typename T>
void permuteArray(uint_t nDims, const uint_t* perm,
                  const uint_t* dimSizesA, const T* A, T* B)
{
    uint_t* dimSizesB = new uint_t[nDims];
    for (uint_t d = 0; d < nDims; ++d) {
        dimSizesB[d] = dimSizesA[perm[d]];
    }
    
    long_uint_t K = arrayProduct(nDims, dimSizesA);
    uint_t* coordsA = new uint_t[nDims];
    uint_t* coordsB = new uint_t[nDims];
    for (long_uint_t k = 0; k < K; ++k) {
        getCoords(k, nDims, dimSizesA, coordsA);
        for (uint_t d = 0; d < nDims; ++d) {
            coordsB[d] = coordsA[perm[d]];
        }
        long_uint_t kB = linearizeCoords(nDims, dimSizesB, coordsB);
        B[kB] = A[k];
    }
    delete[] coordsA;
    delete[] coordsB;
    delete[] dimSizesB;
}

template <typename T>
vec_t<T> permuteVector(uint_t nDims, const uint_t* perm,
                             const vec_t<T>& vec)
{
    vec_t<T> permutedVec(vec);
    for (uint_t d = 0; d < nDims; ++d) {
        permutedVec[d] = vec[perm[d]];
    }
    return permutedVec;
}

template <typename T>
void permuteArray(uint_t nDims, const vec_t<uint_t>& perm,
                  const vec_t<uint_t>& dimSizesA, const T* A, T* B)
{
    permuteArray<T>(nDims, perm.data(), dimSizesA.data(), A, B);
}



/* Miscellaneous linear algebra for working with Vandermonde matrices ********/


/* Compute LU decomposition of a Vandermonde matrix.
 * NOTE: This implementation is inefficient, and assumes that
 * the input is a square n-by-n Vandermonde matrix.
 */
template <typename T>
void luDecomposeVandermonde(scalar_t n, const T* vm, T* lower, T* upper)
{
    // check that first element of each row is 1
    for (scalar_t r = 0; r < n; ++r) {
        assert (vm[r*n] == T::one);
    }
    
    arrayZero(n*n, lower);
    arrayCopy(n*n, vm, upper);
    
    for (scalar_t i = 0; i < n; ++i) {
        // TODO: (dis)prove: never need to pivot because Vandermonde?
        assert (upper[i*n + i] != T::zero);
        lower[i*n + i] = T::one;
        for (scalar_t r = i+1; r < n; ++r) {
            T ratio = upper[r*n + i] / upper[i*n + i];
            lower[r*n + i] = ratio;
            for (scalar_t c = i; c < n; ++c)
                upper[r*n + c] = upper[r*n + c] - ratio * upper[i*n + c];
        }
    }
}

/* tri: n-by-n lower triangular matrix with nonzero diagonal
 * y: array of n elements
 * solve for x: tri * x = y
 */
template <typename T>
void solveLowerTriangular(scalar_t n, const T* tri, const T* y, T* x)
{
    for (scalar_t i = 0; i < n; ++i) {
        assert (tri[i*n + i] != T::zero);
        T sum = T::zero;
        for (scalar_t c = 0; c < i; ++c) {
            sum += x[c] * tri[i*n + c];
        }
        x[i] = (y[i] - sum) / tri[i*n + i];
    }
}

/* tri: n-by-n upper triangular matrix with nonzero diagonal
 * y: array of n elements
 * solve for x: tri * x = y
 */
template <typename T>
void solveUpperTriangular(scalar_t n, const T* tri, const T* y, T* x)
{
    for (index_t i = n-1; i >= 0; --i) {
        assert (tri[i*n + i] != T::zero);
        T sum = T::zero;
        for (index_t c = n-1; c > i; --c) {
            sum += x[c] * tri[i*n + c];
        }
        x[i] = (y[i] - sum) / tri[i*n + i];
    }
}

// NOTE: Inefficient
template <typename T>
void invertVandermonde(scalar_t n, const T* vm, T* vminv, bool transpose=false)
{
    T* lower = new T[n*n];
    T* upper = new T[n*n];
    luDecomposeVandermonde(n, vm, lower, upper);
    
    T* uVminvT = new T[n*n];
    T* y = new T[n];
    for (scalar_t c = 0; c < n; ++c) {
        arrayZero(n, y);
        y[c] = T::one;
        solveLowerTriangular<T>(n, lower, y, uVminvT + c*n);
    }
    
    T* vminvT = new T[n*n];
    for (scalar_t c = 0; c < n; ++c) {
        solveUpperTriangular<T>(n, upper, uVminvT + n*c, vminvT + n*c);
    }
    
    if (!transpose) {
        scalar_t perm[2] {1,0};
        scalar_t dimSizes[2] {n,n};
        permuteArray<T>(2, perm, dimSizes, vminvT, vminv);
    }
    else {
        arrayCopy(n*n, vminvT, vminv);
    }
    
    delete[] lower;
    delete[] upper;
    delete[] uVminvT;
    delete[] y;
    delete[] vminvT;
}

/* Compute n-by-n Vandermonde matrix from n points */
template <typename T>
void genVandermonde(scalar_t n, const T* points, T* vm)
{
    for (scalar_t r = 0; r < n; ++r) {
        T pPowC = T::one;
        T p = points[r];
        for (scalar_t c = 0; c < n; ++c) {
            vm[n*r + c] = pPowC;
            pPowC *= p;
        }
    }
}

/* Compute inverse of n-by-n Vandermonde matrix from n points */
template <typename T>
void genInvVandermonde(scalar_t n, const T* points, T* invvm)
{
    T* vm = new T[n*n];
    genVandermonde(n, points, vm);
    invertVandermonde(n, vm, invvm);
    delete[] vm;
}

template <typename T>
void genRandVandermonde(scalar_t n, T* vm)
{
    vec_t<T> points(n);
    for (scalar_t i = 0; i < n; ++i) {
        T p = T::rand();
        while (isIn(points, p))
            p = T::rand();
        points[i] = p;
    }
    genVandermonde(n, points.data(), vm);
}



#endif  // __IDIGM_UTIL_H__
