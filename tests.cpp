#include <iostream>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <type_traits>
#include "montgomery.cpp"
#include "util.cpp"
#include "factor.cpp"
#include "factorgraph.cpp"
#include "poly.cpp"
#include "eval_algo.cpp"
#include "crt.cpp"
#include "fileio.cpp"
#include "integration.cpp"
#include "primes_list.h"
#include "testbench.cpp"


using MP = MP2147483629;


template <typename MP>
scalar_t get_mod_inv_R()
{
    // R = 1L<<32
    scalar_t mod_inv_R = mod_reduce<MP>(
        first_bezout_factor(1L<<32, MP::modulus)
    );
    // R * mod_inv_R === 1 (mod MP::modulus)
    assert ( ((long_scalar_t)mod_inv_R << 32) % MP::modulus == 1 );
    
    return mod_inv_R;
}


template <typename MP>
void test_Montgomery_arithmetic(bool verbose=false, size_t nSamples=200)
{
    scalar_t mod_inv_R = get_mod_inv_R<MP>();
    
    // mod_reduce
    if (verbose)
        std::cout << "Testing mod_reduce<" << typeid(MP).name() << "> ...\n";
    {
        long_signed_t r {0};
        long_signed_t m {0};
        for (size_t t = 0; t < nSamples; ++t) {
            r = (long_signed_t) rand();
            m = mod_reduce<MP>(r);
            assert ( m >= 0 );
            assert ( m < MP::modulus );
            if (r < 0)
                assert ( (m - r) % MP::modulus == 0);
            else
                assert ( (r - m) % MP::modulus == 0);
        }
    }
    
    // mod_rand
    if (verbose)
        std::cout << "Testing mod_rand<" << typeid(MP).name() << "> ...\n";
    {
        double acc {0};
        size_t n {1000};
        long_signed_t r {0};
        for (size_t i = 0; i < n; ++i) {
            r = mod_rand<MP>();
            assert ( r >= 0 );
            assert ( r < MP::modulus );
            acc += (double) r;
        }
        acc /= n;
        assert ( acc > 0.3 * (double)MP::modulus );
        assert ( acc < 0.7 * (double)MP::modulus );
    }
    
    // mod_neg
    if (verbose)
        std::cout << "Testing mod_neg<" << typeid(MP).name() << "> ...\n";
    {
        for (size_t t = 0; t < nSamples; ++t) {
            scalar_t r = mod_rand<MP>();
            scalar_t nr = mod_neg<MP>(r);
            assert ( ((long_scalar_t)r + nr) % MP::modulus == 0 );
        }
    }
    
    // mod_inv
    if (verbose)
        std::cout << "Testing mod_inv<"<< typeid(MP).name() <<"> ...\n";
    {
        // multiplicative inverses guaranteed to exist iff MP::modulus is prime
        // TODO: proper check for primality
        assert ( MP::modulus % 2 != 0 );
        assert ( MP::modulus % 3 != 0 );
        assert ( MP::modulus % 5 != 0 );
        assert ( MP::modulus % 7 != 0 );
        for (size_t t = 0; t < nSamples; ++t) {
            scalar_t r = mod_rand<MP>();
            // zero does not have a multiplicative inverse
            if (r == 0) continue;
            scalar_t ri = mod_inv<MP>(r);
            assert ( ((long_scalar_t)r * ri) % MP::modulus == 1 );
        }
    }
    
    // montgomery_reduce
    if (verbose) {
        std::cout << "Testing montgomery_reduce<"
        << typeid(MP).name() <<"> ...\n";
    }
    {
        for (size_t t = 0; t < nSamples; ++t) {
            long_mont_t x = ((long_mont_t) mod_rand<MP>()) * mod_rand<MP>();
            mont_t m { montgomery_reduce<MP>(x) };
            assert ( m == ((x % MP::modulus) * mod_inv_R) % MP::modulus );
        }
    }
    
    // to_montgomery and from_montgomery
    if (verbose) {
        std::cout << "Testing to_montgomery and from_montgomery<"
        << typeid(MP).name() <<"> ...\n";
    }
    {
        for (size_t t = 0; t < nSamples; ++t) {
            scalar_t x = mod_rand<MP>();
            mont_t m = to_montgomery<MP>(x);
            assert ( m == ((long_scalar_t)x * (1UL<<32)) % MP::modulus );
            scalar_t y = from_montgomery<MP>(m);
            assert ( y == x );
        }
    }
}

template <typename MP>
void test_Zp_interface(bool verbose=false, size_t nSamples=200)
{
    scalar_t mod_inv_R = get_mod_inv_R<MP>();
    
    // Zp<MP>::rand
    if (verbose)
        std::cout << "Testing Zp<" << typeid(MP).name() << ">::rand ...\n";
    
    double acc {0};
    size_t n {1000};
    long_signed_t r {0};
    for (size_t i = 0; i < n; ++i) {
        r = Zp<MP>::rand().raw();
        assert ( r >= 0 );
        assert ( r < MP::modulus );
        acc += (double) r;
    }
    acc /= n;
    assert ( acc > 0.3 * (double)MP::modulus );
    assert ( acc < 0.7 * (double)MP::modulus );
    
    
    // Zp<MP> constructors and getters
    if (verbose) {
        std::cout << "Testing Zp<" << typeid(MP).name() << ">"
            << " constructors and getters ...\n";
    }
    
    for (size_t s = 0; s < nSamples; ++s) {
        
        scalar_t r {mod_rand<MP>()};
        
        Zp<MP> z1(r);
        assert ( z1.value() == r );
        assert ( z1.raw() == (((long_scalar_t)r) << 32) % MP::modulus );
        
        Zp<MP> z2(z1);
        assert ( z2.value() == r );
        assert ( z2.raw() == (((long_scalar_t)r) << 32) % MP::modulus );
        
        Zp<MP> z3(r, true);
        assert ( z3.raw() == r );
    }
    
    
    // Zp<MP> inversion and negation
    if (verbose) {
        std::cout << "Testing Zp<" << typeid(MP).name() << ">"
            << " inversion and negation ...\n";
    }
    
    for (size_t s = 0; s < nSamples; ++s) {
        
        Zp<MP> z = Zp<MP>::rand();
        
        // zero does not have a multiplicative inverse
        if (!(z.raw() == 0)) {
            Zp<MP> zi = z.inv();
            long_mont_t zi_r {zi.raw()};
            long_mont_t z_r {z.raw()};
            assert (
                (((((zi_r * z_r) % MP::modulus)
                * mod_inv_R) % MP::modulus) * mod_inv_R) % MP::modulus
                == 1
            );
            assert ( (z.inv() * z).value() == 1 );
        }
        
        Zp<MP> zn = -z;
        assert (
            (((zn.raw() + z.raw()) % MP::modulus) *  mod_inv_R) % MP::modulus
            == 0
        );
        assert ( (zn + z).value() == 0 );
    }

    
    // Zp<MP> assignment operators 
    if (verbose) {
        std::cout << "Testing Zp<" << typeid(MP).name() << ">"
            << " assignment operators ...\n";
    }
    
    for (size_t s = 0; s < nSamples; ++s) {
        
        Zp<MP> z1 = Zp<MP>::rand();
        scalar_t z1r = z1.raw();
        scalar_t z1v = z1.value();
        Zp<MP> x = z1;
        assert ( x.raw() == z1r );
        assert ( x.value() == z1v );
        
        Zp<MP> z2 = Zp<MP>::rand();
        scalar_t z2r = z2.raw();
        
        x += z2;
        assert ( x.raw() == ((long_mont_t)z1r + z2r) % MP::modulus );
        // check that (+=) did not affect z1 or z2
        assert ( z1.raw() == z1r );
        assert ( z1.value() == z1v );
        assert ( z2.raw() == z2r );
        
        x = z1;
        x -= z2;
        mont_t z2r_neg = (z2r == 0) ? 0 : MP::modulus - z2r;
        assert (
            ((long_mont_t)z1r + z2r_neg) % MP::modulus
            == x.raw()
        );
        assert ( (x + z2).raw() == z1r );
        // check that (-=) did not affect z1 or z2
        assert ( z1.raw() == z1r );
        assert ( z1.value() == z1v );
        assert ( z2.raw() == z2r );
        
        x = z1;
        x *= z2;
        assert (
            ((((long_mont_t)z1r*z2r) % MP::modulus) * mod_inv_R) % MP::modulus
            == x.raw()
        );
        assert ( z1.raw() == z1r );
        assert ( z1.value() == z1v );
        assert ( z2.raw() == z2r );
        
        // zero does not have a multiplicative inverse
        if (z2.raw() == 0) continue;
        x = z1;
        x /= z2;
        assert ( (x * z2).raw() == z1r );
        assert ( z1.raw() == z1r );
        assert ( z1.value() == z1v );
        assert ( z2.raw() == z2r );
    }
    
    
    // Zp<MP> arithmetic operators
    if (verbose) {
        std::cout << "Testing Zp<" << typeid(MP).name() << ">"
            << " arithmetic operators ...\n";
    }
    
    for (size_t s = 0; s < nSamples; ++s) {
        
        scalar_t r1 = mod_reduce<MP>(rand());
        scalar_t r2 = mod_reduce<MP>(rand());
        Zp<MP> z1(r1), z2(r2);
        
        assert (
            (z1 + z2).value()
            == ((long_scalar_t)r1 + r2) % MP::modulus
        );
        
        assert (
            (z1 - z2).value()
            == ((long_scalar_t)r1 + mod_neg<MP>(r2)) % MP::modulus
        );
        
        assert (
            (z1 * z2).value()
            == ((long_scalar_t)r1 * r2) % MP::modulus
        );
        
        // zero does not have a multiplicative inverse
        if (r2 == 0) continue;
        assert (
            (z1 / z2).value()
            == ((long_scalar_t)r1 * mod_inv<MP>(r2)) % MP::modulus
        );
    }
    
    
    // Zp<MP> comparison operators
    if (verbose) {
        std::cout << "Testing Zp<" << typeid(MP).name() << ">"
            << " comparison operators ...\n";
    }
    
    for (size_t s = 0; s < nSamples; ++s) {
        
        scalar_t r1 = mod_reduce<MP>(rand());
        scalar_t r2 = mod_reduce<MP>(rand());
        Zp<MP> z1(r1), z2(r2);
        
        assert ( z1 == z1 );
        Zp<MP> z3(z1);
        assert ( z1 == z3 ); assert ( z3 == z1 );
        Zp<MP> z4(r1);
        assert ( z1 == z4 ); assert ( z4 == z1 );
        
        assert ( (r1 != r2) == (z1 != z2) );
    }
    
    
    // Zp<MP> additive and multiplicative identity elements
    if (verbose) {
        std::cout << "Testing Zp<" << typeid(MP).name() << ">"
            << " additive and multiplicative identity elements ...\n";
    }
    
    for (size_t s = 0; s < nSamples; ++s) {
        
        Zp<MP> z = Zp<MP>::rand();
        assert ( z * Zp<MP>::one == z );
        assert ( z * Zp<MP>::zero == Zp<MP>::zero );
        assert ( z + Zp<MP>::zero == z );
    }
    
    
    // Zp<MP>::inv_2_pow_k
    if (verbose) {
        std::cout << "Testing Zp<" << typeid(MP).name()
            << ">::inv_2_pow_k...\n";
    }
    
    for (size_t s = 0; s < nSamples; ++s) {
        long_signed_t k = rand() % 11;
        Zp<MP> inv2pk = Zp<MP>::inv_2_pow_k(k);
        Zp<MP> m2pk = Zp<MP>::one;
        for (scalar_t i = 0; i < k; ++i)
            m2pk = m2pk * Zp<MP>(2U);
        assert ( m2pk * inv2pk == Zp<MP>::one );
    }
}


void test_arrayProduct() {
    size_t n = 5;
    int array[] {1,2,3,-4,5};
    assert ( arrayProduct<int>(n, array) == -120 );
}

void test_getCoords() {
    // TODO? More thorough tests
    scalar_t varSizes[] {3,2,4};
    constexpr size_t nVars {3};
    scalar_t coords[nVars];
    scalar_t expectedCoords[] {2, 0, 1};
    long_scalar_t k {17};
    
    getCoords(k, nVars, varSizes, coords);
    for (size_t d = 0; d < nVars; ++d)
        assert ( coords[d] == expectedCoords[d] );
    
    k = 13;
    scalar_t expectedCoords2[] = {1, 1, 1};
    getCoords(k, nVars, varSizes, coords);
    for (size_t d = 0; d < nVars; ++d)
        assert ( coords[d] == expectedCoords2[d] );
}

void test_arraySliceAtDim_fixed_example()
{
    constexpr scalar_t d0 {2}, d1 {3}, d2 {4};
    scalar_t dimSizes[] {d0, d1, d2};
    int arr[d0 * d1 * d2] {
        111,112,113,114,
        121,122,123,124,
        131,132,133,134,
        
        211,212,213,214,
        221,222,223,224,
        231,232,233,234
    };
    int output[d0 * d2];
    scalar_t sliceDim = 1;
    scalar_t sliceInd = 2;
    arraySliceAtDim(sliceDim, sliceInd, 3, dimSizes, arr, output);
    
    for (scalar_t i0 = 0; i0 < d0; ++i0) {
        for (scalar_t i2 = 0; i2 < d2; ++i2) {
            assert ( output[i0*4 + i2] == arr[i0*d1*d2 + sliceInd*d2 + i2] );
        }
    }
}

void test_arraySlice_fixed_examples()
{
    // Example 1
    
    vec_t<scalar_t> dimSizes {2,3,4,5,6};
    map_t<scalar_t, scalar_t> dimsToInds {{1,1}, {3,2}};
    scalar_t arrVol = vectorProduct(dimSizes);
    int* arr = new int[arrVol];
    for (scalar_t i = 0; i < arrVol; ++i)
        arr[i] = (int)i;
    constexpr scalar_t outVol = 2*4*6;
    int* output = new int[outVol];
    
    arraySlice<int>(dimsToInds, dimSizes, arr, output);
    
    int expectedOutput[outVol] = {
        132, 133, 134, 135, 136, 137,
        162, 163, 164, 165, 166, 167,
        192, 193, 194, 195, 196, 197,
        222, 223, 224, 225, 226, 227,

        492, 493, 494, 495, 496, 497,
        522, 523, 524, 525, 526, 527,
        552, 553, 554, 555, 556, 557,
        582, 583, 584, 585, 586, 587
    };
    for (scalar_t i = 0; i < outVol; ++i)
        assert ( output[i] == expectedOutput[i] );
    
    // Example 2
    
    map_t<scalar_t, scalar_t> dimsToInds2 = {{0,0}, {2,1}, {4,5}};
    constexpr scalar_t outVol2 = 3*5;
    int* output2 = new int[outVol2];
    
    arraySlice<int>(dimsToInds2, dimSizes, arr, output2);
    
    int expectedOutput2[outVol2] {
        35,  41,  47,  53,  59,
        155, 161, 167, 173, 179,
        275, 281, 287, 293, 299
    };
    for (scalar_t i = 0; i < outVol2; ++i)
        assert ( output2[i] == expectedOutput2[i] );
    
    // Example 3 : edge case when idsToInds is empty
    
    map_t<scalar_t, scalar_t> dimsToInds3 = {};
    constexpr scalar_t outVol3 = 2*3*4*5*6;
    assert ( outVol3 == arrVol );
    int* output3 = new int[outVol3];
    
    arraySlice<int>(dimsToInds3, dimSizes, arr, output3);

    for (scalar_t i = 0; i < outVol3; ++i)
        assert ( output3[i] == arr[i] );
    
    
    delete[] arr;
    delete[] output;
    delete[] output2;
    delete[] output3;
}

void test_permuteArray_fixed_examples()
{
    size_t nVars1 = 2;
    scalar_t perm1[] {1,0};  // transpose
    scalar_t varSizes1[] {3,2};
    int arrA1[6] {11,12, 21,22, 31,32};
    int arrB1[6] {0,0, 0,0, 0,0};
    
    permuteArray(nVars1, perm1, varSizes1, arrA1, arrB1);
    for (size_t r = 0; r < 3; ++r) {
        for (size_t c = 0; c < 2; ++c) {
            assert ( arrA1[r*2 + c] == arrB1[c*3 + r] );
        }
    }
    
    size_t nVars2 = 1;
    scalar_t perm2[] {0};  // identity
    scalar_t varSizes2[] {5,};
    int arrA2[5] {5,2,7,4,9};
    int arrB2[5] {0,0,0,0,0};
    
    permuteArray(nVars2, perm2, varSizes2, arrA2, arrB2);
    for (size_t c = 0; c < 5; ++c) {
        assert ( arrA2[c] == arrB2[c] );
    }
}


void test_perm_generation_and_inversion(bool verbose=false, size_t nSamples=200)
{
    if (verbose)
        std::cout << "Testing permutation generation and inversion...\n";
    
    size_t maxPermLength {20};
    size_t permLength {1};
    for (size_t s = 0; s < nSamples; ++s) {
        permLength = ((size_t)rand() % maxPermLength) + 1;
        scalar_t* perm = new scalar_t[permLength];
        generateRandPerm(permLength, perm);
        
        // Check that perm appears to be a valid permutation
        if (verbose) {
            std::cout << "perm: ";
            for (size_t i = 0; i < permLength; ++i)
                std::cout << perm[i] << " ";
            std::cout << "\n";
        }
        scalar_t sum {0};
        scalar_t expectedSum {0};
        for (size_t i = 0; i < permLength; ++i) {
            assert ( perm[i] < permLength );
            // assert ( perm[i] >= 0 );
            sum += perm[i];
            expectedSum += i;
        }
        assert ( sum == expectedSum );
        
        // test inversion
        scalar_t* inverse_perm = new scalar_t[permLength];
        inversePermutation(permLength, perm, inverse_perm);
        int* rand_arr = new int[permLength];
        int* permd_rand_arr = new int[permLength];
        int* invd_permd_rand_arr = new int[permLength];
        for (size_t i = 0; i < permLength; ++i)
            rand_arr[i] = rand();
        for (size_t i = 0; i < permLength; ++i)
            permd_rand_arr[i] = rand_arr[perm[i]];
        for (size_t i = 0; i < permLength; ++i)
            invd_permd_rand_arr[i] = permd_rand_arr[inverse_perm[i]];
        for (size_t i = 0; i < permLength; ++i)
            assert ( invd_permd_rand_arr[i] == rand_arr[i] );
        
        delete[] perm;
        delete[] inverse_perm;
        delete[] rand_arr;
        delete[] permd_rand_arr;
        delete[] invd_permd_rand_arr;
    }
}

template <typename F>
void test_permuteArray_perm_and_inverse(bool verbose=false, size_t nSamples=20)
{
    if (verbose) {
        std::cout << "Testing permuteArray: checking that "
            <<"inverse(P)(P(array)) == array, for random permutations P...\n";
    }
    size_t maxNVars {7};
    size_t nVars {1};
    size_t maxVarSize {8};
    
    for (size_t s = 0; s < nSamples; ++s) {
        
        nVars = ((size_t)rand() % maxNVars) + 1;
        scalar_t* perm = new scalar_t[nVars];
        generateRandPerm(nVars, perm);
        scalar_t* inverse_perm = new scalar_t[nVars];
        inversePermutation(nVars, perm, inverse_perm);
        
        scalar_t* varSizesA = new scalar_t[nVars];
        for (size_t d = 0; d < nVars; ++d)
            varSizesA[d] = ((scalar_t)rand() % maxVarSize) + 1;
        
        scalar_t* varSizesB = new scalar_t[nVars];
        for (size_t d = 0; d < nVars; ++d)
            varSizesB[d] = varSizesA[perm[d]];
        
        long_scalar_t volume = arrayProduct<scalar_t>(nVars, varSizesA);
        F* arrA = new F[volume];
        for (long_scalar_t k = 0; k < volume; ++k) {
            arrA[k] = F::rand();
        }
        F* arrB = new F[volume];
        F* arrC = new F[volume];
        
        permuteArray(nVars, perm, varSizesA, arrA, arrB);
        permuteArray(nVars, inverse_perm, varSizesB, arrB, arrC);
        
        for (long_scalar_t k = 0; k < volume; ++k) {
            assert ( arrC[k] == arrA[k] );
        }
        
        delete[] perm;
        delete[] inverse_perm;
        delete[] varSizesA;
        delete[] varSizesB;
        delete[] arrA;
        delete[] arrB;
        delete[] arrC;
    }
}


void test_permuteVector(bool verbose=false, size_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing permuteVector ...\n";

    size_t maxLength {20};
    size_t length {1};
    
    for (size_t s = 0; s < nSamples; ++s) {
        length = ((size_t)rand() % maxLength) + 1;
        vec_t<int> vec (length, 0);
        for (size_t d = 0; d < length; ++d) {
            vec[d] = rand() % 10;
        }
        scalar_t* perm = new scalar_t[length];
        generateRandPerm(length, perm);
        vec_t<int> permdVec = permuteVector(length, perm, vec);
        if (verbose) {
            std::cout << "original vector: "; printVector(vec);
            std::cout << "permutation: "; printArray(length, perm);
            std::cout << "permuted vector: "; printVector(permdVec);
        }
        
        assert ( permdVec.size() == vec.size() );
        for (size_t i = 0; i < length; ++i) {
            assert ( permdVec[i] == vec[perm[i]] );
        }
        
        delete[] perm;
    }
}


void test_printNDArray_fixed_examples()
{
    std::cout << "Testing printNDArray by printing fixed examples...\n";
    size_t nVars = 3;
    scalar_t* varSizes = new scalar_t[3] {2,3,4};
    int* arr = new int[24] {
        111,112,113,114,
        121,122,123,124,
        131,132,133,134,
        
        211,212,213,214,
        221,222,223,224,
        231,232,233,234
    };
    printNDArray<int>(nVars, varSizes, arr);
    delete[] varSizes;
    delete[] arr;
}

void test_indexOf_fixed_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing indexOf with fixed examples...\n";
    
    vec_t<int> v {13,-66,91,6,444,91};
    assert ( indexOf(91, v) == 2 );
    assert ( indexOf(13, v) == 0 );
    assert ( indexOf(66, v) == -1 );
    assert ( indexOf(99999, v) == -1 );
}


void test_isSuperset_fixed_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing isSuperset with fixed examples...\n";
    
    assert ( isSuperset(
                vec_t<int>{1,2,3,4}, 
                vec_t<int>{1,2,3} ));
    
    assert ( isSuperset(
                vec_t<int>{4,3,2,1}, 
                vec_t<int>{1,2,3} ));
    
    assert ( isSuperset(
                vec_t<int>{4,3,2,1}, 
                vec_t<int>{4,1,2,3} ));
    
    assert ( isSuperset(
                vec_t<int>{4,3,2,1}, 
                vec_t<int>{3,4,1,3,3,3,3} ));
    
    assert ( ! isSuperset(
                vec_t<int>{4,3,2,1}, 
                vec_t<int>{5,1,2,3} ));
    
    assert ( ! isSuperset(
                vec_t<int>{4,3,2,1}, 
                vec_t<int>{5,1,2,4,3} ));
}


void test_moveToEnd_fixed_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing moveToEnd with fixed examples...\n";
    vec_t<int> host {1,2,3,4,5,6};
    vec_t<int> tail {4,1,2};
    vec_t<int> m = moveToEnd(host, tail);
    if (verbose) {
        std::cout << "moveToEnd(\n    "; printVector(host);
        std::cout << "  , "; printVector(tail); std::cout << ") :\n    ";
        printVector(m);
    }
    assert ( m[0] == 3 );
    assert ( m[1] == 5 );
    assert ( m[2] == 6 );
    assert ( m[3] == 4 );
    assert ( m[4] == 1 );
    assert ( m[5] == 2 );
    
    //NOTE: moveToEnd deletes duplicates
    vec_t<int> host2 {1,4,4,2,4,3,4,5,6};
    vec_t<int> tail2 {4,1,2};
    vec_t<int> m2 = moveToEnd(host2, tail2);
    if (verbose) {
        std::cout << "moveToEnd(\n    "; printVector(host2);
        std::cout << "  , "; printVector(tail2); std::cout << ") :\n    ";
        printVector(m);
    }
    assert ( m2[0] == 3 );
    assert ( m2[1] == 5 );
    assert ( m2[2] == 6 );
    assert ( m2[3] == 4 );
    assert ( m2[4] == 1 );
    assert ( m2[5] == 2 );
}


void test_getPermutation_fixed_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing getPermutation with fixed examples...\n";
    
    vec_t<int> v1 {1,2,3,4,5};
    vec_t<int> w1 {5,2,3,4,1};
    vec_t<scalar_t> p = getPermutation(v1, w1);
    if (verbose) {
        std::cout << "getPermutation(\n    "; printVector(v1);
        std::cout << "  , "; printVector(w1); std::cout << ") :\n    ";
        printVector(p);
    }
    assert ( p[0] == 4 );
    assert ( p[1] == 1 );
    assert ( p[2] == 2 );
    assert ( p[3] == 3 );
    assert ( p[4] == 0 );
    
    vec_t<int> v2 {1,2,3,4,5};
    vec_t<int> w2 {3,1,5,2,4};
    vec_t<scalar_t> p2 = getPermutation(v2, w2);
    if (verbose) {
        std::cout << "getPermutation(\n    "; printVector(v2);
        std::cout << "  , "; printVector(w2); std::cout << ") :\n    ";
        printVector(p2);
    }
    assert ( p2[0] == 2 );
    assert ( p2[1] == 0 );
    assert ( p2[2] == 4 );
    assert ( p2[3] == 1 );
    assert ( p2[4] == 3 );
    
}


template <typename T>
void test_naiveMatrixMul_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing naiveMatrixMul<" << typeid(T).name()
        << "> with a few examples...\n";
    
    scalar_t na {3};
    scalar_t nb {4};
    scalar_t nc {2};
    T* matA = new T[na*nc] {
        T(1), T(0),
        T(0), T(1),
        T(1), T(1)
    };
    T* matB = new T[nb*nc] {
        T(1), T(2),
        T(3), T(4),
        T(5), T(6),
        T(7), T(8)
    };
    T* matC = new T[na*nb];
    naiveMatrixMul(na, nc, nb, matA, matB, matC);
    
    if (verbose) {
        scalar_t* varSizes = new scalar_t[2] {na, nc};
        std::cout << "matA: "; printNDArray(2, varSizes, matA);
        varSizes[0] = nb;
        std::cout << "matB: "; printNDArray(2, varSizes, matB);
        varSizes[0] = na; varSizes[1] = nb;
        std::cout << "matC: "; printNDArray(2, varSizes, matC);
        delete[] varSizes;
    }
    assert ( matC[0] == matA[0]*matB[0] + matA[1]*matB[1] );
    assert ( matC[1] == matA[0]*matB[2] + matA[1]*matB[3] );
    assert ( matC[2] == matA[0]*matB[4] + matA[1]*matB[5] );
    assert ( matC[7] == matA[2]*matB[6] + matA[3]*matB[7] );
    assert ( matC[8] == matA[4]*matB[0] + matA[5]*matB[1] );
    assert ( matC[10] == matA[4]*matB[4] + matA[5]*matB[5] );
    
    delete[] matA;
    delete[] matB;
    delete[] matC;
    
    scalar_t max {100};
    
    T* matAr = new T[na*nc] {
        T((scalar_t)rand()%max), T((scalar_t)rand()%max),
        T((scalar_t)rand()%max), T((scalar_t)rand()%max),
        T((scalar_t)rand()%max), T((scalar_t)rand()%max)
    };
    T* matBr = new T[nb*nc] {
        T((scalar_t)rand()%max), T((scalar_t)rand()%max),
        T((scalar_t)rand()%max), T((scalar_t)rand()%max),
        T((scalar_t)rand()%max), T((scalar_t)rand()%max),
        T((scalar_t)rand()%max), T((scalar_t)rand()%max)
    };
    T* matCr = new T[na*nb];
    naiveMatrixMul(na, nc, nb, matAr, matBr, matCr);
    
    if (verbose) {
        scalar_t* varSizes = new scalar_t[2] {na, nc};
        std::cout << "matAr: "; printNDArray(2, varSizes, matAr);
        varSizes[0] = nb;
        std::cout << "matBr: "; printNDArray(2, varSizes, matBr);
        varSizes[0] = na; varSizes[1] = nb;
        std::cout << "matCr: "; printNDArray(2, varSizes, matCr);
        delete[] varSizes;
    }
    
    assert ( matCr[0] == matAr[0]*matBr[0] + matAr[1]*matBr[1] );
    assert ( matCr[1] == matAr[0]*matBr[2] + matAr[1]*matBr[3] );
    assert ( matCr[2] == matAr[0]*matBr[4] + matAr[1]*matBr[5] );
    assert ( matCr[7] == matAr[2]*matBr[6] + matAr[3]*matBr[7] );
    assert ( matCr[8] == matAr[4]*matBr[0] + matAr[5]*matBr[1] );
    assert ( matCr[10] == matAr[4]*matBr[4] + matAr[5]*matBr[5] );
    
    delete[] matAr;
    delete[] matBr;
    delete[] matCr;
}


template <typename T>
void test_randSubset(bool verbose=false, uint_t nSamples=20)
{
    if (verbose) std::cout << "Testing randSubset...\n";
    
    scalar_t nv {20};
    vec_t<T> vec (nv, 0);
    for (scalar_t i = 0; i < nv; ++i)
        vec[i] = (T) rand();
    
    if (verbose) {
        std::cout << "superset vector: "; printVector(vec);
    }
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        scalar_t nr = (scalar_t)rand() % nv + 1;
        vec_t<T> randvec = randSubset(vec, nr);
        if (verbose) {
            std::cout << nr << "-sample of superset vector: ";
            printVector(randvec);
        }
        assert ( randvec.size() == nr );
        for (scalar_t i = 0; i < nr; ++i) {
            assert ( isIn(vec, randvec[i]) );
            for (scalar_t j = 0; j < nr; ++j)
                assert ( (j==i) || randvec[j] != randvec[i] );
        }
    }
}

void test_genRandFloatInRange(bool verbose=false, uint_t nSamples=5)
{
    if (verbose) std::cout << "Testing genRandFloatInRange...\n";
    
    uint_t N = 100;
    float minFloat = std::numeric_limits<float>::lowest();
    float maxFloat = std::numeric_limits<float>::max();

    for (uint_t s = 0; s < nSamples; ++s) {
        double sum {0.0};
        for (uint_t n = 0; n < N; ++n) {
            float r = genRandFloatInRange();
            sum += r;
        }
        float avg = (float)(sum / N);
        assert ( avg < 0.25*maxFloat - 0.25*minFloat );
        assert ( avg > - 0.25*maxFloat + 0.25*minFloat );
    }

    uint_t maxUint = std::numeric_limits<uint_t>::max();
    float min {0}, max{0};
    for (uint_t s = 0; s < nSamples; ++s) {
        min = genRandScalar<float>();
        double frac = (double)((uint_t)(rand() << 1)) / (double)maxUint;
        max = (float)(min + frac * ((double)maxFloat - (double)min));
        double sum {0.0};
        for (uint_t n = 0; n < N; ++n) {
            float r = genRandFloatInRange(min, max);
            sum += r;
            assert ( r <= max && r >= min );
        }
        double avg = sum / N;
        double expectation = (float)(((double)max - (double)min)/2.0 + min);
        assert ( avg < expectation + 0.25*max - 0.25*min );
        assert ( avg > expectation - 0.25*max + 0.25*min );
    }
}

void test_vectorProduct(bool verbose=false, uint_t nSamples=200)
{
    if (verbose)
        std::cout << "Testing vectorProduct...\n";
    
    scalar_t maxLength {15};
    for (scalar_t s = 0; s < nSamples; ++s) {
        scalar_t length = ((scalar_t)rand() % maxLength);
        vec_t<long int> vec (length, 0);
        for (scalar_t d = 0; d < length; ++d) {
            vec[d] = rand() % 10;
        }
        // vectors v of length 0 should have vectorProduct(v) == 1
        if (length == 0)
            assert ( vectorProduct(vec) == 1 );
        //assert ...
    }
}


template <typename T>
void assertFactorsEqual(Factor<T> f1, Factor<T> f2)
{
    assert ( f2.label == f1.label );
    assert ( isSuperset(f1.getVars(), f2.getVars()) );
    assert ( isSuperset(f2.getVars(), f1.getVars()) );
    assert ( f2.getVolume() == f1.getVolume() );
    long_uint_t vol = f1.getVolume();
    for (long_uint_t k = 0; k < vol; ++k) {
        assert ( f1.valueAt(k) == f2.valueAt(k) );
    }
}

template <typename T>
void assertEqualUpTillPerm(Factor<T> f1, Factor<T> f2,
                           bool requireLabelsEq=false,
                           double tol=-1)
{
    if (requireLabelsEq)
        assert ( f2.label == f1.label );
    assert ( isSuperset(f1.getVars(), f2.getVars()) );
    assert ( isSuperset(f2.getVars(), f1.getVars()) );
    assert ( f2.getVolume() == f1.getVolume() );
    Factor<T> f1Permd = permuteFactor<T>(f2.getVars(), f1);
    long_uint_t vol = f1Permd.getVolume();
    if (tol <= 0.0) {
        for (long_uint_t k = 0; k < vol; ++k)
            assert ( f1Permd.valueAt(k) == f2.valueAt(k) );
    }
    else {
        for (long_uint_t k = 0; k < vol; ++k) {
            assert ( f1Permd.valueAt(k) + tol >= f2.valueAt(k) );
            assert ( f1Permd.valueAt(k) - tol <= f2.valueAt(k) );
        }
    }
}

void test_sumOverVar_fixed_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing sumOverVar with fixed examples...\n";
    
    // Example 1
    
    vec_t<Var> vars {{4441,2}, {4442,3}, {4443,4}};
    int* arr = new int[24] {
        111,112,113,114,
        121,122,123,124,
        131,132,133,134,
        
        211,212,213,214,
        221,222,223,224,
        231,232,233,234
    };
    
    Var sumVar = vars[1];
    
    Factor<int> f(14831, vars, arr);
    if (verbose) {
        std::cout << "Original factor: \n";
        printFactor(f);
        std::cout << "Var over which to sum: " << sumVar << "\n";
    }
    delete[] arr;
    
    Factor<int> g = sumOverVar(sumVar, f, 54545);
    if (verbose) {
        std::cout << "resulting factor: \n";
        printFactor(g);
    }
    assert ( g.valueAt(0) == 363 );
    assert ( g.valueAt(1) == 366 );
    assert ( g.valueAt(2) == 369 );
    assert ( g.valueAt(3) == 372 );
    assert ( g.valueAt(4+0) == 663 );
    assert ( g.valueAt(4+1) == 666 );
    assert ( g.valueAt(4+2) == 669 );
    assert ( g.valueAt(4+3) == 672 );
    
    // Example 2
    
    sumVar = vars[0];
    Factor<int> h = sumOverVar(sumVar, f, 54545);
    if (verbose) {
        std::cout << "Var over which to sum: " << sumVar << "\n";
        std::cout << "resulting factor: \n";
        printFactor(h);
    }
    assert ( h.valueAt(0) == 322 );
    assert ( h.valueAt(1) == 324 );
    assert ( h.valueAt(2) == 326 );
    assert ( h.valueAt(3) == 328 );
    assert ( h.valueAt(4+0) == 342 );
    assert ( h.valueAt(4+1) == 344 );
    assert ( h.valueAt(4+2) == 346 );
    assert ( h.valueAt(4+3) == 348 );
    assert ( h.valueAt(2*4+0) == 362 );
    assert ( h.valueAt(2*4+1) == 364 );
    assert ( h.valueAt(2*4+2) == 366 );
    assert ( h.valueAt(2*4+3) == 368 );
}


template <typename T>
void test_contractFactors_fixed_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing contractFactors with fixed examples...\n";

    // Example 1
    
    label_t labelF {1234}, labelG {5678};
    
    vec_t<Var> varsF {{12,2}, {14,2}, {11,2}, {13,2}};
    Zp<MP67>* dataF = new Zp<MP67>[2*2*2*2] {
        30,17, 42,2, 14,9, 24,58, 10,6, 59,55, 66,62, 53,47
    };
    Factor<Zp<MP67>> f1 (labelF, varsF, dataF);
    delete[] dataF;
    
    vec_t<Var> varsG {{13,2}, {11,2}, {14,2}, {15,2}};
    Zp<MP67>* dataG = new Zp<MP67>[2*2*2*2] {
        59,17, 41,61, 27,62, 9,16, 66,56, 66,55, 25,43, 26,50
    };
    Factor<Zp<MP67>> g1 (labelG, varsG, dataG);
    delete[] dataG;
    
    vec_t<Var> internVars {{13,2}, {14,2}, {15,2}};
    Factor<Zp<MP67>> result = contractFactors(f1, g1, internVars, 9999);
    
    if(verbose) {
        std::cout << "Factors to be contracted:\n";
        printFactor(f1); printFactor(g1);
        std::cout << "Contraction-internal vars: ";
        printVector(internVars);
        std::cout << "Resulting Factor:\n";
        printFactor(result);
    }

    Zp<MP67>* dataExpected = new Zp<MP67>[2*2] {37,48, 38,19};
    for (uint_t i = 0; i < 2*2; ++i) {
        assert ( result.valueAt(i) == dataExpected[i] );
    }
    delete[] dataExpected;
    
    // Example 2
    
    vec_t<Var> varsF2 {{14,3},{11,2},{13,3}};
    Zp<MP>* dataF2 = new Zp<MP>[3*2*3] {
        529904427, 596849873, 1660430184,
        156526728, 814132196, 29794095,

        33611298, 1266736565, 823726943,
        1517811822, 607214211, 439451303,

        1622159876, 2066965198, 774630728,
        1776201090, 1471416373, 289450985
    };
    Factor<Zp<MP>> f2 (labelF, varsF2, dataF2);
    delete[] dataF2;
    
    vec_t<Var> varsG2 {{14,3}, {11,2}, {12,2}};
    Zp<MP>* dataG2 = new Zp<MP>[3*2*2] {
        1386982170, 1572912529,
        1755285402, 994264123,
        
        1098839211, 1225456759,
        1447161273, 1628743639,

        1822306633, 960107809,
        1785270367, 488955181    
    };
    Factor<Zp<MP>> g2 (labelG, varsG2, dataG2);
    delete[] dataG2;
    
    vec_t<Var> internVars2 {{12,2}, {13,3}, {11,2}};
    Factor<Zp<MP>> result2 = contractFactors(f2, g2, internVars2, 9999);
    
    if(verbose) {
        std::cout << "Factors to be contracted:\n";
        printFactor(f2); printFactor(g2);
        std::cout << "Contraction-internal vars: ";
        printVector(internVars2);
        std::cout << "Resulting Factor:\n";
        printFactor(result2);
    }

    assert ( result2.nVars() == 1 );
    assert ( result2.vars[0].label == 14 );
    assert ( result2.vars[0].size == 3 );
    
    Zp<MP>* dataExpected2 = new Zp<MP>[3] {
        Zp<MP>(1686714793), Zp<MP>(1628096808), Zp<MP>(971519856)
    };
    for (uint_t i = 0; i < 3; ++i) {
        assert ( result2.valueAt(i) == dataExpected2[i]);
    }
    delete[] dataExpected2;
}

template <typename T>
void print_contractFactors_random()
{
    std::cout << "Generating two random Factors and "
        << "printing out results of contracting them...\n";
    
    vec_t<label_t> varLabels {11,12,13,14,15,16};
    scalar_t maxNVars {5}, maxVarSize {4};
    
    std::pair<Factor<T>, map_t<label_t,scalar_t>>
    fp = genRandFactor<T>(maxNVars, maxVarSize, varLabels, {}, 1);
    Factor<T> f = fp.first;
    
    std::pair<Factor<T>, map_t<label_t,scalar_t>>
    gp = genRandFactor<T>(
        maxNVars, maxVarSize, varLabels, fp.second, 1
    );
    Factor<T> g = gp.first;
    
    scalar_t nInternVars = (scalar_t)rand() % (varLabels.size() + 1);
    vec_t<label_t> internVarIds = randSubset(varLabels, nInternVars);
    
    std::cout << "Factors to be contracted:\n";
    printFactor(f); printFactor(g);
    std::cout << "Contraction-internal varIds: ";
    printVector(internVarIds);
    
    std::cout << "Resulting Factor:\n";
    Factor<T> result = contractFactors(f, g, internVarIds, 99992);
    printFactor(result);
}


template <typename T>
void test_contractFactors_edge_cases(bool verbose=false, uint_t nSamples=10)
{
    if (verbose) {
        std::cout << "Testing some aspects of a few edge cases of "
            << "contractFactors...\n";
    }
    
    if (verbose)
        std::cout << "Case: factors have no common variables/dimensions\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        
        vec_t<label_t> varLabelsF {11,12,13,14,15};
        scalar_t maxNVars {5}, maxVarSize {4};
        
        std::pair<Factor<T>, map_t<label_t,scalar_t>>
        fp = genRandFactor<T>(
            maxNVars, maxVarSize, varLabelsF, {}, 1
        );
        Factor<T> f = fp.first;
        
        vec_t<label_t> varLabelsG {16,17,18,19,20};
        std::pair<Factor<T>, map_t<label_t,scalar_t>>
        gp = genRandFactor<T>(
            maxNVars, maxVarSize, varLabelsG, fp.second, 1
        );
        Factor<T> g = gp.first;
        
        vec_t<label_t> varLabelsAll = varLabelsF;
        varLabelsAll.insert(varLabelsAll.end(),
                            varLabelsG.begin(), varLabelsG.end());
        scalar_t nInternVars = (scalar_t)rand() % (varLabelsAll.size() + 1);
        vec_t<label_t> internVarIds = randSubset(
            varLabelsAll, nInternVars);
        
        if (verbose) {
            std::cout << "Factors to be contracted:\n";
            printFactor(f, false); printFactor(g, false);
            std::cout << "Contraction-internal varIds: ";
            printVector(internVarIds);
        }
        
        Factor<T> result = contractFactors(f, g, internVarIds, 99992);
        
        if (verbose) {
            std::cout << "Resulting Factor:\n";
            printFactor(result, false);
        }
        
        long_scalar_t volExpected {1};
        for (scalar_t d = 0; d < f.nVars(); ++d) {
            if ( ! isIn(internVarIds, f.vars[d].label))
                volExpected *= f.vars[d].size;
        }
        for (scalar_t d = 0; d < g.nVars(); ++d) {
            if (!isIn(internVarIds, g.vars[d].label))
                volExpected *= g.vars[d].size;
        }
        assert ( result.volume == volExpected );
    }
    
    if (verbose)
        std::cout << "Case: all vars are internal to the contraction\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        
        vec_t<label_t> varLabels {11,12,13,14,15,16};
        scalar_t maxNVars {5}, maxVarSize {4};
        
        std::pair<Factor<T>, map_t<label_t,scalar_t>>
        fp = genRandFactor<T>(
            maxNVars, maxVarSize, varLabels, map_t<label_t,scalar_t>(), 1
        );
        Factor<T> f = fp.first;
        
        std::pair<Factor<T>, map_t<label_t,scalar_t>>
        gp = genRandFactor<T>(
            maxNVars, maxVarSize, varLabels, fp.second, 1
        );
        Factor<T> g = gp.first;
        
        vec_t<label_t> internVarIds = varLabels;
        
        if (verbose) {
            std::cout << "Factors to be contracted:\n";
            printFactor(f); printFactor(g);
            std::cout << "Contraction-internal varIds: ";
            printVector(internVarIds);
        }
        
        Factor<T> result = contractFactors(f, g, internVarIds, 99992);
        if (verbose) {
            std::cout << "Resulting Factor:\n";
            printFactor(result);
        }
        
        assert ( result.volume == 1 );
        assert ( result.data );
        assert ( result.nVars() == 0 );
        assert ( result.nVars() == 0 );
    }
    
    if (verbose)
        std::cout << "Case: one of the Factors represents a scalar\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        
        vec_t<label_t> varLabels {11,12,13,14,15,16};
        scalar_t maxNVars {5}, maxVarSize {4};
        
        std::pair<Factor<T>, map_t<label_t,scalar_t>>
        fp = genRandFactor<T>(
            maxNVars, maxVarSize, varLabels, map_t<label_t,scalar_t>(), 1
        );
        Factor<T> f = fp.first;
        
        T scalarVal = (T) (rand() % 1000);
        T* dataG = new T[1] {scalarVal};
        Factor<T> g(564142, {}, dataG);
        delete[] dataG;
        
        scalar_t nInternVars = (scalar_t)rand() % (varLabels.size() + 1);
        vec_t<label_t> internVarIds = randSubset(varLabels, nInternVars);
        
        if (verbose) {
            std::cout << "Factors to be contracted:\n";
            printFactor(f); printFactor(g);
            std::cout << "Contraction-internal varIds: ";
            printVector(internVarIds);
        }
        
        Factor<T> result = contractFactors(f, g, internVarIds, 99992);
        if (verbose) {
            std::cout << "Resulting Factor:\n";
            printFactor(result);
        }
        
        vec_t<label_t> internsF = intersect(f.varIds(), internVarIds);
        Factor<T> fs = f;
        for (label_t vId : internsF)
            fs = sumOverVar(vId, fs, fs.label);
        
        for (long_scalar_t k = 0; k < result.volume; ++k) {
            assert ( result.valueAt(k) == fs.valueAt(k) * scalarVal );
        }
    }
}

template <typename Z>
void test_Factor_slice(bool verbose=false, uint_t nSamples=20)
{
    if (verbose)
        std::cout << "Testing Factor.slice ...\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        vec_t<label_t> varLabels {11,12,13,14,15,16,17};
        Factor<Z> f = genRandFactor<Z>(6, 5, varLabels, 2);
        if (verbose) {
            std::cout << "Factor to be sliced:\n"; printFactor(f);
        }
        scalar_t nSliceVars = (scalar_t)rand() % (f.nVars() + 1);
        vec_t<Var> sliceVars = randSubset(f.vars, nSliceVars);
        map_t<Var, scalar_t> varsToInds;
        for (Var v : sliceVars) {
            uint_t sliceInd = (uint_t)rand() % v.size;
            varsToInds[v] = sliceInd;
        }
        if (verbose) {
            std::cout << "slice parameter varsToInds:\n"; printMap(varsToInds);
        }
        Factor<Z> g = f.slice(varsToInds, 999);
        
        assert ( g.nVars() == f.nVars() - nSliceVars );
        long_scalar_t expectedVol = 1;
        for (Var v : f.vars) {
            if (isIn(sliceVars, v))
                assert ( ! isIn(g.vars, v) );
            else {
                assert ( isIn(g.vars, v) );
                expectedVol *= v.size;
            }
        }
        assert ( g.volume == expectedVol );
    }
    
}

template <typename T>
void test_FactorGraph_add_remove(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing FactorGraph.insertFactor and .removeFactor...\n";
    
    scalar_t nFactors {10};
    vec_t<label_t> varLabels {11,12,13,14,15,16,17,18,19};
    scalar_t maxNVars {5}, maxVarSize {4};
    FactorGraph<int> fg = genRandFactorGraph<T>(
        nFactors, maxNVars, maxVarSize, varLabels, false, nFactors);
    
    assert ( fg.nFactors() == nFactors );
    for (label_t n = 0; n < nFactors; ++n) {
        assert ( isIn(fg.getFactorIds(), n) );
        for (Var v : fg.getFactorWithId(n).vars) {
            assert ( isIn(fg.getVars(), v) );
        }
    }
    for (label_t n = 0; n < nFactors; ++n) {
        fg.removeFactor(n);
        assert ( !(isIn(fg.getFactorIds(), n)) );
    }
    assert ( fg.nFactors() == 0 );
}

template <typename T>
void test_FactorGraph_getInternalVarIds_examples(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing FactorGraph.getInternalVarIds with "
            << "a few fixed examples...\n";
    
    vec_t<T> dataA(2*2*2*2, (T)rand());
    Factor<T> a (7111, {{11,2},{12,2},{13,2},{14,2}}, dataA.data());
    vec_t<T> dataB(2*2*2*2, (T)rand());
    Factor<T> b (7222, {{11,2},{12,2},{13,2},{15,2}}, dataB.data());
    vec_t<T> dataC(2, (T)rand());
    Factor<T> c (7333, {{14,2}}, dataC.data());
    
    set_t<Var> boundary {{13,2},{15,2}};
    set_t<Factor<T>> factors {a,b,c};
    FactorGraph<T> fg(factors, boundary);
    
    vec_t<label_t> internVarIds = fg.getInternalVarIds(7111, 7222);
    if (verbose) {
        std::cout << "FactorGraph fg consists of following Factors:\n";
        printFactor(a); printFactor(b); printFactor(c);
        std::cout << "boundary: "; printVector(setToVector(boundary));
        std::cout << "fg.getInternalVarIds(7111, 7222): ";
        printVector(internVarIds);
    }
    assert ( isIn(internVarIds, (label_t)11) );
    assert ( isIn(internVarIds, (label_t)12) );
    assert ( ! isIn(internVarIds, (label_t)13) );
    assert ( ! isIn(internVarIds, (label_t)14) );
    assert ( ! isIn(internVarIds, (label_t)15) );
}

template <typename T>
void test_FactorGraph_contractInOrder(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing FactorGraph.contractInOrder with "
            << "a random example...\n";
    
    scalar_t nFactors {7};
    vec_t<label_t> varLabels {11,12,13,14,15,16,17,18,19};
    scalar_t maxNVars {5}, maxVarSize {4};
    FactorGraph<int> fg = genRandFactorGraph<T>(
        nFactors, maxNVars, maxVarSize, varLabels, false, nFactors);
    scalar_t nBoundaryVars {3};
    set_t<Var> boundary = vectorToSet(randSubset(fg.getVars(), nBoundaryVars));
    fg.setBoundary(boundary);
    
    scalar_t nVarsBefore = fg.nVars();
    
    vec_t<IdTriple> ordering {
        {0, 1, 100}, {2, 3, 101}, {100, 101, 103}
    };
    fg = fg.contractInOrder(ordering);
    vec_t<label_t> remainingFactorIds = fg.getFactorIds();
    assert ( isIn(remainingFactorIds, (label_t)4) );
    assert ( isIn(remainingFactorIds, (label_t)5) );
    assert ( isIn(remainingFactorIds, (label_t)6) );
    assert ( ! isIn(remainingFactorIds, (label_t)0) );
    assert ( ! isIn(remainingFactorIds, (label_t)1) );
    assert ( ! isIn(remainingFactorIds, (label_t)2) );
    assert ( ! isIn(remainingFactorIds, (label_t)3) );
    assert ( ! isIn(remainingFactorIds, (label_t)100) );
    assert ( ! isIn(remainingFactorIds, (label_t)101) );
    assert ( isIn(remainingFactorIds, (label_t)103) );
    
    assert ( fg.nVars() <= nVarsBefore );
    
    ordering = {
        {4, 103, 0}, {5, 0, 1}, {6, 1, 0}
    };
    fg = fg.contractInOrder(ordering);
    assert ( isIn(fg.getFactorIds(), (label_t)0) );
    for (label_t i = 1; i < 109; ++i)
        assert ( ! isIn(fg.getFactorIds(), (label_t)i) );
    
    for (Var bvar: boundary)
        assert ( isIn(fg.getVars(), bvar) );
    
    assert ( fg.nVars() == nBoundaryVars );
    assert ( fg.nFactors() == 1 );
}

template <typename T>
void test_FactorGraph_contractAll(bool verbose=false, uint_t nSamples=100)
{
    if (verbose)
        std::cout << "Testing FactorGraph.contractAll with "
            << "semi-random examples...\n";
        
    for (scalar_t s = 0; s < nSamples; ++s) {
        scalar_t maxNFactors {10};
        vec_t<label_t> varLabels {11,12,13,14,15,16,17,18,19,20};
        scalar_t maxNVars {5}, maxVarSize {4};
        FactorGraph<int> fg = genRandFactorGraph<T>(
            maxNFactors, maxNVars, maxVarSize, varLabels);
        scalar_t nBoundaryVars = (scalar_t)rand() % (fg.nVars() + 1);
        set_t<Var> boundary = randSubset(
            vectorToSet(fg.getVars()), nBoundaryVars);
        fg.setBoundary(boundary);
                
        fg = fg.contractAll();
        
        assert ( fg.nVars() == nBoundaryVars );
        assert ( isSuperset(boundary, fg.getBoundary()) );
        assert ( fg.nFactors() == 1 );
    }
}

template <typename T>
void test_FactorGraph_slice(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing FactorGraph.slice...\n";
        
    for (scalar_t s = 0; s < nSamples; ++s) {
        vec_t<label_t> varLabels {11,12,13,14,15,16,17,18,19};
        FactorGraph<T> fg = genRandFactorGraph<T>(8, 5, 5, varLabels);
        uint_t nSliceVars = (uint_t)rand() % (fg.nVars() + 1);
        vec_t<Var> sliceVars = randSubset(fg.getVars(), nSliceVars);
        map_t<Var, scalar_t> varsToInds;
        for (Var v : sliceVars) {
            uint_t sliceInd = (label_t)rand() % v.size;
            varsToInds[v] = sliceInd;
        }
        FactorGraph<T> fgs = fg.slice(varsToInds);
        
        assert ( fgs.nFactors() == fg.nFactors() );
        assert ( fgs.nVars() == fg.nVars() - nSliceVars );
    }
}

template <typename F>
void test_contractionCost(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing contractionCost (with fixed examples)...\n";
    
    Var v1(11, 2), v2(12, 3), v3(13, 4);
    std::shared_ptr<F[]> dummyData;
    Factor<F> f1(1, {v1,v2}, dummyData);
    FactorGraph<F> fg({f1}, {});
    vec_t<IdTriple> order {};
    assert ( contractionCost(fg, order) == 0 );
    
    Factor<F> f2(2, {v2,v3}, dummyData);
    fg = FactorGraph<F>({f1,f2}, {});
    order = {{1,2,100}};
    assert ( contractionCost(fg, order) == 24 );
    
    fg.setBoundary({v1});
    assert ( contractionCost(fg, order) == 24 );
    fg.setBoundary({v2});
    assert ( contractionCost(fg, order) == 24 );
    
    Var v4(14,5), v5(15,6);
    f2 = Factor<F>(2, {v2,v3,v4,v5}, dummyData);
    Factor<F> f3(3, {v4,v5}, dummyData);
    fg = FactorGraph<F>({f1,f2,f3}, {});
    order = {{1,2,100}, {100,3,101}};
    assert ( contractionCost(fg, order) == 2*3*4*5*6 + 5*6 );
    
    order = {{2,3,100}, {100,1,101}};
    assert ( contractionCost(fg, order) == 3*4*5*6 + 2*3 );
    
    order = {{1,3,100}, {100,2,101}};
    assert ( contractionCost(fg, order) == 2*3*5*6 + 3*4*5*6 );
    
    fg.setBoundary({v3});
    order = {{1,2,100}, {100,3,101}};
    assert ( contractionCost(fg, order) == 2*3*4*5*6 + 4*5*6 );
}

bool isContractionOf(label_t a, label_t b, IdTriple c)
{
    return ( (c.first == a && c.second == b)
              || (c.first == b || c.second == a) );
}

template <typename F>
void test_optimalContractionOrder(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing optimalContractionOrder "
                  << "(with fixed examples)...\n";
    
    Var v1(11,2), v2(12,2), v3(13,2), v4(14,2), v5(15,2), v6(16,2);
    std::shared_ptr<F[]> dummyData;
    Factor<F> f1(1, {v1,v2,v3}, dummyData);
    Factor<F> f2(2, {v1,v2,v3,v4}, dummyData);
    Factor<F> f3(3, {v4,v5,v6}, dummyData);
    Factor<F> f4(4, {v5,v6}, dummyData);
    FactorGraph<F> fg({f1,f2,f3,f4}, {});
    map_t< size_t, long_uint_t > costs;
    map_t< size_t, vec_t<IdTriple> > orders;
    vec_t<label_t> fIds {1,2,3,4};
    
    vec_t<IdTriple> order = optimalContractionOrder(fg, fIds, costs, orders);
    assert ( order.size() == 3 );
    assert ( isContractionOf(1, 2, order[0]) );
    assert ( isContractionOf(3, 4, order[1]) );
    assert ( costs[labelSetHash(fIds)] == (2*2*2)*2 + (2*2)*2 + 2 );
    
    fIds = {2,3,4};
    costs = {}; orders = {};
    order = optimalContractionOrder(fg, fIds, costs, orders);
    assert ( order.size() == 2 );
    assert ( isContractionOf(3, 4, order[0]) );
    assert ( isContractionOf(2, order[0].newLabel, order[1]) );
    assert ( order[0].first != 1 && order[0].second != 1
             && order[1].first != 1 && order[1].second != 1 );
    assert ( costs[labelSetHash(fIds)] == (2*2)*2 + 2*(2*2*2) );
    
    fIds = {1,2,3,4};
    costs = {}; orders = {};
    fg.setBoundary({v1,v2,v3});
    order = optimalContractionOrder(fg, fIds, costs, orders);
    assert ( order.size() == 3 );
    assert ( isContractionOf(3, 4, order[0]) );
    assert ( isContractionOf(2, order[0].newLabel, order[1]) );
    assert ( isContractionOf(1, order[1].newLabel, order[2]) );
    assert ( costs[labelSetHash(fIds)] == (2*2)*2 + 2*(2*2*2) + 2*2*2 );
    
}


template <typename Z>
void test_Poly_rand(bool verbose=false, uint_t nSamples=200)
{
    if (verbose)
        std::cout << "Testing Poly<" << typeid(Z).name() << ">::rand...\n";
    
    scalar_t maxNTerms {8};
    scalar_t nTerms {1};
    long_scalar_t acc {0};
    Z maxZ = Z::zero - Z::one;
    for (scalar_t s = 0; s < nSamples; ++s) {
        nTerms = (scalar_t)rand() % maxNTerms + 1;
        Poly<Z> p = Poly<Z>::rand(nTerms);
        assert ( p.degree() <= nTerms - 1 );
            // p.degree() < (nTerms - 1) if lead coeff is Z::zero
        for (scalar_t i = 0; i < nTerms; ++i) {
            assert ( p[i] >= Z::zero );
            assert ( p[i] <= maxZ );
            acc += p[i].raw();
        }
    }
    long_scalar_t modulus = maxZ.value() + 1;
    long_scalar_t avg = acc / (nSamples * (maxNTerms / 2));
    assert ( avg > 0.3 * modulus );
    assert ( avg < 0.7 * modulus );
}


template <typename Z>
void test_Poly_eval(bool verbose=false, uint_t nSamples=20)
{
    if (verbose)
        std::cout << "Testing Poly<" << typeid(Z).name()
            << ">::operator() ...\n";
    
    index_t maxNTerms {20};
    index_t nTerms {1};
    for (scalar_t s = 0; s < nSamples; ++s) {
        nTerms = rand() % maxNTerms + 1;
        Poly<Z> p = Poly<Z>::rand(nTerms);
        Z x = Z::rand();
        Z y = p(x);
        Z expectedY = p[nTerms-1];
        for (index_t i = nTerms-2; i >= 0; --i) {
            expectedY = expectedY * x + p[i];
        }
        assert ( y == expectedY );
    }
}


template <typename Z>
void test_Poly_arithmetic_plus(bool verbose=false, uint_t nSamples=20)
{
    if (verbose)
        std::cout << "Testing Poly<" << typeid(Z).name()
            << ">::operator+ and Poly::operator- ...\n";
    
    index_t maxNTerms {20};
    index_t nTermsA {1}, nTermsB {1};
    for (scalar_t s = 0; s < nSamples; ++s) {
        nTermsA = rand() % maxNTerms + 1;
        nTermsB = rand() % maxNTerms + 1;
        Poly<Z> a = Poly<Z>::rand(nTermsA);
        Poly<Z> b = Poly<Z>::rand(nTermsB);
        Poly<Z> cSum = a + b;
        Poly<Z> cDiff = a - b;
        index_t ad = a.degree();
        index_t bd = b.degree();
        index_t deg = ad > bd ? ad : bd;
        for (index_t i = 0; i < deg; ++i) {
            Z aCoeffI = ad < i ? Z::zero : a[i];
            Z bCoeffI = bd < i ? Z::zero : b[i];
            assert ( cSum[i] == aCoeffI + bCoeffI );
            assert ( cDiff[i] == aCoeffI - bCoeffI );
        }
        Z r = Z::rand();
        assert ( cSum(r) == a(r) + b(r) );
        assert ( cDiff(r) == a(r) - b(r) );
    }
}

template <typename Z>
void test_Poly_arithmetic_mul(bool verbose=false, uint_t nSamples=100)
{
    if (verbose)
        std::cout << "Testing Poly<" << typeid(Z).name()
            << ">::operator* ...\n";
    
    index_t maxNTerms {10};
    index_t nTermsA {1}, nTermsB {1};
    for (scalar_t s = 0; s < nSamples; ++s) {
        nTermsA = rand() % maxNTerms + 1;
        nTermsB = rand() % maxNTerms + 1;
        Poly<Z> a = Poly<Z>::rand(nTermsA);
        Poly<Z> b = Poly<Z>::rand(nTermsB);
        Poly<Z> c = a * b;
        index_t ad = a.degree();
        index_t bd = b.degree();
        if (ad == -1 || bd == -1)  // a is all zeros || b is all zeros
            assert ( c.degree() == -1 );
        else
            assert ( c.degree() == ad + bd );
        Z r = Z::rand();
        assert ( c(r) == a(r) * b(r) );
    }
}

template <typename Z>
void test_Poly_arithmetic_quorem(bool verbose=false, uint_t nSamples=100)
{
    if (verbose)
        std::cout << "Testing Poly<" << typeid(Z).name()
            << ">::operator/ and Poly::operator% ...\n";
    
    index_t maxNTerms {6};
    index_t nTermsA {1}, nTermsB {1};
    for (scalar_t s = 0; s < nSamples; ++s) {
        nTermsA = rand() % maxNTerms + 1;
        nTermsB = rand() % maxNTerms + 1;
        Poly<Z> a = Poly<Z>::rand(nTermsA);
        Poly<Z> b = Poly<Z>::rand(nTermsB);
        while (b.degree() < 0)  // cannot divide by zero-polynomial
            b = Poly<Z>::rand(nTermsB);
        Poly<Z> quo = a / b;
        Poly<Z> rem = a % b;
        index_t ad = a.degree();
        index_t bd = b.degree();
        /*
        std::cout << "a: " << a << "\n";
        std::cout << "b: " << b << "\n";
        std::cout << "quo: " << quo << "\n";
        std::cout << "rem: " << rem << "\n";
        std::cout << "a.degree(): " << a.degree() << "\n";
        std::cout << "b.degree(): " << b.degree() << "\n";
        std::cout << "quo.degree(): " << quo.degree() << "\n";
        std::cout << "rem.degree(): " << rem.degree() << "\n";
        */
        assert ( quo.degree() == (ad >= bd ? ad - bd : -1) );
        assert ( rem.degree() < bd );
        Z r = Z::rand();
        assert ( quo(r)*b(r) + rem(r) == a(r) );
    }
}

template <typename Z>
void test_Poly_interpolate(bool verbose=false, uint_t nSamples=200)
{
    if (verbose)
        std::cout << "Testing Poly<" << typeid(Z).name()
            << ">::interpolate ...\n";
    
    index_t maxNTerms {20};
    index_t nTerms {1};
    for (scalar_t s = 0; s < nSamples; ++s) {
        nTerms = rand() % maxNTerms + 1;
        Poly<Z> p = Poly<Z>::rand(nTerms);
        vec_t<Z> inputs;
        vec_t<Z> targets;
        for (scalar_t i = 0; i < nTerms; ++i) {
            Z u = Z::rand();
            while (isIn(inputs, u))
                u = Z::rand();
            inputs.push_back(u);
            targets.push_back(p(u));
        }
        Poly<Z> q = Poly<Z>::interpolate(
            nTerms, inputs.data(), targets.data()
        );
        assert ( q.degree() <= nTerms - 1 );
        assert ( q.degree() == p.degree() );
        for (index_t i = 0; i < nTerms; ++i) {
            assert ( q[i] == p[i] );
        }
        Z r = Z::rand();
        assert ( q(r) == p(r) );
    }
}


template <typename Z>
void test_Poly_batchEval(bool verbose=false, uint_t nSamples=100)
{
    if (verbose)
        std::cout << "Testing Poly<" << typeid(Z).name()
            << ">::batchEval ...\n";
    
    index_t maxNTerms {20};
    index_t nTerms = rand() % maxNTerms + 1;
    Poly<Z> p = Poly<Z>::rand(nTerms);
    Z* inputs = new Z[nSamples];
    for (scalar_t i = 0; i < nSamples; ++i)
        inputs[i] = Z::rand();
    Z* evals = new Z[nSamples];
    
    p.batchEval(nSamples, inputs, evals);
    
    for (index_t i = 0; i < nSamples; ++i)
        assert ( p(inputs[i]) == evals[i] );
    
    delete[] inputs;
    delete[] evals;
}

template <typename Z>
void test_rsDecodeXY(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing rsDecodeXY<" << typeid(Z).name() << ">...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        index_t maxDeg {20};
        index_t deg = rand() % maxDeg + 1;
        index_t e = deg + 1 + (rand() % 20);
        Z* inputs = new Z[e];
        for (scalar_t i = 0; i < e; ++i)
            inputs[i] = Z(i);
        
        Poly<Z> p = Poly<Z>::rand(deg);
        Z* evals = new Z[e];
        p.batchEval(e, inputs, evals);
        
        corrupt<Z>(deg, e, evals);
        
        Poly<Z> decodedP;
        bool success = rsDecodeXY(deg, e, inputs, evals, decodedP);
        assert ( success );
        
        for (index_t i = 0; i < e; ++i)
            assert ( p(inputs[i]) == decodedP(inputs[i]) );
        
        delete[] inputs;
        delete[] evals;
    }
}

template <typename Z>
void test_invertVandermonde(bool verbose=false, uint_t nSamples=200)
{
    if (verbose)
        std::cout << "Testing invertVandermonde<" << typeid(Z).name()
            << "> ...\n";
    
    index_t maxN {10};
    index_t n {1};
    for (scalar_t s = 0; s < nSamples; ++s) {
        n = rand() % maxN + 2;
        Z* vm = new Z[n*n];
        genRandVandermonde(n, vm);
        Z* vminv = new Z[n*n];
        invertVandermonde(n, vm, vminv, true);
        Z* product = new Z[n*n];
        naiveMatrixMul(n,n,n, vm, vminv, product);
        /*if (verbose) {
            scalar_t dimSizes[2] {(scalar_t)n, (scalar_t)n};
            std::cout << "Vandermonde matrix: \n";
            printNDArray(2, dimSizes, vm);
            std::cout << "Transpose of inverse Vandermonde matrix: \n";
            printNDArray(2, dimSizes, vminv);
            std::cout << "Product of the two above: \n";
            printNDArray(2, dimSizes, product);
        }*/
        for (index_t r = 0; r < n; ++r) {
            for (index_t c = 0; c < n; ++c) {
                if (r == c)
                    assert ( product[n*r + c] == Z::one );
                else
                    assert ( product[n*r + c] == Z::zero );
            }
        }
        delete[] vm;
        delete[] vminv;
        delete[] product;
    }
}

void test_enumerateVectors(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing enumerateVectors with a few examples ...\n";
    
    // Example 1
    
    vec_t<uint_t> ranges {2,3,4};
    vec_t<vec_t<int>> result;
    enumerateVectors<int>(ranges, result);
    
    assert ( result.size() == 24 );
    for (vec_t<int> v : result)
        assert ( v.size() == 3 );
    
    assert ( result[0][0] == 0 && result[0][1] == 0 && result[0][2] == 0);
    assert ( result[1][0] == 0 && result[1][1] == 0 && result[1][2] == 1);
    assert ( result[2][0] == 0 && result[2][1] == 0 && result[2][2] == 2);
    assert ( result[3][0] == 0 && result[3][1] == 0 && result[3][2] == 3);
    assert ( result[4][0] == 0 && result[4][1] == 1 && result[4][2] == 0);
    assert ( result[5][0] == 0 && result[5][1] == 1 && result[5][2] == 1);
    // ...
    assert ( result[11][0] == 0 && result[11][1] == 2 && result[11][2] == 3);
    assert ( result[12][0] == 1 && result[12][1] == 0 && result[12][2] == 0);
    assert ( result[13][0] == 1 && result[13][1] == 0 && result[13][2] == 1);
    // ...
    assert ( result[23][0] == 1 && result[23][1] == 2 && result[23][2] == 3);
    
    // Example 2
    
    ranges = {3};
    result = {};
    enumerateVectors<int>(ranges, result);
    
    assert ( result.size() == 3 );
    for (vec_t<int> v : result)
        assert ( v.size() == 1 );
    
    assert ( result[0][0] == 0 );
    assert ( result[1][0] == 1 );
    assert ( result[2][0] == 2 );
    
    // Example 3
    
    ranges = {};
    result = {};
    enumerateVectors<int>(ranges, result);
    
    assert ( result.size() == 0 );
}

// Helper for test_nextCombination
// Translate vector to number in variable base
template <typename T>
T vecToNum(vec_t<T> c, vec_t<T> bases)
{
    uint_t n = c.size();
    assert ( n == bases.size() );
    T result = 0;
    T pow = 1;
    for (uint_t i = 1; i <= n; ++i) {
        assert ( c[n-i] < bases[n-i] || bases[n-i] == 0 );
        result += pow * c[n-i];
        pow *= bases[n-i];
    }
    return result;
}

void test_nextCombination(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing nextCombination\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        uint_t n = 1 + (uint_t)rand() % 7;
        int max = 7;
        vec_t<int> ranges;
        for (uint_t i = 0; i < n; ++i)
            ranges.push_back( 1 + rand() % max );
        vec_t<int> combination(n, 0);
        int count = 0;
        int c0 = -1;
        do {
            count++;
            int c1 = vecToNum(combination, ranges);
            assert ( c1 > c0 );  // each combination is unique
            c0 = c1;
        } while (nextCombination(combination, ranges));

        assert ( count == vectorProduct(ranges) );
    }
}

template <typename Z>
void test_transformVarVandermonde(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing some aspects of transformVarVandermonde<"
            << typeid(Z).name() << "> ...\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        FactorGraph<Z> fg = genRandFactorGraph<Z>(6, 5, 5);
        scalar_t nFactors = fg.nFactors();
        scalar_t nVars = fg.nVars();
        Var v = fg.getVars()[(uint_t)rand() % nVars];
        vec_t<label_t> varNeighbors = fg.adjacentFactorIds(v);
        scalar_t nNeighbors = varNeighbors.size();
        map_t<scalar_t, Z> eta;
        for (scalar_t i = 0; i < v.size; ++i)
            eta.insert({i, Z(10+i)});
        Z evalPoint = Z::rand();
        transformVarVandermonde(fg, v, evalPoint, eta);
        
        assert ( ! isIn(fg.getVars(), v) );
        assert ( fg.nFactors() == nFactors + nNeighbors * 2 );
        assert ( fg.nVars() == nVars - 1 + nNeighbors * 2 );
        assert ( fg.getVarIds().size() == nVars - 1 + nNeighbors * 2 );
    }
}

template <typename Z>
vec_t<Var> randVarSubset(const FactorGraph<Z>& fg,
                         uint_t subsetMaxSize, uint_t subsetMinSize=0)
{
    uint_t subsetSize = subsetMinSize;
    subsetSize += (uint_t)rand() % ((subsetMaxSize - subsetMinSize) + 1);
    subsetSize = subsetSize > fg.nVars() ? fg.nVars() : subsetSize;
    vec_t<Var> subset = randSubset(fg.getVars(), subsetSize);
    std::sort(subset.begin(), subset.end());
    return subset;
}

template <typename Z>
void test_evaluateProofPolynomialVandermonde(bool verbose=false,
                                             uint_t nSamples=2)
{
    if (verbose)
        std::cout << "Testing evaluateProofPolynomialVandermonde<"
            << typeid(Z).name() << "> ...\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        // Construct a random FactorGraph and cutset
        FactorGraph<Z> fg = genRandFactorGraph<Z>(6, 4, 4);
        vec_t<Var> cutset = randVarSubset(fg, 3);
        uint_t cutsetSize = cutset.size();
        uint_t boundarySize = (uint_t)rand() % (cutsetSize + 1);
        vec_t<Var> boundary = randSubset(cutset, boundarySize);
        fg.setBoundary(vectorToSet(boundary));
        
        // Construct injective map eta_i for i in cutset
        map_t<Var, map_t<uint_t,Z>> etas = constructEtas<Z>(cutset);
        
        // Construct injective map tau
        map_t<vec_t<scalar_t>, Z> tau = constructTau<Z>(cutset);

        // For any zeta = tau(w) in the range of tau, it should hold that
        // proofPolynomial(zeta) == contractAll(fg.slice(w))
        for (auto keyVal : tau) {
            vec_t<scalar_t> w = keyVal.first;
            Z zeta = tau.at(w);
            Z result = evaluateProofPolynomialVandermonde(
                fg, cutset, etas, tau, zeta);
            
            map_t<Var, scalar_t> varsToInds;
            for (uint_t i = 0; i < cutsetSize; ++i)
                varsToInds[cutset[i]] = w[i];
            FactorGraph<Z> fgv = fg.slice(varsToInds);
            fgv = fgv.contractAll();
            
            assert ( fgv.nFactors() == 1 );
            Factor<Z> f = fgv.getFactors()[0];
            assert ( f.volume == 1 );
            assert ( result == f.valueAt(0) );
        }
        
        // the proof polynomial should also be evaluable at points
        // zeta that are outside the range of tau
        uint_t nRandPoints {20};
        for (uint_t i = 0; i < nRandPoints; ++i) {
            Z zeta = Z::rand();
            evaluateProofPolynomialVandermonde(
                fg, cutset, etas, tau, zeta);
        }
    }
}

template <typename Z>
void test_evaluateProofPolynomial(bool verbose=false, uint_t nSamples=4)
{
    if (verbose)
        std::cout << "Testing evaluateProofPolynomial<"
            << typeid(Z).name() << "> ...\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        // Construct a random FactorGraph and cutset
        FactorGraph<Z> fg = genRandFactorGraph<Z>(6, 4, 4);
        vec_t<Var> cutset = randVarSubset(fg, 3);
        uint_t cutsetSize = cutset.size();
        uint_t boundarySize = (uint_t)rand() % (cutsetSize + 1);
        vec_t<Var> boundary = randSubset(cutset, boundarySize);
        fg.setBoundary(vectorToSet(boundary));
        
        // Construct injective map eta_i for i in cutset
        map_t<Var, map_t<uint_t,Z>> etas = constructEtas<Z>(cutset);
        
        // Construct injective map tau
        map_t<vec_t<scalar_t>, Z> tau = constructTau<Z>(cutset);
        
        // For any zeta = tau(w) in the range of tau, it should hold that
        // proofPolynomial(zeta) == contractAll(fg.slice(w))
        for (auto keyVal : tau) {
            vec_t<scalar_t> w = keyVal.first;
            Z zeta = tau.at(w);
            Z result = evaluateProofPolynomial(fg, cutset, etas, tau, zeta);
            
            map_t<Var, scalar_t> varsToInds;
            for (uint_t i = 0; i < cutsetSize; ++i)
                varsToInds[cutset[i]] = w[i];
            FactorGraph<Z> fgv = fg.slice(varsToInds);
            fgv = fgv.contractAll();
            
            assert ( fgv.nFactors() == 1 );
            Factor<Z> f = fgv.getFactors()[0];
            assert ( f.volume == 1 );
            assert ( f.isScalar() );
            assert ( result == f.scalarValue() );
        }
        
        // the proof polynomial should also be evaluable at points
        // zeta that are outside the range of tau
        uint_t nRandPoints {20};
        for (uint_t i = 0; i < nRandPoints; ++i) {
            Z zeta = Z::rand();
            evaluateProofPolynomial(fg, cutset, etas, tau, zeta);
        }
    }
}

template <typename Z>
void test_evaluateProofPolynomial_compare(bool verbose=false,
                                          uint_t nSamples=3)
{
    if (verbose)
        std::cout << "Testing evaluateProofPolynomial<"
            << typeid(Z).name() << "> by comparing to outputs of "
            << "evaluateProofPolynomialVandermonde...\n";
    
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        // Construct a random FactorGraph and cutset
        FactorGraph<Z> fg = genRandFactorGraph<Z>(6, 4, 4);
        vec_t<Var> cutset = randVarSubset(fg, 3);
        uint_t cutsetSize = cutset.size();
        uint_t boundarySize = (uint_t)rand() % (cutsetSize + 1);
        vec_t<Var> boundary = randSubset(cutset, boundarySize);
        fg.setBoundary(vectorToSet(boundary));
        
        // Construct injective map eta_i for i in cutset
        map_t<Var, map_t<uint_t,Z>> etas = constructEtas<Z>(cutset);
        
        // Construct injective map tau
        map_t<vec_t<scalar_t>, Z> tau = constructTau<Z>(cutset);

        uint_t nRandPoints {20};
        for (uint_t i = 0; i < nRandPoints; ++i) {
            Z zeta = Z::rand();
            Z r = evaluateProofPolynomial(fg, cutset, etas, tau, zeta);
            Z rv = evaluateProofPolynomialVandermonde(
                fg, cutset, etas, tau, zeta);
            assert ( r == rv );
        }
    }
}

template <typename Z>
void test_verifyProofPolynomial(bool verbose=false,
                                uint_t nSamples=3)
{
    if (verbose)
        std::cout << "Testing verifyProofPolynomial<"
            << typeid(Z).name() << "> ...\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        // Construct a random FactorGraph and cutset
        FactorGraph<Z> fg = genRandFactorGraph<Z>(5, 4, 4);
        vec_t<Var> cutset = randVarSubset(fg, 3, 1);
        uint_t cutsetSize = cutset.size();
        uint_t boundarySize = (uint_t)rand() % (cutsetSize + 1);
        vec_t<Var> boundary = randSubset(cutset, boundarySize);
        fg.setBoundary(vectorToSet(boundary));
        
        // Construct injective map eta_i for i in cutset
        map_t<Var, map_t<uint_t,Z>> etas = constructEtas<Z>(cutset);
        
        // Construct injective map tau
        map_t<vec_t<scalar_t>, Z> tau = constructTau<Z>(cutset);
        
        // compute evaluations needed for inferring
        // the proof polynomial's coefficients
        long_uint_t d = proofPolynomialDegreeBound(fg, cutset);
        long_uint_t e = d+1;
        vec_t<Z> evalPoints(e);
        for (long_uint_t i = 0; i < e; ++i)
            evalPoints[i] = Z(i);
        vec_t<Z> evals(e);
        for (long_uint_t i = 0; i < e; ++i) {
            evals[i] = evaluateProofPolynomial(
                fg, cutset, etas, tau, evalPoints[i]);
        }
        
        // recover proof polynomial's coefficients
        Poly<Z> pp;
        assert (d < 1UL<<62); index_t sd = (index_t)d, se = (index_t)e;
        assert ( rsDecodeXY(sd, se, evalPoints.data(), evals.data(), pp) );
        
        // pp should now be a valid proof
        assert ( verifyProofPolynomial(pp, fg, cutset, etas, tau) );
        
        // after corrupting pp, verification should fail
        // (after multiple checks, provided d is much smaller than |Z|)
        uint_t corruptInd = (uint_t)rand() % pp.degree();
        Z corruptVal = Z::rand();
        while (corruptVal == pp[corruptInd])
            corruptVal = Z::rand();
        pp[corruptInd] = corruptVal;
        assert ( ! verifyProofPolynomial(pp, fg, cutset, etas, tau, 10) );
    }
}

template <typename Z>
void test_boundaryMapFromProof(bool verbose=false,
                                uint_t nSamples=3)
{
    if (verbose)
        std::cout << "Testing boundaryMapFromProof<"
            << typeid(Z).name() << "> ...\n";
    
    for (scalar_t s = 0; s < nSamples; ++s) {
        // Construct a random FactorGraph and cutset
        FactorGraph<Z> fg = genRandFactorGraph<Z>(4, 3, 3);
        vec_t<Var> cutset = randVarSubset(fg, 3, 0);
        uint_t boundarySize = (uint_t)rand() % (cutset.size() + 1);
        vec_t<Var> boundary = randSubset(cutset, boundarySize);
        fg.setBoundary(vectorToSet(boundary));
        
        // Construct injective map eta_i for i in cutset
        map_t<Var, map_t<uint_t,Z>> etas = constructEtas<Z>(cutset);
        
        // Construct injective map tau
        map_t<vec_t<scalar_t>, Z> tau = constructTau<Z>(cutset);
        
        // Compute evaluations needed for inferring
        // the proof polynomial's coefficients
        long_uint_t d = proofPolynomialDegreeBound(fg, cutset);
        long_uint_t e = d+1;
        vec_t<Z> evalPoints(e);
        for (long_uint_t i = 0; i < e; ++i)
            evalPoints[i] = Z(i);
        vec_t<Z> evals(e);
        for (long_uint_t i = 0; i < e; ++i) {
            evals[i] = evaluateProofPolynomial(
                fg, cutset, etas, tau, evalPoints[i]);
        }
        
        // Recover proof polynomial's coefficients, verify them
        Poly<Z> pp;
        assert (d < 1UL<<62); index_t sd = (index_t)d, se = (index_t)e;
        assert ( rsDecodeXY(sd, se, evalPoints.data(), evals.data(), pp) );
        assert ( verifyProofPolynomial(pp, fg, cutset, etas, tau) );
        
        Factor<Z> marginal = boundaryMapFromProof(
            pp, cutset, boundary, tau);
        FactorGraph<Z> fgContracted = fg.contractAll();
        vec_t<Factor<Z>> factorsContracted = fgContracted.getFactors();
        assert ( factorsContracted.size() == 1 );
        Factor<Z> expectedMarginal = factorsContracted[0];
        
        assertEqualUpTillPerm(marginal, expectedMarginal);
    }
}

template <typename Z>
void test_contractViaProof(bool verbose=false,
                           uint_t nSamples=5)
{
    if (verbose)
        std::cout << "Testing contractViaProof<"
            << typeid(Z).name() << "> ...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        FactorGraph<Z> fg = genRandFactorGraph<Z>(4, 3, 3);
        vec_t<Var> cutset = genRandCutset(fg);
        assertEqualUpTillPerm(
            contractViaProof<Z>(fg, cutset),
            fg.contractAll().getFactors()[0]
        );
    }
}


void test_MPInt(bool verbose=false, uint_t nSamples=100)
{
    if (verbose) std::cout << "Testing class MPInt...\n";
    
    // test constructors and operator=
    for (uint_t s = 0; s < nSamples; ++s) {
        uint_t r = (uint)rand();
        mpz_t m; mpz_init(m); mpz_set_ui(m, r);
        
        MPInt x(r);
        mpz_t mpzX; x.copyToMPZ(mpzX);
        assert ( mpz_cmp(mpzX, m) == 0 );
        
        MPInt y(m);
        mpz_t mpzY; y.copyToMPZ(mpzY);
        assert ( mpz_cmp(mpzY, m) == 0 );
        
        MPInt z(x);
        mpz_t mpzZ; z.copyToMPZ(mpzZ);
        assert ( mpz_cmp(mpzX, mpzZ) == 0 );
        
        MPInt w;
        w = y;
        mpz_t mpzW; w.copyToMPZ(mpzW);
        assert ( mpz_cmp(mpzW, mpzY) == 0 );
        
        mpz_clears(m, mpzX, mpzY, mpzZ, mpzW, NULL);
    }
    
    // test rand
    int lower {-666666}, upper {999999};
    MPInt sum(0L);
    for (uint_t s = 0; s < nSamples; ++s) {
        MPInt x = MPInt::randInRange(lower, upper);
        assert ( x.toLongInt() > lower );
        assert ( x.toLongInt() < upper );
        sum += x;
    }
    MPInt expectation = (upper - lower) / 2;
    assert ( sum / nSamples > 0.7 * expectation );
    assert ( sum / nSamples < 1.3 * expectation );
    
    // test arithmetic
    for (uint_t s = 0; s < nSamples; ++s) {
        uint_t r1 = (uint)rand();
        uint_t r2 = (uint)rand();
        mpz_t m1; mpz_init(m1); mpz_set_ui(m1, r1);
        mpz_t m2; mpz_init(m2); mpz_set_ui(m2, r2);
        MPInt x1(r1), x2(r2);
        mpz_t result; mpz_init(result);
        
        mpz_add(result, m1, m2);
        MPInt sum = x1 + x2;
        mpz_t sumCopy; sum.copyToMPZ(sumCopy);
        assert ( mpz_cmp(sumCopy, result) == 0 );
        
        mpz_sub(result, m1, m2);
        MPInt diff = x1 - x2;
        mpz_t diffCopy; diff.copyToMPZ(diffCopy);
        assert ( mpz_cmp(diffCopy, result) == 0 );
        
        mpz_mul(result, m1, m2);
        MPInt prod = x1 * x2;
        mpz_t prodCopy; prod.copyToMPZ(prodCopy);
        assert ( mpz_cmp(prodCopy, result) == 0 );
        
        MPInt x3 = x1;
        x3 += x2;
        mpz_t x3Copy; x3.copyToMPZ(x3Copy);
        assert ( mpz_cmp(x3Copy, sumCopy) == 0 );
        
        MPInt x4 = x1;
        x4 *= x2;
        mpz_t x4Copy; x4.copyToMPZ(x4Copy);
        assert ( mpz_cmp(x4Copy, prodCopy) == 0 );
        
        MPInt negx1 = -x1;
        MPInt x0 = x1 + negx1;
        mpz_t x0Copy; x0.copyToMPZ(x0Copy);
        mpz_set_ui(result, 0);
        assert ( mpz_cmp(x0Copy, result) == 0 );
        
        MPInt x5 = MPInt::randInRange(-999999, 999999);
        MPInt x6 = MPInt::randInRange(-100000, 100000);
        while (x6.toInt() == 0) x6 = MPInt::randInRange(-100000, 100000);
        assert ( (x5 / x6).toLongInt() == x5.toLongInt() / x6.toLongInt() );
        assert ( (x5 % x6).toLongInt() == x5.toLongInt() % x6.toLongInt() );
        
        uint_t modulus = (uint_t)rand();
        MPInt x7(modulus);
        long int rem = x5.toLongInt() % modulus;
        assert ( mod(x5, x7).toLongInt() == (rem < 0 ? modulus + rem : rem) );
        
        // modular inverses
        MPInt mpi67(67);
        MPInt x8 = MPInt::rand();
        while (x8.toInt() % 67 == 0) x8 = MPInt::rand();
        MPInt x8InvMod67 = mpi67.modInverse(x8);
        assert ( x8InvMod67.toInt() > 0 );
        assert ( (x8InvMod67.toInt() * (x8.toInt() % 67)) % 67 == 1 );
        
        mpz_clears(m1, m2, result, sumCopy, diffCopy, prodCopy, NULL);
        mpz_clears(x3Copy, x4Copy, x0Copy, NULL);
    }
}

void test_floatParsing(bool verbose=false)
{
    if (verbose) std::cout << "Testing functions for parsing floats...\n";
    
    float f {0};
    assert ( getFloatSign(f) == 0 );
    assert ( getFloatExpField(f) == 0 );
    assert ( getFloatSignificand(f) == 0 );
    assert ( getLSBExp(f) == 0 );
    
    f = 1;
    assert ( getFloatSign(f) == 0 );
    assert ( getFloatExpField(f) == 127 );
    assert ( getFloatExp(f) == 0 );
    assert ( getFloatSignificand(f) == (1<<23) );
    assert ( getLSBExp(f) == 0 );
    
    f = -1097.7415771484375;
    assert ( getFloatSign(f) == 1 );
    assert ( getFloatExpField(f) == 137 );
    assert ( getFloatExp(f) == 10 );
    assert ( getFloatSignificand(f) == 8992699 );
    assert ( getLSBExp(f) == -13 );
    
    f = 13492640368413930160128.0;
    assert ( getFloatSign(f) == 0 );
    assert ( getFloatExpField(f) == 200 );
    assert ( getFloatExp(f) == 73 );
    assert ( getFloatSignificand(f) == 11983872 );
    assert ( getLSBExp(f) == 60 );
    
}

void test_floatToMPInt(bool verbose=false)
{
    if (verbose) std::cout << "Testing floatToMPInt...\n";
    
    float f = -1097.7415771484375;
    MPInt x = floatToMPInt(f, 13);
    assert ( x.toInt() == -8992699 );
    x = floatToMPInt(f, 15);
    assert ( x.toLongInt() == -8992699 * 4 );
    
    f = 13492640368413930160128.0;
    /* The significand of f is 0110110 11011100 00000000;
     * the interpreted exponent of f is 73;
     * the i'th bit of the significand is interpreted as 2^(-i);
     * thus the last bit of f corresponds to 2^(73 - 23) = 2^(50),
     * and thus the last non-zero bit of f corresponds to 2^(60).
     * Thus: f is already an integer, and (f divided by 2^(60))
     * should still be an integer, equal to f / 2^(10). */
    x = floatToMPInt(f, -60);
    assert ( x.toLongInt() * 1024 == 11983872 );
    
    f = -1;
    x = floatToMPInt(f, 6);
    assert ( x.toInt() == -64 );
}

void test_intConversionExp(bool verbose=false)
{
    if (verbose) std::cout << "Testing intConversionExp...\n";
    
    Var v1((label_t)rand(), 2), v2((label_t)rand(), 3);
    vec_t<float> data {0.5, 9, -0.25, 0.25, 0.0625, -9};
    Factor<float> f((label_t)rand(), {v1,v2}, data.data());
    assert ( intConversionExp(f, true) == 4 );
    
    data = {4, 16, -32, 8, 0.0, 32};
    f = Factor<float>((label_t)rand(), {v1,v2}, data.data());
    assert ( intConversionExp(f, true, -150) == 0 );
    assert ( intConversionExp(f, false, -150) == -2 );
    
    data[1] = 15;  // in binary, 15 has a '1' at position 0
    f = Factor<float>((label_t)rand(), {v1,v2}, data.data());
    assert ( intConversionExp(f, true, -150) == 0 );
    assert ( intConversionExp(f, false, -150) == 0 );
    
    data = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    f = Factor<float>((label_t)rand(), {v1,v2}, data.data());
    /* if f is all zeros, intConversionExp should return ifZeroReturn
     * regardless of whether nonnegative is true or false. */
    assert ( intConversionExp(f, true, -150) == -150 );
    assert ( intConversionExp(f, false, -150) == -150 );
}

void test_mpIntToFloat(bool verbose=false, uint_t nSamples=1000)
{
    if (verbose) std::cout << "Testing mpIntToFloat...\n";
    
    vec_t<float> edgeCases {
        0.0,
        std::numeric_limits<float>::min(),
        -std::numeric_limits<float>::min(),
        std::numeric_limits<float>::lowest(),
        -std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::denorm_min(),
        -std::numeric_limits<float>::denorm_min(),
        std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max()
    };
    for (float f : edgeCases) {
        int ce = -getLSBExp(f);
        MPInt i = floatToMPInt(f, ce);
        float g = mpIntToFloat(i, -ce);
        assert ( g == f );
    }
    for (uint_t s = 0; s < nSamples; ++s) {
        float f = genRandScalar<float>();
        int ce = -getLSBExp(f);
        MPInt i = floatToMPInt(f, ce);
        float g = mpIntToFloat(i, -ce);
        assert ( g == f );
    }
}

void test_floatToMPIntFG(bool verbose=false, uint_t nSamples=20)
{
    if (verbose) std::cout << "Testing floatToMPIntFG...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        
        FactorGraph<float> gF = genRandFactorGraph<float>(
            5, 4, 4, {}, false, 1, 2, &genRandScalar<float>);
        
        map_t<label_t, int> conversionExps;
        map_t<label_t, int> upperBoundExps;
        FactorGraph<MPInt> gI = floatToMPIntFG(
            gF, conversionExps, upperBoundExps, false);
        
        assert ( gI.nFactors() == gF.nFactors() );
        for (Factor<float> f : gF.getFactors()) {
            assert ( isIn(conversionExps, f.label) );
            assert ( intConversionExp(f, false) == conversionExps[f.label] );
            assert ( isIn(gI.getFactorIds(), f.label) );
            Factor<MPInt> fI = gI.getFactorWithId(f.label);
            assert ( fI.volume == f.volume );
            double upperBound = std::pow(2, upperBoundExps[f.label]);
            for (uint_t k = 0; k < f.volume; ++k) {
                assert (
                    mpIntToFloat(fI.valueAt(k), -conversionExps[fI.label])
                    == f.valueAt(k)
                );
                assert ( f.valueAt(k) <= upperBound );
            }
        }
    }
}

template <typename P>
void test_mpIntToMod(bool verbose=false, uint_t nSamples=200)
{
    if (verbose)
        std::cout << "Testing mpIntToMod<"
            << typeid(P).name() << "> ...\n";
    
    vec_t<uint_t> edgeCases {
        std::numeric_limits<uint_t>::min(),
        std::numeric_limits<uint_t>::max()
    };
    for (uint_t x : edgeCases) {
        MPInt mpIntX(x);
        Zp<P> zpX = mpIntToMod<P>(mpIntX);
        assert ( zpX.value() == x % P::modulus );
    }
    for (uint_t s = 0; s < nSamples; ++s) {
        MPInt x = MPInt::rand() * MPInt::rand();
        Zp<P> zpX = mpIntToMod<P>(x);
        uint_t correctValue = (uint_t)(x % P::modulus).toLongInt();
        assert ( zpX.value() == correctValue );
    }
}

template <typename P>
void test_mpIntToModFactor(bool verbose=false, uint_t nSamples=10)
{
    /* NOTE: If MPInt::random() is modified to produce values larger than
    * uint_t, then this test will fail, and need to be replaced/rewritten. */
    if (verbose)
        std::cout << "Testing mpIntToModFactor<"
            << typeid(P).name() << "> ...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        Factor<MPInt> fMPInt = genRandFactor<MPInt>(
            5, 4, {1,2,3,4,5}, {}, 2, &(MPInt::rand)).first;
        Factor<Zp<P>> fMod = mpIntToModFactor<P>(fMPInt);
        assert ( fMod.volume == fMPInt.volume );
        assert ( fMod.nVars() == fMPInt.nVars() );
        for (uint_t k = 0; k < fMod.volume; ++k) {
            assert (fMod.valueAt(k).value()
                    == (uint_t)fMPInt.valueAt(k).toInt() % P::modulus);
        }
    }
}

void test_floatToMPIntToFloatFactor(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing mpIntToFloatFactor and floatToMPIntFactor...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        // float -> MPInt -> float
        Factor<float> floatFactor1 = genRandFactor<float>(4,3,{1,2,3,4,5},2);
        int cExp = intConversionExp(floatFactor1);
        Factor<MPInt> mpIntFactor = floatToMPIntFactor(floatFactor1, cExp);
        Factor<float> floatFactor2 = mpIntToFloatFactor(mpIntFactor, -cExp);
        assertFactorsEqual(floatFactor1, floatFactor2);
    }
}

void test_crtReconstruct(bool verbose=false, uint_t nSamples=50)
{
    if (verbose)
        std::cout << "Testing crtReconstruct on random volume-1 Factors...\n";
    
    vec_t<uint_t> moduli {5,7,11,13,17,19,23,29,31,37,41,43,47,53,59};
    vec_t<MPInt> primes {5,7,11,13,17,19,23,29,31,37,41,43,47,53,59};
    MPInt N = vectorProduct(primes);
    MPInt M = N/2 - 1;
    
    for (uint_t s = 0; s < nSamples; ++s) {
        MPInt x = MPInt::randInRange(-M, M);
        label_t fId = (label_t)rand();
        Factor<MPInt> f(fId, {}, &x);
        vec_t<Factor<uint_t>> gs;
        for (MPInt p : primes) {
            uint_t gVal = (uint_t)mod(x, p).toInt();
            gs.push_back(Factor<uint_t>((label_t)rand(), {}, &gVal));
        }
        Factor<MPInt> fReconstructed = crtReconstruct(gs, moduli, fId, true);
        assert ( fReconstructed.isScalar() );
        assert ( fReconstructed.scalarValue() == x );
    }
    
    if (verbose)
        std::cout << "Testing crtReconstruct on random Factors...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        Factor<MPInt> fMPInt = genRandFactor<MPInt>(3, 3, {1,2,3,4,5}, 2);
        uint_t vol = fMPInt.volume;
        for (uint_t k = 0; k < vol; ++k)
            *fMPInt.dataAt(k) = MPInt::randInRange(-M, M);
        
        vec_t<Factor<uint_t>> gs;
        for (MPInt p : primes) {
            vec_t<uint_t> data(vol);
            for (uint_t k = 0; k < vol; ++k)
                data[k] = (uint_t)(mod(fMPInt.valueAt(k), p).toInt());
            gs.push_back(
                Factor<uint_t>((label_t)rand(), fMPInt.vars, data.data())
            );
        }
        Factor<MPInt> fReconstructed = crtReconstruct(
            gs, moduli, fMPInt.label, true);
        assertFactorsEqual(fReconstructed, fMPInt);
    }
}

void test_primesToReach(bool verbose=false, uint_t nSamples=50)
{
    if (verbose)
        std::cout << "Testing primesToReach...\n";
    
    vec_t<uint_t> allPrimes{5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
    for (uint_t s = 0; s < nSamples; ++s) {
        MPInt x = MPInt::randInRange(0L, 1L<<63);
        vec_t<uint_t> primes = primesToReach(x, allPrimes);
        uint_t nPrimes = primes.size();
        MPInt product{1};
        for (uint_t i = 0; i < nPrimes; ++i)
            product *= MPInt(primes[i]);
        assert ( product >= x );
        uint_t lastPrime = allPrimes[nPrimes - 1];
        assert ( product / lastPrime  <= x  );
    }
}

template <typename Z>
void test_fileio_rw_misc(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing writing/loading Factor<"
            << typeid(Z).name() << "> and numeric params to/from file...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        char filename[] = "temp_fileio_test_file_XXXXXX";
        FILE* fptr = fdopen(mkstemp(filename), "r+");
        assert ( fptr != nullptr );
        FileHandler* fh1 = new FileHandler(fptr);
        
        Factor<Z> f1 = genRandFactor<Z>(4, 4, {101,102,103,104}, 2);
        if (verbose) {
            std::cout << "First Factor to be written to file:\n";
            printFactor(f1);
        }
        writeFactorToFile<Z>(*fh1, f1);
        
        uint_t param1 = (uint_t)rand();
        str_t param1Name = "modulus";
        writeParamToFile<uint_t>(*fh1, param1Name, param1);
        
        Factor<Z> f2 = genRandFactor<Z>(4, 4, {101,102,103,104}, 2);
        if (verbose) {
            std::cout << "Second Factor to be written to file:\n";
            printFactor(f2);
        }
        writeFactorToFile<Z>(*fh1, f2);
        
        long param2 = rand()<<1;
        str_t param2Name = "conversionExponent";
        writeParamToFile<long>(*fh1, param2Name, param2);
        
        delete fh1;
        
        FileHandler fh2 = FileHandler(str_t(filename));
        
        Factor<Z> loadedF1;
        loadObject<Z>(fh2, loadedF1);
        if (verbose) {
            std::cout << "First Factor loaded from file:\n";
            printFactor(loadedF1);
        }
        assertFactorsEqual(f1, loadedF1);
        
        // Move just past the first Factor
        fh2.toBeginning(); fh2.toNextOccurrence("Factor"); fh2.nextToken();
        Factor<Z> loadedF2;
        loadObject<Z>(fh2, loadedF2);
        if (verbose) {
            std::cout << "Second Factor loaded from file:\n";
            printFactor(loadedF2);
        }
        assertFactorsEqual(f2, loadedF2);
        
        uint_t loadedParam1;
        assert ( loadParamFromFile(fh2, param1Name, loadedParam1) );
        if (verbose) {
            std::cout << "Expecting first loaded param "<< loadedParam1
                << " to be equal to first written param " << param1 << "\n";
        }
        assert ( loadedParam1 == param1 );
        
        long loadedParam2;
        assert ( loadParamFromFile<long>(fh2, param2Name, loadedParam2) );
        if (verbose) {
            std::cout << "Expecting second loaded param "<< loadedParam2
                << " to be equal to second written param " << param2 << "\n";
        }
        assert ( loadedParam2 == param2 );
        
        fh2.close();
        remove(filename);
    }
}

template <typename Z>
void test_fileio_factorgraph_rw(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing writing/loading FactorGraph<"
            << typeid(Z).name() << "> to/from file...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        FactorGraph<Z> g = genRandFactorGraph<Z>(4, 4, 3);
        if (verbose) {
            std::cout << "FactorGraph to be written to file:\n";
            printFactorGraph(g);
        }
        char filename[] = "temp_fileio_test_file_XXXXXX";
        FILE* fptr = fdopen(mkstemp(filename), "r+");
        assert ( fptr != nullptr );
        FileHandler* fh1 = new FileHandler(fptr);
        writeFGToFile<Z>(*fh1, g);
        delete fh1;
        
        FileHandler fh2 = FileHandler(str_t(filename));
        FactorGraph<Z> loadedG;
        loadFactorGraph<Z>(fh2, loadedG);
        if (verbose) {
            std::cout << "FactorGraph loaded from file:\n";
            printFactorGraph(loadedG);
        }
        
        vec_t<Var> gVars = g.getVars();
        vec_t<Var> loadedVars = loadedG.getVars();
        assert ( isSuperset(gVars, loadedVars) );
        assert ( isSuperset(loadedVars, gVars) );
        for (Factor<Z> f : g.getFactors()) {
            Factor<Z> loadedF = loadedG.getFactorWithId(f.label);
            assert ( f.getVolume() == loadedF.getVolume() );
            for (long_uint_t k = 0; k < f.getVolume(); ++k)
                assert ( f.valueAt(k) == loadedF.valueAt(k) );
        }
        
        fh2.close();
        remove(filename);
    }
}

template <typename Z>
void test_fileio_map_rw(bool verbose=false, uint_t nSamples=10)
{
    if (verbose)
        std::cout << "Testing writing/loading map_t to/from file...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        char filename[] = "temp_fileio_test_file_XXXXXX";
        FILE* fptr = fdopen(mkstemp(filename), "r+");
        assert ( fptr != nullptr );
        FileHandler* fh1 = new FileHandler(fptr);
        
        uint_t n = (uint_t)rand() % 100;
        map_t<int, Z> randmap;
        for (uint_t i = 0; i < n; ++i)
            randmap.insert({rand(), Z((uint_t)rand())});
        writeObject(*fh1, randmap);
        delete fh1;
        
        FileHandler fh2 = FileHandler(str_t(filename));
        map_t<int, Z> loadedMap;
        loadObject(fh2, loadedMap);
        
        for (auto it = randmap.begin(); it != randmap.end(); ++it)
            assert ( loadedMap.at(it->first) == it->second );
        
        fh2.close();
        remove(filename);
    }
}



// NOTE: This test is only intended for types float, MPInt, or Zp
template <typename T>
void test_integration_generate(bool verbose=false, uint_t modulus=0,
                               uint_t nSamples=5)
{
    if (verbose)
        std::cout << "Testing integration / generate, with type "
                  << typeid(T).name() << "...\n";

    str_t type;
    if constexpr (std::is_same_v<T,float>)
        type = "float";
    else if constexpr (std::is_same_v<T,MPInt>)
        type = "MPInt";
    else {
        type = "modular";
        assert ( isIn(IDIGM_PRIMES, modulus) );
    }

    for (uint_t s = 0; s < nSamples; ++s) {
        
        // Define parameters for generating a random factor graph
        uint_t maxNFactors = 2 + (uint_t)rand() % 5;
        uint_t maxNVarsPerFactor = 2 + (uint_t)rand() % 5;
        uint_t maxVarSize = 4;
        vec_t<label_t> varIds;
        for (uint_t i = 0; i < maxNVarsPerFactor * 3; ++i)
            varIds.push_back((label_t)rand());
        uint_t minNFactors = 2;
        uint_t minVarSize = 3;
        
        // Write the parameters to a file
        char paramsFile[] = "temp_test_file_params-for-gen_XXXXXX";
        FILE* fptr = fdopen(mkstemp(paramsFile), "r+");
        assert ( fptr != nullptr );
        FileHandler hParams(fptr);
        writeParamToFile(hParams, "type", type);
        if (type == "modular")
            writeParamToFile(hParams, "modulus", modulus);
        writeParamToFile(hParams, "maxNFactors", maxNFactors);
        writeParamToFile(hParams, "maxNVarsPerFactor", maxNVarsPerFactor);
        writeParamToFile(hParams, "maxVarSize", maxVarSize);
        writeParamToFile(hParams, "varIds", varIds);
        writeParamToFile(hParams, "minVarSize", minVarSize);
        writeParamToFile(hParams, "minNFactors", minNFactors);
        writeParamToFile(hParams, "randomFactorLabels", true);
        hParams.close();
        
        // Generate a random factor graph
        str_t fgFile = generate(str_t(paramsFile));
        assert ( fgFile != "" );
        
        // Load the generated factor graph;
        // check that it is consistent with the parameters
        FileHandler hFG(fgFile);
        FactorGraph<T> fg;
        assert ( loadFactorGraph<T>(hFG, fg) );
        set_t<Var> cutset;
        assert ( loadParamFromFile(hFG, "cutset", cutset) );
        hFG.close();
        assert ( isSuperset(cutset, fg.getBoundary()) );
        assert ( isSuperset(vectorToSet(fg.getVars()), cutset) );
        assert ( fg.nFactors() <= maxNFactors );
        assert ( fg.nFactors() >= minNFactors );
        for (Factor<T> f : fg.getFactors()) {
            assert ( f.nVars() <= maxNVarsPerFactor ) ;
            for (Var v : f.getVars()) {
                assert ( v.size <= maxVarSize );
                assert ( v.size >= minVarSize );
                assert ( isIn(varIds, v.label) );
            }
        }
        
        remove(paramsFile);
        remove(fgFile.c_str());
    }
}

//NOTE: Only for types float and Zp
template <typename T>
void test_integration_contract(bool verbose=false, uint_t modulus=0,
                               uint_t nSamples=5)
{
    str_t type;
    if constexpr (std::is_same_v<T,float>)
        type = "float";
    else {
        type = "modular";
        assert ( isIn(IDIGM_PRIMES, modulus) );
    }
    
    if (verbose)
        std::cout << "Testing integration / contract, for type "
                  << type << " ...\n";
    
    for (uint_t s = 0; s < nSamples; ++s) {
        
        // Generate a random factor graph
        uint_t maxNFactors = 2 + (uint_t)rand() % 4;
        uint_t maxNVarsPerFactor = 2 + (uint_t)rand() % 4;
        uint_t maxVarSize = 3;
        vec_t<label_t> varIds;
        for (uint_t i = 0; i < maxNVarsPerFactor * 3; ++i)
            varIds.push_back((label_t)rand());
        uint_t minNFactors = 1;
        uint_t minVarSize = 2;
        uint_t maxBoundarySize = 3;
        FactorGraph<T> fg = genRandFactorGraph<T>(
            maxNFactors, maxNVarsPerFactor, maxVarSize, varIds,
            true, minNFactors, minVarSize, nullptr, 0, maxBoundarySize
        );
        
        // Write the FactorGraph to file
        char fgFile[] = "temp_test_file_factorgraph_XXXXXX";
        FILE* fptr = fdopen(mkstemp(fgFile), "r+");
        assert ( fptr != nullptr );
        FileHandler hFG(fptr);
        writeFGToFile(hFG, fg);
        writeParamToFile(hFG, "type", type);
        if (type == "modular")
            writeParamToFile(hFG, "modulus", modulus);
        writeParamToFile(hFG, "cutset", genRandCutset(fg, maxBoundarySize+2));
        hFG.close();
        
        // Contract using integration::contract
        Options opts(false, {});
        str_t resultFile = contract(fgFile, opts);
        assert ( resultFile != "" );
        FileHandler hResult(resultFile);
        Factor<T> result;
        if (type == "float")
            assert ( loadObject(hResult, result) );
        else {
            Factor<uint_t> uintResult;
            assert ( loadObject(hResult, uintResult) );
            result = uintToModFactor<T>(uintResult);
        }
        
        // Compare to result of directly contracting fg
        vec_t<Factor<T>> factors = fg.contractAll().getFactors();
        assert ( factors.size() == 1 );
        Factor<T> expectedResult = factors[0];
        if (type == "float")
            assertEqualUpTillPerm(result, expectedResult, false, 0.0001);
        else
            assertEqualUpTillPerm(result, expectedResult, false, 0);
        remove(fgFile);
        remove(resultFile.c_str());
    }
}


void test_tb_solvePermanentBase_int(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing testbench/solvePermanentBase<int>, "
                  << "with fixed examples...\n";
    
    vec_t<int> elems {1,2,3,4,5,6,7,8,9};
    Array2D matrix(3, 3, elems.data());
    int permanent = solvePermanentBase(matrix);
    assert ( permanent == 450 );
    
    elems = {0,5,9,5,3,0,2,2,2};
    matrix = Array2D(3, 3, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == 194 );
    
    elems = {6,-5,-2,8, -5,-4,6,-6, -8,-1,-8,5, 0,10,-1,5};
    matrix = Array2D(4, 4, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == 4093 );
    
    elems = {0,-2, 5,-7};
    matrix = Array2D(2, 2, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == -10 );
    
    elems = {rand()};
    matrix = Array2D(1, 1, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == elems[0] );
}

template <typename Z>
void test_tb_solvePermanent(bool verbose=false)
{
    if (verbose)
        std::cout << "Testing solvePermanentDFGC and solvePermanentBase "
                  << "with type <" << typeid(Z).name()
                  << "> with fixed examples...\n";
    
    vec_t<Z> elems {1,2,3,4,5,6,7,8,9};
    Array2D matrix(3, 3, elems.data());
    Z permanent = solvePermanentBase(matrix);
    assert ( permanent == Z(450) );
    permanent = solvePermanentDFGC(matrix);
    assert ( permanent == Z(450) );
    
    elems = {0,5,9, 5,3,0, 2,2,2};
    matrix = Array2D(3, 3, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == Z(194) );
    permanent = solvePermanentDFGC(matrix);
    assert ( permanent == Z(194) );
    
    elems = {7,2,9,10, 7,3,6,7, 4,8,6,7, 5,10,2,6};
    matrix = Array2D(4, 4, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == Z(36997) );
    permanent = solvePermanentDFGC(matrix);
    assert ( permanent == Z(36997) );
    
    elems = {9,4, 5,2};
    matrix = Array2D(2, 2, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == Z(38) );
    permanent = solvePermanentDFGC(matrix);
    assert ( permanent == Z(38) );
    
    elems = {Z((uint_t)rand())};
    matrix = Array2D(1, 1, elems.data());
    permanent = solvePermanentBase(matrix);
    assert ( permanent == elems[0] );
    permanent = solvePermanentDFGC(matrix);
    assert ( permanent == elems[0] );
}



int main(int argc, char* argv[])
{
    bool verbose {false};
    time_t seed = time(NULL);
    
    for (int i = 0; i < argc; ++i) {
        if (!strcmp(argv[i], "verbose"))
            verbose = true;
        if (!strcmp(argv[i], "seed")) {
            assert ( argc > i+1 );
            seed = atoi(argv[i+1]);
        }
    }
    srand(seed);
    if (verbose)
        std::cout << "Using random seed " << seed << "\n";

    // montgomery.cpp
    test_Montgomery_arithmetic<MP>(verbose);
    test_Zp_interface<MP>(verbose);
    test_Montgomery_arithmetic<MP67>(verbose, 1000);
    test_Zp_interface<MP67>(verbose, 1000);
    
    // utils.cpp
    test_arrayProduct();
    test_getCoords();
    test_arraySliceAtDim_fixed_example();
    test_arraySlice_fixed_examples();
    test_vectorProduct(verbose);
    if (verbose) test_printNDArray_fixed_examples();
    test_indexOf_fixed_examples(verbose);
    test_isSuperset_fixed_examples(verbose);
    test_moveToEnd_fixed_examples(verbose);
    test_getPermutation_fixed_examples(verbose);
    test_naiveMatrixMul_examples<int>(verbose);
    test_naiveMatrixMul_examples<Zp<MP67>>(verbose);
    test_randSubset<int>(verbose);
    test_genRandFloatInRange(verbose);
    test_permuteArray_fixed_examples();
    test_perm_generation_and_inversion(verbose,10);
    test_permuteArray_perm_and_inverse<Zp<MP>>(verbose);
    test_permuteVector(verbose);
    
    // factor.cpp
    test_sumOverVar_fixed_examples(verbose);
    test_contractFactors_fixed_examples<Zp<MP>>(verbose);
    if (verbose) print_contractFactors_random<Zp<MP67>>();
    test_contractFactors_edge_cases<int>(verbose);
    test_Factor_slice<Zp<MP67>>(verbose);
    
    // factorgraph.cpp
    test_FactorGraph_add_remove<int>(verbose);
    test_FactorGraph_getInternalVarIds_examples<int>(verbose);
    test_FactorGraph_contractInOrder<int>(verbose);
    test_FactorGraph_contractAll<int>(verbose);
    test_FactorGraph_slice<Zp<MP67>>(verbose);
    test_contractionCost<Zp<MP67>>(verbose);
    test_optimalContractionOrder<Zp<MP>>(verbose);

    // class Poly
    test_Poly_rand<Zp<MP67>>(verbose);
    test_Poly_eval<Zp<MP67>>(verbose);
    test_Poly_arithmetic_plus<Zp<MP67>>(verbose);
    test_Poly_arithmetic_mul<Zp<MP67>>(verbose);
    test_Poly_arithmetic_quorem<Zp<MP67>>(verbose);
    test_Poly_interpolate<Zp<MP67>>(verbose);
    test_Poly_batchEval<Zp<MP67>>(verbose);
    test_rsDecodeXY<Zp<MP67>>(verbose);
    
    // Vandermonde
    test_invertVandermonde<Zp<MP67>>(verbose);
    
    // eval_algo
    test_enumerateVectors(verbose);
    test_nextCombination(verbose);
    test_transformVarVandermonde<Zp<MP67>>(verbose);
    test_evaluateProofPolynomialVandermonde<Zp<MP>>(verbose);
    test_evaluateProofPolynomial<Zp<MP67>>(verbose);
    test_evaluateProofPolynomial_compare<Zp<MP>>(verbose);
    test_verifyProofPolynomial<Zp<MP>>(verbose);
    test_boundaryMapFromProof<Zp<MP>>(verbose);
    test_contractViaProof<Zp<MP>>(verbose);

    // Chinese remaindering, conversions to/from MPInt
    test_MPInt(verbose);
    test_floatParsing(verbose);
    test_floatToMPInt(verbose);
    test_intConversionExp(verbose);
    test_mpIntToFloat(verbose);
    test_floatToMPIntFG(verbose);
    test_mpIntToMod<MP67>(verbose);
    test_mpIntToModFactor<MP67>(verbose);
    test_floatToMPIntToFloatFactor(verbose);
    test_crtReconstruct(verbose);
    test_primesToReach(verbose);

    // fileio.cpp
    test_fileio_rw_misc<Zp<MP67>>(verbose);
    test_fileio_rw_misc<float>(verbose);
    test_fileio_rw_misc<MPInt>(verbose);
    test_fileio_factorgraph_rw<Zp<MP67>>(verbose);
    test_fileio_factorgraph_rw<float>(verbose);
    test_fileio_factorgraph_rw<MPInt>(verbose);
    test_fileio_map_rw<float>(verbose);
    test_fileio_map_rw<Zp<MP67>>(verbose);

    // integration.cpp
    test_integration_generate<float>(verbose);
    test_integration_generate<MPInt>(verbose);
    test_integration_generate<Zp<MP>>(verbose, MP::modulus);
    test_integration_contract<float>(verbose);
    test_integration_contract<Zp<MP>>(verbose, MP::modulus);

    // testbench.cpp
    test_tb_solvePermanentBase_int(verbose);
    test_tb_solvePermanent<Zp<MP>>(verbose);

    return 0;
}


