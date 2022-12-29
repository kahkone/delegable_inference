#ifndef __IDIGM_EVAL_ALGO_H__
#define __IDIGM_EVAL_ALGO_H__

#include <map>
#include <vector>
#include <algorithm>
#include "util.cpp"
#include "factorgraph.cpp"
#include "poly.cpp"


//NOTE: Constructs (\prod_i ranges[i]) vectors in memory at once
template <typename T>
void enumerateVectors(
    vec_t<scalar_t> ranges,
    vec_t<vec_t<T>>& result,
    vec_t<T>& partialElement,
    scalar_t level=0)
{
    if (ranges.size() == 0)
        return;
    assert (0 <= level && level < ranges.size());
    if (level == 0)
        result.clear();
    
    for (scalar_t i = 0; i < ranges[level]; ++i) {
        vec_t<T> element = vec_t<T>(partialElement);
        element.push_back(T(i));
        if (level == ranges.size()-1)
            result.push_back(element);
        else
            enumerateVectors<T>(ranges, result, element, level+1);
    }
}

template <typename T>
void enumerateVectors(vec_t<scalar_t> ranges, vec_t<vec_t<T>>& result)
{
    vec_t<T> partialElement = {};
    enumerateVectors<T>(ranges, result, partialElement, 0);
}

template <typename T>
bool nextCombination(vec_t<T>& combination, const vec_t<T>& ranges)
{
    assert ( ranges.size() == combination.size() );
    assert ( combination.size() < 1UL<<63 );
    index_t n = (index_t)combination.size();
    if (n == 0) 
        return false;
    
    #if DEBUG_ASSERTS
    for (uint_t i = 0; i < ranges.size(); ++i) {
        assert ( ranges[i] > 0 );
        assert ( combination[i] <= ranges[i] - 1 );
        assert ( combination[i] >= 0 );
    }
    #endif
    
    for (index_t pos = n-1; pos >= 0; --pos) {
        uint_t p = (uint_t)pos;
        if (combination[p] < ranges[p] - 1) {
            combination[p] += T(1);
            return true;
        }
        else
            combination[p] = T(0);
    }
    return false;
}

template <typename F>
vec_t<F> powersOf(F x, uint_t n)
{
    vec_t<F> powers(n);
    F xPow = F::one;
    for (scalar_t p = 0; p < n; ++p) {
        powers[p] = xPow;
        xPow *= x;
    }
    return powers;
}

template <typename F>
void transformVarVandermonde(
    FactorGraph<F>& fg, const Var& v,
    F evalPoint, const map_t<scalar_t,F>& eta)
{
    vec_t<label_t> originalFactorIds {fg.getFactorIds()};
    for (label_t fId : originalFactorIds) {
        
        if ( ! fg.adjacent(fId, v) )
            continue;
        
        Var varQ = fg.createVar(v.size);
        fg.replaceVarInFactor(fId, v, varQ);
        
        vec_t<F> vandermondePoints(v.size);
        for (scalar_t d = 0; d < v.size; ++d)
            vandermondePoints[d] = eta.at(d);
        F* invVandermonde = new F[v.size * v.size];
        genInvVandermonde(v.size, vandermondePoints.data(), invVandermonde);
        Var varP = fg.createVar(v.size);
        fg.createFactor({varP, varQ}, invVandermonde);
        delete[] invVandermonde;
        
        vec_t<F> evalPointPowers = powersOf(evalPoint, v.size);
        fg.createFactor({varP}, evalPointPowers.data());
    }
    bool removed = fg.removeVarIfDetached(v);
    assert (removed);
}

//TODO: Refactor to use nextCombination
template <typename F>
map_t<Var, Poly<F>> interpolateElls(
    vec_t<Var> cutset,
    map_t<Var, map_t<scalar_t,F>> etas,
    map_t<vec_t<scalar_t>, F> tau)
{
    #if TRACK_TIME
    logTimeStart("interpolateElls");
    #endif
    // Construct interpolation inputs
    vec_t<vec_t<scalar_t>> cutsetDomain;
    enumerateVectors(varsToSizes(cutset), cutsetDomain);
    vec_t<F> ellDomain;
    for (vec_t<scalar_t> cutsetElement : cutsetDomain)
        ellDomain.push_back(tau.at(cutsetElement));
    assert ( ellDomain.size() < 1UL<<62 );
    
    map_t<Var, Poly<F>> ells;
    scalar_t nC = cutset.size();
    assert (etas.size() == nC);
    for (scalar_t i = 0; i < nC; ++i) {
        // Construct interpolation targets
        map_t<scalar_t,F> etaMap = etas.at(cutset[i]);
        vec_t<F> ellTargets;
        for (vec_t<scalar_t> cutsetElement : cutsetDomain)
            ellTargets.push_back(etaMap.at(cutsetElement[i]));
        // Interpolate \hat{\ell}_i
        Poly<F> ell = Poly<F>::interpolate(
            (index_t)ellDomain.size(), ellDomain.data(), ellTargets.data());
        ells[cutset[i]] = ell;
    }
    
    #if TRACK_TIME
    logTimeEnd("interpolateElls");
    #endif
    return ells;
}

/* NOTE: Extremely inefficient */
template <typename F>
F evaluateProofPolynomialVandermonde(
    const FactorGraph<F>& fgIn,
    const vec_t<Var>& cutset,
    const map_t<Var, map_t<scalar_t,F>>& etas,
    const map_t<vec_t<scalar_t>, F>& tau,
    F zeta)
{
    assert ( isSuperset(vectorToSet(cutset), fgIn.getBoundary()) );
    FactorGraph<F> fg(fgIn);
    
    // Interpolate \hat{\ell}_i for i in cutset
    map_t<Var, Poly<F>> ells = interpolateElls(cutset, etas, tau);
    
    // Perform local graph transformations for each cutset var
    for (Var v : cutset) {
        F ellZeta = ells[v](zeta);
        transformVarVandermonde(fg, v, ellZeta, etas.at(v));
    }
    
    // Contract the transformed factor graph; obtain scalar result
    fg = fg.contractAll();
    vec_t<Factor<F>> factors = fg.getFactors();
    assert ( factors.size() == 1 );
    assert ( factors[0].isScalar() );

    return factors[0].scalarValue();
}

bool contracts(label_t factorId, IdTriple t)
{
    return (t.first == factorId || t.second == factorId);
}

bool contracts(label_t factorId, const vec_t<IdTriple>& contractions)
{
    for (IdTriple c : contractions) {
        if (contracts(factorId, c))
            return true;
    }
    return false;
}

bool contractsAll(const vec_t<label_t>& factorIds,
                  const vec_t<IdTriple>& contractions)
{
    for (label_t fId : factorIds) {
        if ( ! contracts(fId, contractions))
            return false;
    }
    return true;
}

/* NOTE: Inefficient */
template <typename F>
F evaluateProofPolynomial(
    const FactorGraph<F>& fgIn,
    const vec_t<Var>& cutset,
    const map_t<Var, map_t<uint_t,F>>& etas,
    const map_t<Var, Poly<F>>& ells,
    const F zeta,
    const vec_t<IdTriple>& contractionOrder={})
{
    #if TRACK_TIME >= 2
    logTimeStart("evaluateProofPolynomial");
    #endif
    assert ( isSuperset(vectorToSet(cutset), fgIn.getBoundary()) );
    assert ( contractionOrder.size() == 0
             || contractsAll(fgIn.getFactorIds(), contractionOrder) );

    FactorGraph<F> fg(fgIn);
    // As cutset variables are eliminated, fg may end up with singly-connected
    // cutset variables; we do not want to sum over those
    fg.setBoundary(vectorToSet(cutset));
    
    vec_t<label_t> originalFactorIds = fg.getFactorIds();
    for (label_t fId : originalFactorIds) {
        
        for (Var v : cutset) {
            Factor<F> f = fg.getFactorWithId(fId);
            if ( ! f.adjacentToVar(v))
                continue;
            
            // Permute f so that v is last
            #if TRACK_TIME >= 2
            str_t volStr = std::to_string(f.getVolume());
            logTimeStart("evaluateProofPolynomial:permuteToEnd", volStr);
            f = f.permuteToEnd({v});
            logTimeEnd("evaluateProofPolynomial:permuteToEnd", volStr);
            #else
            f = f.permuteToEnd({v});
            #endif
            
            #if TRACK_TIME >= 2
            logTimeStart("evaluateProofPolynomial:interpolateAndTransform");
            #endif
            
            vec_t<uint_t> otherVarSizes = f.varSizesWithout({v});
            // slice will go through all possible slices of other vars
            vec_t<uint_t> slice(otherVarSizes.size(), 0);
            
            // For each value-combination (slice) of the other vars,
            // interpolate a corresponding polynomial over var v
            map_t<uint_t,F> eta = etas.at(v);
            vec_t<F> interpolationPoints(v.size);
            for (uint_t j = 0; j < v.size; ++j)
                interpolationPoints[j] = eta.at(j);
            F* coeffData = new F[f.volume];
            do {
                vec_t<uint_t> coords(slice);
                coords.push_back(0);  // point to beginning of v dimension
                long_uint_t offset = linearizeCoords(
                    f.nVars(), f.varSizes().data(), coords.data());
                Poly<F> p = Poly<F>::interpolate(
                    v.size, interpolationPoints.data(), f.dataAt(offset));
                p.copyCoeffsTo(coeffData + offset);
            } while (nextCombination(slice, otherVarSizes));
            
            // Update f to be the result of the interpolations
            fg.replaceFactor(fId, Factor<F>(fId, f.vars, coeffData));
            delete[] coeffData;
            
            // Replace v by a Var corresponding to powers of ell(zeta)
            Var varP = fg.createVar(v.size);
            fg.replaceVarInFactor(fId, v, varP);
            
            // Construct Factor e consisting of powers of ell(zeta)
            F ellZeta = ells.at(v)(zeta);
            vec_t<F> ellZetaPowers = powersOf(ellZeta, v.size);
            label_t eId = fg.createFactor({varP}, ellZetaPowers.data());
            #if TRACK_TIME >= 2
            logTimeEnd("evaluateProofPolynomial:interpolateAndTransform");
            #endif
            
            // Replace f with the result of contracting e with f
            #if TRACK_TIME >= 2
            logTimeStart("evaluateProofPolynomial:contractE");
            fg = fg.contractInOrder({{fId, eId, fId}});
            logTimeEnd("evaluateProofPolynomial:contractE");
            #else
            fg = fg.contractInOrder({{fId, eId, fId}});
            #endif
        }
    }
    
    #if TRACK_TIME >= 2
    logTimeStart("evaluateProofPolynomial:finalContraction");
    #endif
    if (contractionOrder.size() == 0)
        fg = fg.contractAll();
    else
        fg = fg.contractInOrder(contractionOrder);
    #if TRACK_TIME >= 2
    logTimeEnd("evaluateProofPolynomial:finalContraction");
    #endif
    vec_t<Factor<F>> factors = fg.getFactors();
    assert ( factors.size() == 1 );
    assert ( factors[0].isScalar() );

    #if TRACK_TIME >= 2
    logTimeEnd("evaluateProofPolynomial");
    #endif
    return factors[0].scalarValue();
}

/* NOTE: Inefficient */
template <typename F>
F evaluateProofPolynomial(
    const FactorGraph<F>& fgIn,
    const vec_t<Var>& cutset,
    const map_t<Var, map_t<uint_t,F>>& etas,
    const map_t<vec_t<uint_t>, F>& tau,
    const F zeta,
    const vec_t<IdTriple> order={})
{
    // Interpolate \hat{\ell}_i for i in cutset
    map_t<Var, Poly<F>> ells = interpolateElls(cutset, etas, tau);
    return evaluateProofPolynomial<F>(fgIn, cutset, etas, ells, zeta, order);
}
    
/* See equation 12 of Karimi, Kaski, and Koivisto (2020) */
template <typename F>
long_uint_t proofPolynomialDegreeBound(
    const FactorGraph<F>& fg,
    const vec_t<Var>& cutset)
{
    #if TRACK_TIME
    logTimeStart("proofPolynomialDegreeBound");
    #endif
    long_uint_t sum {0};
    for (Factor<F> f : fg.getFactors()) {
        for (Var v : f.getVars()) {
            if (isIn(cutset, v))
                sum += v.size - 1;
        }
    }
    long_uint_t cutsetVolume = vectorProduct(varsToSizes(cutset));
    #if TRACK_TIME
    logTimeEnd("proofPolynomialDegreeBound");
    #endif
    return (cutsetVolume - 1)*sum;
}

template <typename F>
bool verifyProofPolynomial(
    const Poly<F>& proof,
    const FactorGraph<F>& fg,
    const vec_t<Var>& cutset,
    const map_t<Var, map_t<uint_t,F>>& etas,
    const map_t<vec_t<uint_t>, F>& tau,
    uint_t numOfTrials=1)
{
    #if TRACK_TIME
    str_t degStr = std::to_string(proof.degree());
    logTimeStart("verifyProofPolynomial", degStr);
    #endif
    for (uint_t t = 0; t < numOfTrials; ++t) {
        F zeta = F::rand();
        F trueValue = evaluateProofPolynomial(fg, cutset, etas, tau, zeta);
        if (proof(zeta) != trueValue) {
            #if TRACK_TIME
            logTimeEnd("verifyProofPolynomial");
            #endif
            return false;
        }
    }
    #if TRACK_TIME
    logTimeEnd("verifyProofPolynomial", degStr);
    #endif
    return true;
}

template <typename F>
bool verifyProofPolynomial(
    const Poly<F>& proof,
    const FactorGraph<F>& fg,
    const vec_t<Var>& cutset,
    const map_t<Var, map_t<uint_t,F>>& etas,
    const map_t<vec_t<uint_t>, F>& tau,
    uint_t numOfTrials,
    vec_t<F>& trialInputs,
    vec_t<F>& trueValues,
    vec_t<F>& outputValues)
{
    #if TRACK_TIME
    logTimeStart("verifyProofPolynomial[var2]");
    #endif
    bool success = true;
    for (uint_t t = 0; t < numOfTrials; ++t) {
        F zeta = F::rand();
        trialInputs.push_back(zeta);
        F trueValue = evaluateProofPolynomial(fg, cutset, etas, tau, zeta);
        trueValues.push_back(trueValue);
        F outputValue = proof(zeta);
        outputValues.push_back(outputValue);
        if (outputValue != trueValue) {
            #if TRACK_TIME
            logTimeEnd("verifyProofPolynomial[var2]");
            #endif
            success = false;
        }
    }
    #if TRACK_TIME
    logTimeEnd("verifyProofPolynomial[var2]");
    #endif
    return success;
}

template <typename Z>
map_t<Var, map_t<uint_t, Z>> constructEtas(
    const vec_t<Var>& cutset,
    long_scalar_t maxZ = (Z::zero - Z::one).value())
{
    // Construct injective map eta_i for i in cutset
    map_t<Var, map_t<uint_t, Z>> etas;
    for (Var v : cutset) {
        // In order for the maps to be injective, the cardinality of Z
        // needs to be greater than the cardinality of any var's domain
        assert ( maxZ > v.size );
        map_t<uint_t, Z> eta;
        for (uint_t i = 0; i < v.size; ++i)
            eta[i] = Z(i);
        etas[v] = eta;
    }
    return etas;
}

template <typename Z>
map_t<Var, map_t<uint_t, Z>> constructDefaultEtas(
    const set_t<Var>& cutset,
    long_scalar_t maxZ = (Z::zero - Z::one).value())
{
    vec_t<Var> cutsetVec = setToVector<Var>(cutset);
    std::sort(cutsetVec.begin(), cutsetVec.end());
    return constructEtas<Z>(cutsetVec, maxZ);
}

//NOTE: Not a generator; builds a large map in memory; memory-inefficient.
template <typename Z>
map_t<vec_t<scalar_t>, Z> constructTau(
    const vec_t<Var>& cutset, 
    long_scalar_t maxZ = (Z::zero - Z::one).value())
{
    vec_t<vec_t<scalar_t>> tauDomain;
    enumerateVectors(varsToSizes(cutset), tauDomain);
    // In order for tau to be injective, and for proofs to be verifiable,
    // the cardinality of Z needs to be greater than the cardinality tauDomain
    assert ( maxZ > tauDomain.size() );
    map_t<vec_t<scalar_t>, Z> tau;
    for (scalar_t i = 0; i < tauDomain.size(); ++i) {
        tau[tauDomain[i]] = Z(i);
    }
    return tau;
}

template <typename Z>
map_t<vec_t<scalar_t>, Z> constructDefaultTau(
    const set_t<Var>& cutset, 
    long_scalar_t maxZ = (Z::zero - Z::one).value())
{
    vec_t<Var> cutsetVec = setToVector<Var>(cutset);
    std::sort(cutsetVec.begin(), cutsetVec.end());
    return constructTau<Z>(cutsetVec, maxZ);
}

template <typename F>
vec_t<F> defaultEvalPoints(const FactorGraph<F>& fg,
                           const vec_t<Var>& cutset,
                           uint_t numRedundantPoints=0)
{
    long_uint_t deg = proofPolynomialDegreeBound(fg, cutset);
    long_uint_t n = deg + 1 + numRedundantPoints;
    assert ( F::characteristic > n );
    vec_t<F> points(n);
    for (long_uint_t i = 0; i < n; ++i)
        points[i] = F(i);
    return points;
}

// NOTE: vectors in tau's domain are assumed to be in the same order as cutset
template <typename F>
F sumProofPolynomialOver(
    const vec_t<Var>& sumVars,
    const map_t<Var,scalar_t>& constVarValues,
    const Poly<F>& proofPoly,
    const vec_t<Var>& cutset,
    const map_t<vec_t<scalar_t>, F>& tau)
{
    assert ( isSuperset(cutset, sumVars) );
    if ( cutset.size() == 0 ) {
        assert ( proofPoly.degree() == 0 );
        return proofPoly[0];
    }
    
    vec_t<uint_t> sumVarIndices;
    for (Var v : sumVars)
        sumVarIndices.push_back(indexOf(v, cutset));
    
    vec_t<scalar_t> completeSlice(cutset.size());
    for (uint_t i = 0; i < cutset.size(); ++i) {
        if (isIn(sumVarIndices, i))
            completeSlice[i] = 0;
        else
            completeSlice[i] = constVarValues.at(cutset[i]);
    }
    
    if (sumVars.size() == 0)
        return proofPoly(tau.at(completeSlice));
    
    vec_t<scalar_t> ranges = varsToSizes(sumVars);
    vec_t<scalar_t> sumVarsSlice(sumVars.size(), 0);
    F sum {0};
    do {
        for (uint_t i = 0; i < sumVars.size(); ++i)
            completeSlice[sumVarIndices[i]] = sumVarsSlice[i];
        sum += proofPoly(tau.at(completeSlice));
    } while (nextCombination(sumVarsSlice, ranges));
    
    return sum;
}

template <typename F>
Factor<F> boundaryMapFromProof(
    const Poly<F>& proofPoly,
    const vec_t<Var>& cutset,
    const vec_t<Var>& boundary,
    const map_t<vec_t<scalar_t>, F>& tau,
    label_t factorLabel=0)
{
    #if TRACK_TIME
    str_t degStr = "deg=" + std::to_string(proofPoly.degree());
    logTimeStart("boundaryMapFromProof", degStr);
    #endif
    assert ( isSuperset(cutset, boundary) );
    vec_t<Var> sumVars = complement(cutset, boundary);
    vec_t<scalar_t> boundaryRanges = varsToSizes(boundary);
    vec_t<F> resultData;
    vec_t<scalar_t> boundarySlice(boundary.size(), 0);
    do {
        map_t<Var,scalar_t> boundarySliceMap;
        for (uint_t i = 0; i < boundarySlice.size(); ++i)
            boundarySliceMap[boundary[i]] = boundarySlice[i];
        resultData.push_back(
            sumProofPolynomialOver(
                sumVars, boundarySliceMap, proofPoly, cutset, tau));
    } while (nextCombination(boundarySlice, boundaryRanges));
    
    #if TRACK_TIME
    logTimeEnd("boundaryMapFromProof", degStr);
    #endif
    return Factor<F>(factorLabel, boundary, resultData.data());
}

template <typename F>
Factor<F> contractViaProof(
    const FactorGraph<F>& fg,
    const vec_t<Var>& cutset,
    const map_t<Var, map_t<uint_t, F>>& etas,
    const map_t<vec_t<scalar_t>, F>& tau,
    const map_t<Var, Poly<F>>& ells,
    Poly<F>& proofPoly,
    const vec_t<IdTriple>& order={},
    uint_t numRedundantPoints=0)
{
    #if TRACK_TIME
    logTimeStart("contractViaProof", "", true);
    #endif
    
    // Compute evaluations required to construct proof polynomial
    vec_t<F> evalPoints = defaultEvalPoints(fg, cutset, numRedundantPoints);
    long_uint_t e = evalPoints.size();
    vec_t<F> evals(e);
    for (long_uint_t i = 0; i < e; ++i)
        evals[i] = evaluateProofPolynomial(
            fg, cutset, etas, ells, evalPoints[i], order);
    
    // Decode the proof polynomial from the evaluations
    assert ( e < 1UL<<63 );
    long_uint_t d = proofPolynomialDegreeBound(fg, cutset);
    assert ( d < 1UL<<63 );
    assert ( rsDecodeXY((index_t)d, (index_t)e, 
                        evalPoints.data(), evals.data(), proofPoly) );
    
    // Recover result from proofPoly
    Factor<F> result = boundaryMapFromProof(
        proofPoly, cutset, setToVector(fg.getBoundary()), tau);
    
    #if TRACK_TIME
    logTimeEnd("contractViaProof", "", true);
    #endif
    return result;
}

template <typename F>
Factor<F> contractViaProof(const FactorGraph<F>& fg, const vec_t<Var>& cutset)
{
    map_t<Var, map_t<uint_t, F>> etas = constructEtas<F>(cutset);
    map_t<vec_t<scalar_t>, F> tau = constructTau<F>(cutset);
    map_t<Var, Poly<F>> ells = interpolateElls(cutset, etas, tau);
    Poly<F> proofPoly;
    return contractViaProof(fg, cutset, etas, tau, ells, proofPoly);
}

#endif  // __IDIGM_EVAL_ALGO_H__
