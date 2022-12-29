#ifndef __IDIGM_FACTOR_H__
#define __IDIGM_FACTOR_H__

#include <iostream>
#include <vector>
#include <cstring>
#include <map>
#include <utility>
#include <exception>
#include <memory>
#include "util.cpp"


typedef long unsigned int label_t;
typedef long unsigned int long_scalar_t;


struct Var {
    label_t label;
    uint_t size;
    
    Var() : label {0}, size {0} {}
    Var(label_t l, uint_t s) : label{l}, size{s} {}
    
    friend bool operator==(const Var& l, const Var& r) {
        if (l.label == r.label) {
            assert ( l.size == r.size );
            return true;
        }
        return false;
    }
    
    /* Needed for std::map<Var,_>; not meaningful as a comparison relation */
    friend bool operator<(const Var& l, const Var& r) {
        return l.label < r.label;
    }
};

std::ostream& operator<<(std::ostream& out, const Var& v)
{
    return out << "Var_" << v.label << "[" << v.size << "]";
}

vec_t<uint_t> varsToSizes(vec_t<Var> vars)
{
    uint_t nVars = vars.size();
    vec_t<uint_t> sizes(nVars);
    for (uint_t i = 0; i < nVars; ++i)
        sizes[i] = vars[i].size;
    return sizes;
}

vec_t<label_t> varsToLabels(vec_t<Var> vars)
{
    uint_t nVars = vars.size();
    vec_t<label_t> labels(nVars);
    for (uint_t i = 0; i < nVars; ++i)
        labels[i] = vars[i].label;
    return labels;
}

vec_t<label_t> varsToLabels(set_t<Var> varsSet)
{
    vec_t<Var> vars = setToVector<Var>(varsSet);
    uint_t nVars = vars.size();
    vec_t<label_t> labels(nVars);
    for (uint_t i = 0; i < nVars; ++i)
        labels[i] = vars[i].label;
    return labels;
}


template <typename T>
class Factor {
public:
    label_t label;
    vec_t<Var> vars;
    long_scalar_t volume;
    std::shared_ptr<T[]> data;

    Factor() :
        label{0}, vars {vec_t<Var>()}, volume{0}, data() {};
    
    Factor(label_t plabel, vec_t<Var> pvars, const T* pdata) :
        label {plabel}, vars {pvars},
        volume {vectorProduct(varsToSizes(vars))},
        data(new T[volume])
    {
        #if DEBUG_ASSERTS
        assert ( volume > 0 && volume < 1UL<<63);
        #endif
        arrayCopy((index_t)volume, pdata, data.get());
    }
    
    Factor(label_t plabel, vec_t<Var> pvars, std::shared_ptr<T[]> pdata) :
        label {plabel}, vars {pvars},
        volume {vectorProduct(varsToSizes(vars))},
        data(pdata)
    {
        #if DEBUG_ASSERTS
        assert ( volume > 0 );
        #endif
    }
    
    Factor(const Factor& other) :
        label {other.label}, vars {other.vars}, volume {other.volume},
        data(other.data)
    {
        #if DEBUG_ASSERTS
        assert ( volume > 0 );
        assert ( volume == vectorProduct(varsToSizes(vars)) );
        #endif
    }
    
    ~Factor()
    {
        data.reset();
    }
    
    Factor& operator= (const Factor& other)
    {
        label = other.label;
        vars = other.vars;
        volume = other.volume;
        data = other.data;
        return *this;
    }
    
    friend bool operator< (const Factor& lhs, const Factor& rhs)
    {
        return lhs.label < rhs.label;
    }
    
    const T* getData() const { return data.get(); }
    T* dataAt(long_scalar_t pos) { return data.get() + pos; }
    T valueAt(long_scalar_t pos) const { return *(data.get() + pos); }
    
    bool isScalar() const { return volume == 1; }
    
    T scalarValue() const
    {
        assert (volume == 1);
        return data.get()[0];
    }
    
    label_t getId() const { return label; }
    
    long_scalar_t getVolume() const { return volume; }
    
    uint_t nVars() const { return vars.size(); }
    
    vec_t<Var> getVars() const { return vars; }
    
    vec_t<uint_t> varSizes() const { return varsToSizes(vars); }
    
    vec_t<uint_t> varSizesWithout(vec_t<Var> excludedVars) const
    { return varsToSizes(complement(vars, excludedVars)); }
    
    vec_t<label_t> varIds() const { return varsToLabels(vars); }
    
    Factor<T> renameVar(label_t oldId, label_t newId) const
    {
        vec_t<Var> newVars(vars);
        uint_t n = newVars.size();
        bool success {false};
        for (uint_t i = 0; i < n; ++i) {
            if (newVars[i].label == oldId) {
                newVars[i].label = newId;
                success = true; break;
            }
        }
        assert (success);
        return Factor<T>(label, newVars, data);
    }
    
    uint_t sizeOfVar(label_t vId) const
    {
        for (Var v : vars) {
            if (v.label == vId)
                return v.size;
        }
        throw std::invalid_argument(
            "sizeOfVar: Var " + std::to_string(vId)
            + " not found in Factor " + std::to_string(label) + ".\n");
    }
    
    Var getVarWithId(label_t vId) const
    {
        for (Var v : vars) {
            if (v.label == vId)
                return v;
        }
        throw std::invalid_argument(
            "getVarWithId: Var " + std::to_string(vId)
            + " not found in Factor " + std::to_string(label) + ".\n");
    }
    
    bool adjacentToVar(label_t vId) const
    { return isIn(varsToLabels(vars), vId); }
    
    bool adjacentToVar(Var v) const
    { return isIn(vars, v); }

    Factor<T> slice(map_t<Var, uint_t> varsToInds, label_t newFId)
    {
        map_t<uint_t, uint_t> dimsToInds;
        long_scalar_t resultVol = volume;
        vec_t<Var> resultVars(vars);
        auto it = resultVars.begin();
        for (scalar_t i = 0; i < vars.size(); ++i) {
            if (isIn(varsToInds, vars[i])) {
                dimsToInds[i] = varsToInds[vars[i]];
                resultVol /= vars[i].size;
                it = resultVars.erase(it);
            }
            else
                ++it;
        }
        T* resultData = new T[resultVol];
        arraySlice(dimsToInds, varSizes(), data.get(), resultData);
        Factor<T> result(newFId, resultVars, resultData);
        delete[] resultData;
        return result;
    }
    
    Factor<T> permuteToEnd(vec_t<Var> endVars) const
    {
        vec_t<Var> permutedVars = moveToEnd(vars, endVars);
        vec_t<uint_t> perm = getPermutation(vars, permutedVars);
        return permuteFactor(perm.data(), *this);
    }
};

template <typename T>
void printFactor(const Factor<T>& f, bool printData=true)
{
        std::cout << "Factor[" << f.label << "]:\n";
        std::cout << "  .nVars(): " << f.nVars() << "\n";
        std::cout << "  .vars: "; printVector(f.vars);
        std::cout << "  .volume: " << f.volume << "\n";
        if (printData) {
            std::cout << "  .data: ";
            printNDArray<T>(f.nVars(), f.varSizes().data(), f.data.get());
        }
}

template <typename T>
std::pair<Factor<T>, map_t<label_t, scalar_t>>
genRandFactor(scalar_t maxNVars, scalar_t maxVarSize,
              vec_t<label_t> varLabels,
              map_t<label_t, scalar_t> labelsToSizes={},
              scalar_t minVarSize=2,
              T (*randScalar)() = nullptr)
{
    assert ( maxNVars <= varLabels.size() );
    if (randScalar == nullptr)
        randScalar = &genRandScalar<T>;
    
    uint_t nVars = ((uint_t)rand() % maxNVars) + 1;
    vec_t<label_t> varIds = randSubset(varLabels, nVars);
    vec_t<Var> vars(nVars);
    for (uint_t d = 0; d < nVars; ++d) {
        uint_t sizeOfVarD;
        bool alreadyDefined = isIn(labelsToSizes, varIds[d]);
        if (alreadyDefined)
            sizeOfVarD = labelsToSizes.at(varIds[d]);
        else {
            sizeOfVarD = minVarSize;
            sizeOfVarD += (uint_t)rand() % (maxVarSize - minVarSize + 1);
            labelsToSizes[varIds[d]] = sizeOfVarD;
        }
        vars[d] = Var(varIds[d], sizeOfVarD);
    }
    
    long_scalar_t vol = vectorProduct(varsToSizes(vars));
    T* data = new T[vol];
    for (long_scalar_t k = 0; k < vol; ++k)
        data[k] = randScalar();
    
    label_t label { (label_t) rand() };
    Factor<T> f(label, vars, data);
    delete[] data;
    
    return {f, labelsToSizes};
}

template <typename T>
Factor<T> genRandFactor(scalar_t maxNVars, scalar_t maxVarSize,
                        vec_t<label_t> varLabels, scalar_t minVarSize=2)
{
    return genRandFactor<T>(
        maxNVars, maxVarSize, varLabels, {}, minVarSize).first;
}


template <typename T>
Factor<T> permuteFactor(const scalar_t* perm, const Factor<T>& f)
{
    #if DEBUG_ASSERTS
    // TODO: check that the values in range [perm, perm + f.nVars()] contains
    // exactly the numbers {0, 1, ..., f.nVars()-1} in some order
    #endif
    #if TRACK_TIME >= 4
    str_t vol = std::to_string(f.volume);
    logTimeStart("permuteFactor", vol);
    #endif
    
    vec_t<Var> varsG = permuteVector(f.nVars(), perm, f.getVars());
    T* dataG = new T[f.volume];
    permuteArray<T>(f.nVars(), perm, f.varSizes().data(), f.data.get(), dataG);
    Factor<T> g(f.label, varsG, dataG);
    delete[] dataG;
    
    #if TRACK_TIME >= 4
    logTimeEnd("permuteFactor", vol);
    #endif
    return g;
}

template <typename T>
Factor<T> permuteFactor(const vec_t<Var>& targetVarOrder, const Factor<T>& f)
{
    return permuteFactor<T>(getPermutation(f.vars, targetVarOrder).data(), f);
}


template <typename T>
Factor<T> sumOverVar(Var v, const Factor<T>& f, label_t newFId)
{
    #if TRACK_TIME >= 3
    logTimeStart("sumOverVar", std::to_string(v.size));
    #endif
    // construct dataFPerm: copy of f.data.get() permuted so that v is last
    
    uint_t nVars = f.nVars();
    index_t vIndSigned = indexOf(v, f.vars);
    assert ( (vIndSigned >= 0) && (vIndSigned < nVars) );
    uint_t vInd = (uint_t) vIndSigned;
    
    T* dataFPerm = new T[f.volume];
    
    if (vInd == nVars - 1) {  // v was already last in f
        memcpy((void*)dataFPerm, (void*)f.data.get(), f.volume * sizeof(T));
    }
    else {
        uint_t* perm = new uint_t[nVars];
        for (uint_t i = 0; i < nVars; ++i) {
            if (i == nVars - 1)
                perm[i] = vInd;
            else if (i >= vInd)
                perm[i] = i + 1;
            else // i < vInd
                perm[i] = i;
        }
        permuteArray<T>(
            nVars, perm, f.varSizes().data(), f.data.get(), dataFPerm);
        delete[] perm;
    }
    
    // construct dataG: copy of dataFPerm, summed over last dimension
    
    uint_t vSize {v.size};
    long_scalar_t volG = f.volume / vSize;
    T* dataG = new T[volG];
    for (long_scalar_t k = 0; k < volG; ++k) {
        T sum = 0L;
        for (long_scalar_t j = k*vSize; j < k*vSize + vSize; ++j)
            sum += dataFPerm[j];
        dataG[k] = sum;
    }
    delete[] dataFPerm;
    
    vec_t<Var> varsG (f.vars);
    varsG.erase(varsG.begin() + vInd);
    
    Factor<T> g(newFId, varsG, dataG);
    delete[] dataG;
    #if TRACK_TIME >= 3
    logTimeEnd("sumOverVar", std::to_string(v.size));
    #endif
    return g;
}

template <typename T>
Factor<T> sumOverVar(label_t vId, const Factor<T>& f, label_t newFId)
{
    Var v = f.getVarWithId(vId);
    return sumOverVar(v, f, newFId);
}

template <typename T>
Factor<T> contractFactors(const Factor<T>& f, const Factor<T>& g,
                          vec_t<Var> internalVars,
                          label_t newFactorId)
{
    #if TRACK_TIME >= 3
    str_t cost = std::to_string(
        vectorProduct(varsToSizes(unionVectors(f.vars, g.vars))));
    logTimeStart("contractFactors", cost);
    #endif

    // Find singly-connected contraction-internal vars (SCIVs) for f and g
    
    vec_t<Var> internsF = intersect(f.vars, internalVars);
    vec_t<Var> scivsF = complement(internsF, g.vars);
        
    vec_t<Var> internsG = intersect(g.vars, internalVars);
    vec_t<Var> scivsG = complement(internsG, f.vars);
        
    // For each SCIV s, sum the adjacent Factor over s
    
    Factor<T> fs = f;
    for (Var s : scivsF)
        fs = sumOverVar(s, fs, fs.label);
    
    Factor<T> gs = g;
    for (Var s : scivsG)
        gs = sumOverVar(s, gs, gs.label);
    
    // Permute fs and gs:
    // arrange contraction-internal vars last and in same order,
    // arrange shared contraction-external vars first and in same order
    
    vec_t<Var> commonInterns = intersect(internsF, internsG);
    vec_t<Var> commonExterns = complement(
        intersect(fs.vars, gs.vars),
        internalVars);
    
    vec_t<Var> permdVarsF = moveToEnd(fs.vars, commonInterns);
    permdVarsF = moveToFront(permdVarsF, commonExterns);
    vec_t<uint_t> permF = getPermutation(fs.vars, permdVarsF);
    fs = permuteFactor(permF.data(), fs);
    
    vec_t<Var> permdVarsG = moveToEnd(gs.vars, commonInterns);
    permdVarsG = moveToFront(permdVarsG, commonExterns);
    vec_t<uint_t> permG = getPermutation(gs.vars, permdVarsG);
    gs = permuteFactor(permG.data(), gs);
    
    #if DEBUG_ASSERTS
    //assert ( fs.vars == permdVarsF );
    //assert ( gs.vars == permdVarsG );
    #endif
    
    // Construct vars of Factor h resulting from contraction
    
    vec_t<Var> externsF = complement(fs.vars, internalVars);
    externsF = complement(externsF, commonExterns);
    vec_t<Var> externsG = complement(gs.vars, internalVars);
    externsG = complement(externsG, commonExterns);
    
    vec_t<Var> varsH;
    varsH.insert(varsH.begin(), commonExterns.begin(), commonExterns.end());
    varsH.insert(varsH.end(), externsF.begin(), externsF.end());
    varsH.insert(varsH.end(), externsG.begin(), externsG.end());
    
    // construct data of Factor h via matrix products of fs and gs
    
    vec_t<scalar_t> varSizesH = varsToSizes(varsH);
    long_scalar_t volH { vectorProduct(varSizesH) }; 
    T* dataH = new T[volH];
    vec_t<scalar_t> commonExternSizes = varsToSizes(commonExterns);
    long_scalar_t volCommonExterns { vectorProduct(commonExternSizes) };
    scalar_t nf { vectorProduct(varsToSizes(externsF)) };
    scalar_t ng { vectorProduct(varsToSizes(externsG)) };
    scalar_t nc { vectorProduct(varsToSizes(commonInterns)) };
    
    #if DEBUG_ASSERTS
    assert ( fs.volume / (nf * nc) == volCommonExterns );
    assert ( gs.volume / (ng * nc) == volCommonExterns );
    assert ( volH == volCommonExterns * nf * ng);
    #endif
    
    T* subArrF = fs.dataAt(0);
    T* subArrG = gs.dataAt(0);
    T* subArrH = dataH;
    scalar_t nfnc = nf * nc;
    scalar_t ngnc = ng * nc;
    scalar_t nfng = nf * ng;
    #if TRACK_TIME >= 3
    str_t volStr = std::to_string(volCommonExterns); 
    logTimeStart("contractFactors:matmulloop", volStr);
    #endif
    for (long_scalar_t k = 0; k < volCommonExterns; ++k) {
        naiveMatrixMul(nf, nc, ng, subArrF, subArrG, subArrH);
        subArrF += nfnc;
        subArrG += ngnc;
        subArrH += nfng;
    }
    #if TRACK_TIME >= 3
    logTimeEnd("contractFactors:matmulloop", volStr);
    #endif
    
    Factor<T> h(newFactorId, varsH, dataH);
    delete[] dataH;

    #if TRACK_TIME >= 3
    logTimeEnd("contractFactors", cost);
    #endif
    return h;
}


template <typename T>
Factor<T> contractFactors(const Factor<T>& f, const Factor<T>& g,
                          vec_t<label_t> internalVarIds,
                          label_t newFactorId)
{
    vec_t<Var> internalVars;
    for (label_t vId : internalVarIds) {
        bool inF {false}, inG {false};
        if (f.adjacentToVar(vId)) {
            internalVars.push_back(f.getVarWithId(vId));
            inF = true;
        }
        else if (g.adjacentToVar(vId)) {
            internalVars.push_back(g.getVarWithId(vId));
            inG = true;
        }
        else {
            // TODO? Log a warning
            #if STRICT
            throw std::invalid_argument(
            "contractFactors: internal Var with label " + std::to_string(vId)
            + " not found in Factor " + std::to_string(f.label)
            + " or Factor " + std::to_string(g.label) + ".\n");
            #endif
        }
        if (inF && inG)
            assert ( f.sizeOfVar(vId) == g.sizeOfVar(vId) );
    }
    return contractFactors(f, g, internalVars, newFactorId);
}

template class Factor<int>;

#endif  // __IDIGM_FACTOR_H__

