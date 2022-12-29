#ifndef __IDIGM_FACTORGRAPH_H__
#define __IDIGM_FACTORGRAPH_H__

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <memory>
#include <functional>
#include "montgomery.cpp"
#include "util.cpp"
#include "factor.cpp"

#define MAX_NUM_FACTORS_IN_GRAPH 1UL<<31
#define MAX_NUM_VARS_IN_GRAPH 1UL<<31


struct IdTriple {
    label_t first, second, newLabel;
};

std::ostream& operator<<(std::ostream& out, const IdTriple& t)
{
    out << "{" << t.first << ", " << t.second << ", " << t.newLabel << "}";
    return out;
}


template <typename T>
class FactorGraph {
private:
    map_t<label_t, std::shared_ptr<Factor<T>>> factors;
    set_t<Var> boundary;
    vec_t<Var> vars;
    
public:
    FactorGraph() : factors {}, boundary {}, vars {} {}
    
    FactorGraph(const set_t<Factor<T>> pfactors, set_t<Var> pboundary) :
        factors {}, boundary {pboundary}, vars {}
    {
        for (Factor<T> f : pfactors)
            insertFactor(f);
    }
    
    FactorGraph(const FactorGraph& other) :
        factors(other.factors), boundary(other.boundary), vars(other.vars) {}
    
    FactorGraph& operator=(const FactorGraph& other)
    {
        factors = other.factors;
        boundary = other.boundary;
        vars = other.vars;
        return *this;
    }
    
    uint_t nFactors() const { return factors.size(); }
    
    vec_t<Factor<T>> getFactors() const
    {
        vec_t<Factor<T>> fs;
        for (auto it = factors.begin(); it != factors.end(); ++it)
            fs.push_back(*(it->second));
        return fs;
    }
    
    vec_t<label_t> getFactorIds() const
    {
        vec_t<label_t> fIds;
        for (auto it = factors.begin(); it != factors.end(); ++it)
            fIds.push_back(it->first);
        return fIds;
    }
    
    Factor<T> getFactorWithId(label_t fId) const { return *factors.at(fId); }
    
    label_t createFactor(vec_t<Var> pvars, const T* pdata)
    {
        insertVars(pvars);
        label_t fId = getFreeFactorId();
        std::shared_ptr<Factor<T>> fptr(new Factor<T>(fId, pvars, pdata));
        factors.insert({fId, fptr});
        return fId;
    }
    
    label_t getFreeFactorId() const
    {
        vec_t<label_t> reservedIds = getFactorIds();
        label_t fId {0};
        long_uint_t tries {0}, maxTries {MAX_NUM_FACTORS_IN_GRAPH};
        while (isIn(reservedIds, fId)) {
            fId += 1;
            tries += 1;
            assert ( tries < maxTries );
        }
        return fId;
    }
    
    void insertFactor(const Factor<T>& f)
    {
        // Factor with identical label already exists?
        assert ( !(isIn(factors, f.label)) );
        insertVars(f.vars);
        std::shared_ptr<Factor<T>> fptr(new Factor<T>(f));
        factors.insert({f.label, fptr});
    }
    
    void removeFactor(label_t fId)
    {
        for (auto it = factors.begin(); it != factors.end(); ) {
            if (it->first == fId) {
                vec_t<Var> fVars = (it->second)->vars;
                (it->second).reset();
                it = factors.erase(it);
                for (Var v : fVars)
                    removeVarIfDetached(v);
            }
            else
                ++it;
        }
    }
    
    void replaceFactor(label_t fId, Factor<T> newFactor)
    {
        assert ( isIn(factors, fId) );
        Factor<T> f = *factors.at(fId);
        assert ( isSuperset(f.vars, newFactor.vars) );
        removeFactor(fId);
        insertFactor(newFactor);
    }
    
    set_t<Var> getBoundary() const { return boundary; }
    vec_t<label_t> getBoundaryIds() const { return varsToLabels(boundary); }
    void setBoundary(set_t<Var> b) { boundary = b; }
    
    uint_t nVars() const { return vars.size(); }
    vec_t<Var> getVars() const { return vars; }
    vec_t<label_t> getVarIds() const { return varsToLabels(vars); }
    
    label_t getFreeVarId() const
    {
        vec_t<label_t> reservedIds = getVarIds();
        reservedIds = unionVectors(reservedIds, varsToLabels(boundary));
        label_t vId {0};
        long_uint_t tries {0}, maxTries {MAX_NUM_VARS_IN_GRAPH};
        while (isIn(reservedIds, vId)) {
            vId += 1;
            tries += 1;
            assert ( tries < maxTries );
        }
        return vId;
    }
    
    void insertVars(const vec_t<Var>& pvars)
    {
        for (Var v : pvars) {
            if (isIn(vars, v))
                assert ( v.size == vars[(uint_t)indexOf(v,vars)].size );
            else
                vars.push_back(v);
        }
    }
    
    void insertVar(Var v) { insertVars({v}); }
    
    /* get Vars that are internal to the contraction of fId with gId */
    vec_t<Var> getInternalVars(label_t fId, label_t gId) const
    {
        vec_t<Var> internalVars = unionVectors(
            factors.at(fId) -> vars,
            factors.at(gId) -> vars
        );
        internalVars = complement(internalVars, setToVector(boundary));
        for (auto it = factors.begin(); it != factors.end(); ++it) {
            if ((it->first != fId) && (it->first != gId)) {
                internalVars = complement(internalVars,
                                                (it->second)->vars);
            }
        }
        return internalVars;
    }
    
    vec_t<label_t> getInternalVarIds(label_t fId, label_t gId) const
    { return varsToLabels(getInternalVars(fId, gId)); }
    
    Var createVar(scalar_t varSize)
    {
        Var v(getFreeVarId(), varSize);
        insertVars({v});
        return v;
    }
    
    // false if v is adjacent to a Factor, else true
    bool removeVarIfDetached(Var v)
    {
        for (auto it = factors.begin(); it != factors.end(); ++it) {
            if (isIn((it->second)->vars, v)) {
                return false;
            }
        }
        index_t i = indexOf(v,vars);
        if (i >= 0)  // detached but registered in graph
            vars.erase(vars.begin() + (uint_t)i);
        return true;
    }
    
    void replaceVarInFactor(label_t fId, Var oldVar, Var newVar)
    {
        assert ( oldVar.size == newVar.size );
        assert ( isIn(factors, fId) );
        Factor f = factors.at(fId) -> renameVar(oldVar.label, newVar.label);
        removeFactor(fId);
        insertFactor(f);
    }
    
    bool adjacent(label_t factorId, Var v)
    { return factors.at(factorId) -> adjacentToVar(v); }
    
    vec_t<label_t> adjacentFactorIds(Var v) const
    {
        vec_t<label_t> fIds;
        for (auto it = factors.begin(); it != factors.end(); ++it) {
            if ( isIn((it->second)->vars, v) )
                fIds.push_back(it->first);
        }
        return fIds;
    }
    
    FactorGraph<T> contractInOrder(vec_t<IdTriple> ordering) const
    {
        FactorGraph<T> fg(*this);
        for (IdTriple t : ordering) {
            vec_t<Var> internalVars = fg.getInternalVars(t.first, t.second);
            Factor<T> newF = contractFactors(
                *fg.factors.at(t.first), *fg.factors.at(t.second),
                internalVars, t.newLabel
            );
            fg.removeFactor(t.first);
            fg.removeFactor(t.second);
            fg.insertFactor(newF);
        }
        return fg;
    }
    
    // NOTE: This function is extremely inefficient
    FactorGraph<T> contractAll() const
    {
        FactorGraph<T> fg(*this);
        if (fg.nFactors() == 1) {
            // Add a "dummy" unit factor
            T* dataF = new T[1] {T(1)};
            Factor<T> f(fg.getFreeFactorId(), {fg.createVar(1)}, dataF);
            delete[] dataF;
            fg.insertFactor(f);
        }
        while (fg.nFactors() > 1) {
            label_t newFId = fg.getFreeFactorId();
            vec_t<label_t> existingFIds = fg.getFactorIds();
            IdTriple t {existingFIds[0], existingFIds[1], newFId};
            fg = fg.contractInOrder({t});
        }
        return fg;
    }
    
    FactorGraph<T> slice(map_t<Var, scalar_t> varsToInds)
    {
        FactorGraph<T> g(*this);
        for (auto it = g.factors.begin(); it != g.factors.end(); ++it) {
            (it->second).reset(
                new Factor<T>((it->second) -> slice(varsToInds, it->first)));
        }
        for (auto varIndPair: varsToInds) {
            bool removed = g.removeVarIfDetached(varIndPair.first);
            assert (removed);
        }
        return g;
    }
};

template <typename F>
FactorGraph<F> genRandFactorGraph(
    uint_t maxNFactors,
    uint_t maxNVarsPerFactor,
    uint_t maxVarSize,
    vec_t<label_t> varIds = {},
    bool randomFactorLabels=false,
    uint_t minNFactors=1,
    uint_t minVarSize=2,
    F (*genRandScalar)() = nullptr,
    uint_t minBoundarySize = 0,
    uint_t maxBoundarySize = -1U)
{
    if (varIds.size() == 0) {
        uint_t nVars = maxNVarsPerFactor * (1 + maxNFactors/5);
        for (label_t l = 0; l < nVars; ++l)
            varIds.push_back(l + maxNFactors*1000);
    }
    FactorGraph<F> fg;
    uint_t nFactors = minNFactors;
    nFactors += ((uint_t)rand()) % (maxNFactors - minNFactors + 1);
    map_t<label_t,uint_t> varSizes;
    for (label_t n = 0; n < nFactors; ++n) {
        std::pair<Factor<F>, map_t<label_t,uint_t>>
        fp = genRandFactor<F>(maxNVarsPerFactor, maxVarSize, varIds,
                              varSizes, minVarSize, genRandScalar);
        Factor<F> f = fp.first;
        if (!randomFactorLabels)
            f.label = n;
        fg.insertFactor(f);
        varSizes = fp.second;
    }
    
    vec_t<Var> vars = fg.getVars();
    uint_t nv = vars.size();
    assert (nv >= minBoundarySize);
    uint_t nb = minBoundarySize;
    nb += (uint_t)rand() % (nv - nb + 1);
    nb = nb > maxBoundarySize ? maxBoundarySize : nb;
    fg.setBoundary(vectorToSet(randSubset(vars, nb)));
    
    return fg;
}

// NOTE: cutsetMaxSize is ignored if is is less than the boundary's size;
// since the cutset must be a superset of the boundary.
template <typename T>
vec_t<Var> genRandCutset(const FactorGraph<T>& fg, uint_t cutsetMaxSize=-1U)
{
    vec_t<Var> cutset = setToVector(fg.getBoundary());
    uint_t nc = cutset.size();
    if (nc < fg.nVars() && nc < cutsetMaxSize) {
        vec_t<Var> availableVars = complement(fg.getVars(), cutset);
        uint_t nToAdd = (uint_t)rand() % availableVars.size();
        nToAdd = nToAdd + nc > cutsetMaxSize ? cutsetMaxSize - nc : nToAdd;
        cutset = unionVectors(cutset, randSubset(availableVars, nToAdd));
    }
    return cutset;
}

template <typename F>
void printFactorGraph(const FactorGraph<F>& fg)
{
    std::cout << "Boundary: "; printSet(fg.getBoundary());
    std::cout << "Vars:\n"; printVector(fg.getVars());
    std::cout << "Factors:\n";
    for (Factor<F> f: fg.getFactors()) {
        std::cout << "Factor " << f.label << " : ";
        printFactor(f);
    }
}

template <typename F>
void printFactorAdjacencies(const FactorGraph<F>& fg)
{
    for (auto fPair : fg.factors) {
        std::cout << "Factor " << fPair.first << " : \n    ";
        printVector(fPair.second -> vars);
    }
}


template <typename F>
vec_t<Var> adjacentVars(const FactorGraph<F>& fg,
                        const vec_t<label_t>& factorIds)
{
    vec_t<Var> vars;
    for (label_t fId : factorIds)
        vars = unionVectors(vars, fg.getFactorWithId(fId).getVars());
    return vars;
}

template <typename F>
void simulateContraction(FactorGraph<F>& fg, IdTriple c)
{
    vec_t<Var> internalVars = fg.getInternalVars(c.first, c.second);
    vec_t<Var> remainingVars = complement(
        adjacentVars(fg, {c.first, c.second}),
        internalVars
    );
    fg.removeFactor(c.first);
    fg.removeFactor(c.second);
    Factor<F> dummyFactor(c.newLabel, remainingVars, std::shared_ptr<F[]>());
    fg.insertFactor(dummyFactor);
}

template <typename F>
void simulateContractions(FactorGraph<F>& fg, const vec_t<IdTriple>& order)
{
    for (IdTriple c : order)
        simulateContraction(fg, c);
}

template <typename F>
long_uint_t contractionCost(const FactorGraph<F>& fg,
                            label_t fId1, label_t fId2)
{
    return vectorProduct(varsToSizes(adjacentVars(fg,{fId1,fId2})));
}

template <typename F>
long_uint_t contractionCost(const FactorGraph<F>& fg, vec_t<IdTriple> order)
{
    uint_t n = order.size();
    if (n == 0)
        return 0;
    
    FactorGraph<F> dummyFG = fg;
    long_uint_t cost = 0;
    for (uint_t i = 0; i < n; ++i) {
        cost += contractionCost(dummyFG, order[i].first, order[i].second);
        simulateContraction(dummyFG, order[i]);
    }
    return cost;
}

size_t labelSetHash(vec_t<label_t> labels)
{
    std::sort(labels.begin(), labels.end());
    str_t s = "";
    for (label_t label : labels)
        s += std::to_string(label);
    return std::hash<str_t>()(s);
}

// NOTE: Exponential complexity in number of factors
template <typename F>
vec_t<IdTriple> optimalContractionOrder(
    const FactorGraph<F>& fg,
    vec_t<label_t> factorIds,
    map_t< size_t, long_uint_t >& optimalCosts,
    map_t< size_t, vec_t<IdTriple> >& optimalOrders
)
{
    #if TRACK_TIME
    logTimeStart("optimalContractionOrder");
    #endif
    
    vec_t<label_t> allFactorIds = fg.getFactorIds();
    assert ( allFactorIds.size() < 63 );
    assert ( isSuperset(allFactorIds, factorIds) );
    uint_t n = factorIds.size();
    assert ( n > 0 );
    if (n == 1) {
        size_t hash = labelSetHash(factorIds);
        optimalCosts[hash] = 0;
        optimalOrders[hash] = {};
        #if TRACK_TIME
        logTimeEnd("optimalContractionOrder");
        #endif
        return {};
    }
    if (n == 2) {
        size_t hash = labelSetHash(factorIds);
        vec_t<IdTriple> order {{factorIds[0], factorIds[1], factorIds[1]}};
        if ( ! isIn(optimalCosts, hash)) {
            optimalCosts[hash] = contractionCost(fg, order);
            optimalOrders[hash] = order;
        }
        #if TRACK_TIME
        logTimeEnd("optimalContractionOrder");
        #endif
        return order;
    }
    
    for (index_t i = 1; i < (1L<<n)-1; ++i) {
        FactorGraph<F> dummyFG = fg;
        
        vec_t<label_t> leftSubset = subset(factorIds, i);
        vec_t<label_t> rightSubset = complement(factorIds, leftSubset);
        
        // compute optimal cost and contraction order for these subsets
        optimalContractionOrder(fg, leftSubset, optimalCosts, optimalOrders);
        optimalContractionOrder(fg, rightSubset, optimalCosts, optimalOrders);
        
        size_t leftHash = labelSetHash(leftSubset);
        size_t rightHash = labelSetHash(rightSubset);
        vec_t<IdTriple> leftOrder = optimalOrders[leftHash];
        vec_t<IdTriple> rightOrder = optimalOrders[rightHash];
        simulateContractions(dummyFG, leftOrder);
        simulateContractions(dummyFG, rightOrder);
        #if DEBUG_ASSERTS
        assert ( leftOrder.size() >= 1 || leftSubset.size() == 1 );
        assert ( rightOrder.size() >= 1 || rightSubset.size() == 1 );
        #endif
        label_t leftFId = leftOrder.size() > 0 ?
            leftOrder.back().newLabel : leftSubset[0];
        label_t rightFId = rightOrder.size() > 0 ?
            rightOrder.back().newLabel : rightSubset[0];
        long_uint_t cost = contractionCost(dummyFG, leftFId, rightFId)
                           + optimalCosts[leftHash]
                           + optimalCosts[rightHash];
        size_t hash = labelSetHash(factorIds);
        if ( !isIn(optimalCosts, hash) || optimalCosts.at(hash) > cost) {
            optimalCosts[hash] = cost;
            vec_t<IdTriple> order = leftOrder;
            order.insert(order.end(), rightOrder.begin(), rightOrder.end());
            order.push_back( {leftFId, rightFId, dummyFG.getFreeFactorId()} );
            optimalOrders[hash] = order;
        }
    }
    
    #if TRACK_TIME
    logTimeEnd("optimalContractionOrder");
    #endif
    return optimalOrders[labelSetHash(factorIds)];
}


#endif // __IDIGM_FACTORGRAPH_H__

