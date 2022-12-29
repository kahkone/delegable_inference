#include <iostream>
#include <cstring>
#include "montgomery.cpp"
#include "util.cpp"
#include "fileio.cpp"
#include "factorgraph.cpp"
#include "crt.cpp"
#include "eval_algo.cpp"
#include "poly.cpp"


class Options {
public:
    bool keepIntermediateFiles;
    vec_t<str_t> outputFilenames;
    Options() : keepIntermediateFiles{false}, outputFilenames() {}
    Options(bool k, vec_t<str_t> o) :
        keepIntermediateFiles{k}, outputFilenames(o) {}
};


template <typename F>
str_t generateFG(str_t paramsFile)
{
    FileHandler fh(paramsFile);
    
    // Mandatory parameters
    uint_t maxNFactors;
    loadParamFromFile(fh, "maxNFactors", maxNFactors);
    uint_t maxNVarsPerFactor;
    loadParamFromFile(fh, "maxNVarsPerFactor", maxNVarsPerFactor);
    uint_t maxVarSize;
    loadParamFromFile(fh, "maxVarSize", maxVarSize);
    
    // Optional parameters
    vec_t<label_t> varIds = {};
    fh.toBeginning();
    if (fh.toNextOccurrence("varIds"))
        assert ( loadObject(fh, varIds) );
    bool randomFactorLabels = false;
    loadParamFromFile(fh, "randomFactorLabels", randomFactorLabels);
    uint_t minNFactors = 1;
    loadParamFromFile(fh, "minNFactors", minNFactors);
    uint_t minVarSize = 2;
    loadParamFromFile(fh, "minVarSize", minVarSize);
    uint_t cutsetMaxSize = -1U;
    loadParamFromFile(fh, "cutsetMaxSize", cutsetMaxSize);
    
    FactorGraph<F> g = genRandFactorGraph<F>(
        maxNFactors, maxNVarsPerFactor, maxVarSize, varIds,
        randomFactorLabels, minNFactors, minVarSize
    );
    vec_t<Var> cutset = genRandCutset<F>(g, cutsetMaxSize);
    
    str_t outputFilename = "generated_factorgraph_" + timestamp();
    FileHandler fhOut(outputFilename);
    writeFGToFile<F>(fhOut, g);
    writeParamToFile(fhOut, "cutset", cutset);
    return outputFilename;
}

template <typename P>
str_t crtReduce(str_t mpIntFGFile)
{
    FileHandler fh(mpIntFGFile);
    FactorGraph<MPInt> g;
    assert ( loadFactorGraph<MPInt>(fh, g) );
    scalar_t modulus = P::modulus;
    FactorGraph<Zp<P>> gReduced = mpIntToModFG<P>(g);
    str_t reducedFGFile = "factorgraph_reduced_mod_"
                          + std::to_string(modulus) + "_"
                          + timestamp();
    FileHandler fhr(reducedFGFile);
    writeParamToFile<str_t>(fhr, "type", str_t("modular"));
    writeParamToFile<scalar_t>(fhr, "modulus", modulus);
    writeFGToFile(fhr, gReduced);
    return reducedFGFile;
}

template <typename F>
bool loadOrConstructEtas(FileHandler& h,
                         const str_t& filename,
                         const vec_t<Var>& cutset,
                         map_t<Var, map_t<uint_t, F>>& etas)
{
    bool etasDefined = loadParamFromFile(h, "etas", etas);
    if ( ! etasDefined) {
        std::cout << "Info: Could not find definitions for injective maps "
            << "eta_i in file " << filename << "; generating default maps.\n";
        etas = constructDefaultEtas<F>(vectorToSet(cutset));
    }
    return etasDefined;
}

template <typename F>
bool loadOrConstructTau(FileHandler& h,
                        const str_t& filename,
                        const vec_t<Var>& cutset,
                        map_t<vec_t<uint_t>, F>& tau)
{
    bool tauDefined = loadParamFromFile(h, "tau", tau);
    if ( ! tauDefined) {
        std::cout << "Info: Could not find definition for injective map tau "
            << "in file " << filename << "; generating default tau.\n";
        tau = constructDefaultTau<F>(vectorToSet(cutset));
    }
    return tauDefined;
}

template <typename F>
str_t computeEvaluations(str_t workFile)
{
    FileHandler h(workFile);
    
    FactorGraph<F> fg;
    assert ( loadFactorGraph<F>(h, fg) );
    
    vec_t<Var> cutset;
    assert ( loadParamFromFile(h, "cutset", cutset) );
    std::sort(cutset.begin(), cutset.end());
    
    vec_t<F> evalPoints;
    bool evalPointsDefined = loadParamFromFile(h, "evalPoints", evalPoints);
    if ( ! evalPointsDefined) {
        uint_t nRedundantEvals {0};
        loadParamFromFile(h, "redundantEvaluations", nRedundantEvals);
        evalPoints = defaultEvalPoints(fg, cutset, nRedundantEvals);
        std::cout << "Warning: Could not find evaluation points in file "
                  << workFile << "; generated " << evalPoints.size()
                  << " unique evaulation points (" << nRedundantEvals
                  << " redundant points).\n";
    }
    
    long_uint_t proofPolyDegree;
    bool degreeDefined = loadParamFromFile(h, "degree", proofPolyDegree);
    if ( ! degreeDefined)
        proofPolyDegree = proofPolynomialDegreeBound(fg, cutset);
    
    map_t<Var, map_t<uint_t, F>> etas;
    loadOrConstructEtas(h, workFile, cutset, etas);
    
    map_t<vec_t<uint_t>, F> tau;
    loadOrConstructTau(h, workFile, cutset, tau);
    
    map_t<Var, Poly<F>> ells = interpolateElls<F>(cutset, etas, tau);
    
    vec_t<F> evaluations;
    for (F z : evalPoints)
        evaluations.push_back(
            evaluateProofPolynomial<F>(fg, cutset, etas, ells, z));
    
    str_t evalsFile = "evaluations_" + timestamp();
    FileHandler hOut(evalsFile);
    writeParamToFile(hOut, "degree", proofPolyDegree);
    uint_t modulus = F::characteristic;
    writeParamToFile(hOut, "modulus", modulus);
    writeParamToFile(hOut, "evalPoints", evalPoints);
    writeParamToFile(hOut, "evaluations", evaluations);
    
    return evalsFile;
}

template <typename F>
str_t decodeProofCoefficients(str_t evalsFile)
{
    FileHandler h(evalsFile);
    
    vec_t<F> evals;
    assert ( loadParamFromFile(h, "evaluations", evals) );
    vec_t<F> inputs;
    assert ( loadParamFromFile(h, "evalPoints", inputs) );
    assert ( evals.size() == inputs.size() );
    assert ( inputs.size() < 1UL<<63 );
    long e = (long)evals.size();
    long d;
    assert ( loadParamFromFile(h, "degree", d) );
    
    Poly<F> proofPoly;
    if (!rsDecodeXY(d, e, inputs.data(), evals.data(), proofPoly)) {
        std::cout << "Error: Failed to decode proof polynomial from data in "
            << evalsFile << ".\n";
        return "";
    }
    assert ( proofPoly.capacity() > 0 );
    vec_t<F> coeffs((long_uint_t)proofPoly.capacity());
    proofPoly.copyCoeffsTo(coeffs.data());
    
    str_t coeffsFile = "proof_coefficients_" + timestamp();
    FileHandler hOut(coeffsFile);
    writeParamToFile(hOut, "type", str_t("modular"));
    uint_t m = F::characteristic;
    writeParamToFile(hOut, "modulus", m);
    writeParamToFile(hOut, "coefficients", coeffs);
    
    return coeffsFile;
}

template <typename F>
str_t verifyProof(str_t fgFile, str_t proofFile)
{
    FileHandler h1(proofFile);
    vec_t<F> coeffs;
    assert ( loadParamFromFile(h1, "coefficients", coeffs) );
    Poly<F> proof(coeffs);
    
    FileHandler h2(fgFile);
    
    FactorGraph<F> fg;
    assert ( loadFactorGraph<F>(h2, fg) );
    
    vec_t<Var> cutset;
    assert ( loadParamFromFile(h2, "cutset", cutset) );
    std::sort(cutset.begin(), cutset.end());
    
    map_t<Var, map_t<uint_t, F>> etas;
    loadOrConstructEtas(h2, fgFile, cutset, etas);
    
    map_t<vec_t<uint_t>, F> tau;
    loadOrConstructTau(h2, fgFile, cutset, tau);
    
    uint_t nTrials {1};
    if ( ! loadParamFromFile(h1, "verificationTrials", nTrials) ) {
        if ( ! loadParamFromFile(h2, "verificationTrials", nTrials) ) {
            std::cout << "Warning: Could not find integer parameter "
                << "\"verificationTrials\"; defaulting to 1 trial.\n";
        }
    }
    
    vec_t<F> trialInputs;
    vec_t<F> trueValues;
    vec_t<F> outputValues;
    bool valid = verifyProofPolynomial(
        proof, fg, cutset, etas, tau, nTrials,
        trialInputs, trueValues, outputValues
    );
    
    str_t resultFile = "verification_result_" + timestamp();
    FileHandler hOut(resultFile);
    writeParamToFile(hOut, "proofIsValid", valid);
    writeParamToFile(hOut, "trialInputs", trialInputs);
    writeParamToFile(hOut, "outputValues", outputValues);
    writeParamToFile(hOut, "trueValues", trueValues);
    
    return resultFile;
}

template <typename F>
str_t contractMod(str_t fgFile, Options options=Options())
{
    str_t evalsFile = computeEvaluations<F>(fgFile);
    str_t coeffsFile = decodeProofCoefficients<F>(evalsFile);
    str_t verificationResultFile = verifyProof<F>(fgFile, coeffsFile);
    FileHandler hv(verificationResultFile);
    bool proofIsValid {false};
    if ( ! loadParamFromFile(hv, "proofIsValid", proofIsValid)) {
        std::cout << "Error: Could not determine validity of proof in file "
                  << coeffsFile << " based on verification results in file "
                  << verificationResultFile << ".\n";
        return "";
    }
    if ( ! proofIsValid) {
        std::cout << "Failure: The proof in file " << coeffsFile
                  << " is invalid!\n";
        return "";
    }
    
    FileHandler hc(coeffsFile);
    vec_t<F> coeffs;
    assert ( loadParamFromFile(hc, "coefficients", coeffs) );
    Poly<F> proofPoly(coeffs);
    
    FileHandler hfg(fgFile);
    vec_t<Var> boundary;
    assert ( loadBoundary(hfg, boundary) );
    vec_t<Var> cutset;
    assert ( loadParamFromFile(hfg, "cutset", cutset) );
    std::sort(cutset.begin(), cutset.end());
    map_t<vec_t<uint_t>, F> tau;
    loadOrConstructTau(hfg, fgFile, cutset, tau);
    
    Factor<F> modularResult = boundaryMapFromProof(
        proofPoly, cutset, boundary, tau);
    Factor<uint_t> result = modToUintFactor<F>(modularResult);
    
    str_t resultFile = "contracted_factor_" + timestamp();
    FileHandler hOut(resultFile);
    writeParamToFile(hOut, "type", str_t("uint_t"));
    uint_t modulus = F::characteristic;
    writeParamToFile(hOut, "modulus", modulus);
    writeFactorToFile<uint_t>(hOut, result);
    
    if ( ! options.keepIntermediateFiles) {
        remove(evalsFile.c_str());
        remove(coeffsFile.c_str());
        remove(verificationResultFile.c_str());
    }
    
    return resultFile;
}


template <typename P>
str_t demuxCmd(str_t cmd, vec_t<str_t> workFiles, Options opts=Options())
{
    assert ( workFiles.size() >= 1 );
    if (cmd == "contract")
        return contractMod<Zp<P>>(workFiles[0], opts);
    else if (cmd == "evaluations")
        return computeEvaluations<Zp<P>>(workFiles[0]);
    else if (cmd == "generate")
        return generateFG<Zp<P>>(workFiles[0]);
    else if (cmd == "reduce")
        return crtReduce<P>(workFiles[0]);
    else if (cmd == "decode")
        return decodeProofCoefficients<Zp<P>>(workFiles[0]);
    else if (cmd == "verify") {
        assert ( workFiles.size() == 2 );  // factorgraph, proof
        return verifyProof<Zp<P>>(workFiles[0], workFiles[1]);
    }
    else {
        std::cout << "Error: Unrecognized command \"" << cmd << "\"\n";
        return "";
    }
}


#include "demux_mod.cpp"


str_t generate(str_t paramsFile)
{
    FileHandler hIn(paramsFile);
    str_t type;
    if ( ! loadParamFromFile<str_t>(hIn, "type", type) ) {
        std::cout << "Error: Could not find type parameter in file "
                  << paramsFile <<".\n";
        return "";
    }
    if (type == "float") {
        str_t floatFGFile = generateFG<float>(paramsFile);
        FileHandler h(floatFGFile);
        writeParamToFile(h, "type", str_t("float"));
        return floatFGFile;
    }
    else if (type == "MPInt") {
        str_t mpIntFGFile = generateFG<MPInt>(paramsFile);
        FileHandler h(mpIntFGFile);
        writeParamToFile(h, "type", str_t("MPInt"));
        return mpIntFGFile;
    }
    else if (type == "modular") {
        str_t modFGFile = demuxMod({paramsFile}, "generate");
        FileHandler hOut(modFGFile);
        writeParamToFile(hOut, "type", str_t("modular"));
        copyParam<uint_t>(hIn, hOut, str_t("modulus"));
        return modFGFile;
    }
    else {
        std::cout << "Error: Unsupported FactorGraph type: "<< type << ".\n";
        return "";
    }
}

void copyStandardParams(str_t srcFile, str_t dstFile)
{
    using F = long_uint_t;
    FileHandler srcH(srcFile);
    FileHandler dstH(dstFile);
    copyParam< vec_t<F> >(srcH, dstH, "evalPoints");
    copyParam< uint_t >(srcH, dstH, "redundantEvaluations");
    copyParam< uint_t >(srcH, dstH, "verificationTrials");
    copyParam< vec_t<Var> >(srcH, dstH, "cutset");
    copyParam< int >(srcH, dstH, "conversionExp");
    copyParam< map_t<label_t, int> >(srcH, dstH, "conversionExps");
    copyParam< map_t<label_t, int> >(srcH, dstH, "upperBoundExps");
    copyParam< MPInt >(srcH, dstH, "finalMapMaxAbsValue");
    copyParam< map_t<Var,map_t<uint_t,F>> >(srcH, dstH, "etas");
    copyParam< map_t<vec_t<uint_t>,F> >(srcH, dstH, "tau");
}

str_t floatToMPIntFG(str_t floatFGFile, bool copyParams=true)
{
    FileHandler floatFGFileH(floatFGFile);
    FactorGraph<float> floatFG;
    floatFGFileH.toBeginning();
    assert ( loadFactorGraph<float>(floatFGFileH, floatFG) );
    
    map_t<label_t, int> conversionExps;
    map_t<label_t, int> upperBoundExps;
    FactorGraph<MPInt> mpIntFG = floatToMPIntFG(
        floatFG, conversionExps, upperBoundExps);
    
    // Compute upper bound for values of the fully contracted mpIntFG
    long expSum{0};
    for (auto labelExpPair : conversionExps)
        expSum += labelExpPair.second;
    for (auto labelExpPair : upperBoundExps)
        expSum += labelExpPair.second;
    assert ( expSum >= 0 );
    // Each value of the contracted graph should fit into 2^32 bits.
    // NOTE: If a more performant version of this software is written,
    // the below assertion may need to be changed.
    assert ( expSum < 1L<<32 );
    mpz_t ub; mpz_init_set_si(ub, 1);
    mpz_mul_2exp(ub, ub, (uint_t)expSum);
    for (Var v : mpIntFG.getVars())
        mpz_mul_ui(ub, ub, v.size);
    MPInt upperBound(ub);
    mpz_clear(ub);
    
    str_t mpIntFGFile = "mpInt_factorgraph_" + timestamp();
    FileHandler mpIntFGFileH(mpIntFGFile);
    writeParamToFile(mpIntFGFileH, "type", str_t("MPInt"));
    writeParamToFile(mpIntFGFileH, "conversionExp", sumValues(conversionExps));
    writeParamToFile(mpIntFGFileH, "conversionExps", conversionExps);
    writeParamToFile(mpIntFGFileH, "upperBoundExps", upperBoundExps);
    writeParamToFile(mpIntFGFileH, "finalMapMaxAbsValue", upperBound);
    writeFGToFile<MPInt>(mpIntFGFileH, mpIntFG);
    mpIntFGFileH.close();
    if (copyParams)
        copyStandardParams(floatFGFile, mpIntFGFile);
    
    return mpIntFGFile;
}

vec_t<str_t> splitMPIntFG(str_t mpIntFGFile, bool copyParams=true)
{
    FileHandler fh(mpIntFGFile);
    FactorGraph<MPInt> g;
    assert ( loadFactorGraph<MPInt>(fh, g) );
    MPInt M;
    assert ( loadParamFromFile(fh, "finalMapMaxAbsValue", M) );
    vec_t<uint_t> moduli = primesToReach(2*M);
    vec_t<str_t> crtReducedFileNames;
    for (uint_t modulus : moduli) {
        str_t fgFile = mpIntFGFile + "_mod_" + std::to_string(modulus);
        FileHandler fhm(fgFile);
        writeParamToFile(fhm, "modulus", modulus);
        writeFGToFile(fhm, g);
        str_t crtReducedFile = demuxMod({fgFile}, "reduce");
        crtReducedFileNames.push_back(crtReducedFile);
        fhm.close();
        remove(fgFile.c_str());
        if (copyParams)
            copyStandardParams(mpIntFGFile, crtReducedFile);
    }
    return crtReducedFileNames;
}

str_t crtMerge(vec_t<str_t> factorFiles)
{
    assert ( factorFiles.size() > 0 );
    vec_t<Factor<uint_t>> gs;
    vec_t<uint_t> moduli;
    bool recoverNegatives {true};
    for (str_t fgFile : factorFiles) {
        FileHandler h(fgFile);
        
        str_t type;
        assert ( loadParamFromFile(h, "type", type) );
        assert ( ! strcmp(type.c_str(), "uint_t") );
        
        uint_t modulus;
        assert ( loadParamFromFile(h, "modulus", modulus) );
        moduli.push_back(modulus);
        
        h.toBeginning();
        Factor<uint_t> g;
        assert ( loadObject(h, g) );
        gs.push_back(g);
        
        loadParamFromFile(h, "recoverNegatives", recoverNegatives);
    }
    
    Factor<MPInt> mergedG = crtReconstruct(gs, moduli, 0, recoverNegatives);
    str_t resultFile = "merged_mpInt_factor_" + timestamp();
    FileHandler hOut(resultFile);
    writeFactorToFile(hOut, mergedG);
    writeParamToFile(hOut, "type", str_t("MPInt"));
    
    return resultFile;
}

str_t contract(str_t inputFile, Options options=Options())
{
    FileHandler hIn(inputFile);
    str_t type;
    if ( ! loadParamFromFile<str_t>(hIn, "type", type) ) {
        std::cout << "Error: Could not find type parameter in file "
                  << inputFile <<".\n";
        return "";
    }
    vec_t<Var> cutset;
    if ( ! loadParamFromFile(hIn, "cutset", cutset)) {
        std::cout <<"Error: could not find cutset in file "<< inputFile <<"\n";
        return "";
    }
    
    if (type == "float") {
        str_t mpIntFGFile = floatToMPIntFG(inputFile);
        copyStandardParams(inputFile, mpIntFGFile);
        str_t mpIntFactorFile = contract(mpIntFGFile, options);
        
        FileHandler hMPInt(mpIntFactorFile);
        Factor<MPInt> mpIntFactor;
        assert ( loadObject(hMPInt, mpIntFactor) );
        int conversionExp;
        assert ( loadParamFromFile(hMPInt, "conversionExp", conversionExp) );
        bool allowSubnormOrInf {false};
        loadParamFromFile(hMPInt, "allowSubnormOrInf", allowSubnormOrInf);
        Factor<float> floatFactor = mpIntToFloatFactor(
            mpIntFactor, -conversionExp, allowSubnormOrInf);
        
        str_t resultFile = "merged_float_factor_" + timestamp();
        FileHandler hOut(resultFile);
        writeParamToFile(hOut, "type", str_t("float"));
        writeFactorToFile(hOut, floatFactor);
        
        if ( ! options.keepIntermediateFiles) {
            remove(mpIntFGFile.c_str());
            hMPInt.close();
            remove(mpIntFactorFile.c_str());
        }
        
        return resultFile;
    }
    
    else if (type == "MPInt") {
        vec_t<str_t> crtReducedFGFiles = splitMPIntFG(inputFile);
        vec_t<str_t> contractedFGFiles;
        for (str_t crtReducedFGFile : crtReducedFGFiles) {
            copyStandardParams(inputFile, crtReducedFGFile);
            contractedFGFiles.push_back(
                demuxMod({crtReducedFGFile}, "contract", options));
        }
        str_t mergedMPIntFGFile = crtMerge(contractedFGFiles);
        
        int conversionExp;
        if (loadParamFromFile(hIn, "conversionExp", conversionExp)) {
            FileHandler hOut(mergedMPIntFGFile);
            writeParamToFile(hOut, "conversionExp", conversionExp);
        }
        else {
            std::cout << "Warning: could not find parameter conversionExp "
                      << "in file " << inputFile << ".\n";
        }
        
        if ( ! options.keepIntermediateFiles) {
            for (str_t filename : crtReducedFGFiles)
                remove(filename.c_str());
            for (str_t filename : contractedFGFiles)
                remove(filename.c_str());
        }
        
        return mergedMPIntFGFile;
    }
    
    else if (type == "modular")
        return demuxMod({inputFile}, "contract", options);
    
    else {
        std::cout << "Error: Unsupported FactorGraph type: "<< type << ".\n";
        return "";
    }
}

vec_t<str_t> split(str_t inputFile, Options options=Options())
{
    FileHandler fh(inputFile);
    str_t type;
    if ( ! loadParamFromFile<str_t>(fh, "type", type) ) {
        std::cout << "Error: Could not find type parameter in file "
                  << inputFile <<".\n";
        return {""};
    }
    if (type == "modular") {
        std::cout << "FactorGraph already in modular form; doing nothing.\n";
        return vec_t<str_t> {inputFile};
    }
    else if (type == "MPInt")
        return splitMPIntFG(inputFile, true);
    else if (type == "float") {
        str_t mpIntFGFile = floatToMPIntFG(inputFile);
        vec_t<str_t> resultFiles = splitMPIntFG(mpIntFGFile, true);
        if ( ! options.keepIntermediateFiles)
            remove(mpIntFGFile.c_str());
        return resultFiles;
    }
    else {
        std::cout << "Error: Unsupported type \"" << type << "\"\n";
        return {""};
    }
}


int checkAndReportResults(vec_t<str_t> resultFiles, Options opts)
{
    if ( (resultFiles.size() == 0) || (resultFiles[0] == "") ) {
        std::cout << "An error occurred; some intermediate may have been "
                  << "written to file; final results not available.\n";
        return 1;
    }
    uint_t nOutFiles = opts.outputFilenames.size();
    uint_t nResultFiles = resultFiles.size();
    if (nOutFiles > 0) {
        if (nOutFiles < nResultFiles)
            std::cout << "Warning: Insufficient output files specified; "
                      << "some files will have default names.\n";
        if (nOutFiles > nResultFiles)
            std::cout << "Warning: Too many output files specified; "
                      << "some output file names will be ignored.\n";
        for (uint_t i = 0; i < nResultFiles && i < nOutFiles; ++i) {
            rename(resultFiles[i].c_str(), opts.outputFilenames[i].c_str());
            resultFiles[i] = opts.outputFilenames[i];
        }
    }
    std::cout << "Results written to file(s):\n";
    for (str_t filename : resultFiles)
        std::cout << "\t" << filename << "\n";
    return 0;
}

