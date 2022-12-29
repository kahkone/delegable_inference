#include "testbench.cpp"

bool isAnOption(str_t s)
{
    return (s.size() > 2)  && (s[0] == '-') && (s[1] == '-');
}

str_t listCat(str_t a, str_t b) { return (a == "") ? b : a + ", " + b; }


int main(int argc, char** argv)
{
    #if TRACK_TIME
    setTimeLogPrecision(9);
    #endif
    using MP = MP2147483629;
    using Z = Zp<MP>;
    using IMatMul = vec_t<Array2D<Z>>;  // type of mat mul problem instances
    using IPerm = Array2D<Z>;  // type of matrix permanent problem instances
    using SMatMul = Array2D<Z>; // type of mat mul problem solutions
    
    vec_t<str_t> experimentNames {"matMul", "permanent"};
    
    str_t helpMessage =
        "Usage:\n"
        + str_t("  experiment-name\n")
        + str_t("   --solution-algos name [name ...]\n")
        + str_t("   --instance-params params\n")
        + str_t("   [--number-of-runs num]\n")
        + str_t("   [--time-log filename]\n")
        + str_t("   [--seed num]\n")
        + str_t("   [--verbose num]\n")
    + str_t("The --instance-params option may be specified multiple times.\n")
    + str_t("Available experiments: ")
        + foldl<str_t>(listCat, "", experimentNames) + str_t("\n");
    
    if (argc < 6) {
        std::cout << helpMessage;
        return 1;
    }
    
    str_t experimentName = argv[1];
    vec_t<str_t> algoNames;
    vec_t<ParamsMatMul<Z>> paramsMatMul;
    vec_t<ParamsPermanent> paramsPermanent;
    uint_t nRuns {10};
    str_t outputFile = "";
    time_t seed = time(NULL);
    int verbosity {0};
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "--solution-algos")) {
            ++i;
            while (i < argc && !isAnOption(argv[i])) {
                algoNames.push_back(argv[i]);
                ++i;
            }
            --i;
        }
        else if (!strcmp(argv[i], "--instance-params")) {
            if (experimentName == "matMul") {
                if (argc < i + 6) {
                    std::cout << "Error: matMul --instance-params requires "
                              << "6 numeric parameters. "
                              << "(See class ParamsMatMul.)\n";
                    return 1;
                }
                ParamsMatMul p(
                    (uint_t)atoi(argv[i+1]), (uint_t)atoi(argv[i+2]),
                    (uint_t)atoi(argv[i+3]), (uint_t)atoi(argv[i+4]),
                    Z((uint_t)atoi(argv[i+5])), Z((uint_t)atoi(argv[i+6]))
                );
                paramsMatMul.push_back(p);
                i += 6;
            }
            else if (experimentName == "permanent") {
                if (argc < i + 2) {
                std::cout << "Error: permanent --instance-params requires"
                          << " 2 numeric parameters."
                          << " (See class ParamsPermanent.)\n";
                return 1;
                }
                ParamsPermanent p(
                    (uint_t)atoi(argv[i+1]), (uint_t)atoi(argv[i+2])
                );
                paramsPermanent.push_back(p);
                i += 2;
            }
        }
        else if (!strcmp(argv[i], "--number-of-runs")) {
            assert ( argc > i+1 );
            int r = atoi(argv[i+1]);
            assert (r > 0);
            nRuns = (uint_t)r;
            ++i;
        }
        else if (!strcmp(argv[i], "--time-log")) {
            assert ( argc > i+1 );
            str_t timeStackOutputFile  = str_t(argv[i+1]);
            ++i;
            #if TRACK_TIME
                setTimeLogFile(timeStackOutputFile);
            #else
                std::cout << "Error: The global time-tracking stack "
                          << "is not available; recompile with "
                          << "\"-D TRACK_TIME=1\"\n";
                return 1;
            #endif
        }
        else if (!strcmp(argv[i], "--seed")) {
            assert ( argc > i+1 );
            seed = atoi(argv[i+1]);
            ++i;
        }
        else if (!strcmp(argv[i], "--verbose")) {
            assert ( argc > i+1 );
            verbosity = atoi(argv[i+1]);
            ++i;
        }
        else {
            std::cout << "Error: Unrecognized option: " << argv[i] << "\n";
            return 1;
        }
    }
    srand(seed);
    std::cout << "Using random seed " << seed << "\n";
    
    TBOptions opts {nRuns, verbosity};
    
    if (experimentName == "matMul") {
        vec_t<
            SMatMul (*)(IMatMul)
        > solutionAlgos;
        for (str_t algoName : algoNames) {
            if (algoName == "base")
                solutionAlgos.push_back( solveMatMulBase<Z> );
            else if (algoName == "dfgc")
                solutionAlgos.push_back( solveMatMulDFGC<Z> );
            else {
                std::cout << "Error: Unrecognized algorithm: "
                          << algoName << "\n";
                return 1;
            }
        }
        vec_t<InstanceGenerator<IMatMul, ParamsMatMul<Z>>> generators;
        for (ParamsMatMul<Z> params : paramsMatMul)
            generators.push_back(
                InstanceGenerator(params, instanceGenMatMul<Z>));
        runExperiments<
            IMatMul, ParamsMatMul<Z>, SMatMul
        >(generators, solutionAlgos, opts);
    }
    
    else if (experimentName == "permanent") {
        vec_t<
            Z (*)(IPerm)
        > solutionAlgos;
        for (str_t algoName : algoNames) {
            if (algoName == "base")
                solutionAlgos.push_back( solvePermanentBase );
            else if (algoName == "dfgc")
                solutionAlgos.push_back( solvePermanentDFGC );
            else if (algoName == "dfgc_slow")
                solutionAlgos.push_back( solvePermanentDFGC_slow );
            else {
                std::cout << "Error: Unrecognized algorithm: "
                          << algoName << "\n";
                return 1;
            }
        }
        vec_t<InstanceGenerator<IPerm, ParamsPermanent>> generators;
        for (ParamsPermanent params : paramsPermanent)
            generators.push_back(
                InstanceGenerator(params, instanceGenPermanent<Z>));
        runExperiments<
            IPerm, ParamsPermanent, Z
        >(generators, solutionAlgos, opts);
    }
    
    else {
        std::cout << "Error: Unrecognized experiment: "<< experimentName <<"\n"
                  << helpMessage;
        return 1;
    }
    
    return 0;
}



