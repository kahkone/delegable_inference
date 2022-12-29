#include <iostream>
#include <cstring>
#include <sys/time.h>
#include "montgomery.cpp"
#include "util.cpp"
#include "fileio.cpp"
#include "factorgraph.cpp"
#include "crt.cpp"
#include "eval_algo.cpp"
#include "poly.cpp"
#include "integration.cpp"

bool isAnOption(str_t s)
{
    return (s.size() > 2)  && (s[0] == '-') && (s[1] == '-');
}

int main(int argc, char* argv[])
{
    str_t helpMessage =
        "Usage:\n"
        + str_t("  generate params-file [options]\n")
        + str_t("  contract factorgraph-file [options]\n")
        + str_t("  split factorgraph-file [options]\n")
        + str_t("  decode evaluations-file [options]\n")
        + str_t("  evaluations factorgraph-file [options]\n")
        + str_t("  verify factorgraph-file proof-file [options]\n")
        + str_t("  merge factorgraph-file [factorgraph-file ...] [options]\n")
        + str_t("Possible options are:\n")
        + str_t("  --keepIntermediateFiles {0|1}\n")
        + str_t("  --output output-file [output-file ...]\n");

    if (argc <= 2 || !strcmp(argv[1],"help") || !strcmp(argv[1],"--help")) {
        std::cout << helpMessage;
        return 0;
    }
    
    str_t cmd(argv[1]);
    Options options {false, {}};
    vec_t<str_t> inputFiles;
    for (int i = 2; i < argc; ++i) {
        str_t arg = str_t(argv[i]);
        if (isAnOption(arg)) {
            assert ( argc > i+1 );
            if ( ! strcmp(argv[i],"--keepIntermediateFiles"))
                options.keepIntermediateFiles = atoi(argv[i+1]);
            if ( ! strcmp(argv[i],"--output")) {
                for (int j = i+1; j < argc; ++j) {
                    if (isAnOption(str_t(argv[j])))
                        break;
                    options.outputFilenames.push_back(str_t(argv[j]));
                }
            }
        }
        else
            inputFiles.push_back(arg);
    }
    
    vec_t<str_t> resultFiles;
    
    if (cmd == "generate")
        resultFiles.push_back( generate(inputFiles[0]) );
    else if (cmd == "contract")
        resultFiles.push_back( contract(inputFiles[0], options) );
    else if (cmd == "split")
        resultFiles = split(inputFiles[0], options);
    else if (cmd == "decode")
        resultFiles.push_back( demuxMod(inputFiles, "decode") );
    else if (cmd == "evaluations")
        resultFiles.push_back( demuxMod(inputFiles, "evaluations") );
    else if (cmd == "verify")
        resultFiles.push_back( demuxMod(inputFiles, "verify") );
    else if (cmd == "merge")
        resultFiles.push_back( crtMerge(inputFiles) );
    else {
        std::cout << "Error: Unrecognized command.\n" << helpMessage;
        return 1;
    }

    int returnStatus = checkAndReportResults(resultFiles, options);
    return returnStatus;
}

