import subprocess
import sys
import os
import time
import json
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
plt.rcParams.update({"font.size": 18, "figure.autolayout": True})


def timeTag():
    return str(time.time())


def runExperiment(
    executableName, experimentName, solutionAlgos, instanceParams, outputFile,
    numRuns=10, verbosity=1
    ):
    
    commandParts = [
        "./" + executableName,
        experimentName,
        "--solution-algos", *solutionAlgos,
        "--time-log", outputFile,
        "--verbose", str(verbosity),
        "--number-of-runs", str(numRuns)
        ]
    for ps in instanceParams:
        commandParts.append("--instance-params")
        for p in ps:
            commandParts.append(str(p))
    commandString = ' '.join(commandParts)
    print("About to run command\n\t", commandString)
    #proceed = input("Proceed? [y/n] ")
    #if proceed != 'y':
    #    print("Aborting.\n"); sys.exit()
    startTime = time.time()
    subprocess.run(commandParts)
    endTime = time.time()
    print("Experiment completed in ", endTime - startTime, " seconds.\n")


def catenateWith(words, interstice="_"):
    result = ""
    for i in range(len(words)):
        result += str(words[i])
        if i < len(words) - 1:
            result += interstice
    return result


# Define parameters for the experiments

config = {"matMul": {}, "permanent": {}}

# matMul
numMatrices = 2
config["matMul"]["ns"] = [n for n in range(4,7)]
config["matMul"]["instanceParams"] = {
    n: [numMatrices, numMatrices, n, n, 0, 999]
    for n in config["matMul"]["ns"]
    }
config["matMul"]["solutionAlgos"] = ["dfgc"]
config["matMul"]["numRuns"] = 5

# matrix permanent
config["permanent"]["ns"] = [n for n in range(4,7)]
config["permanent"]["instanceParams"] = {
    n: [n, n] for n in config["permanent"]["ns"]
    }
config["permanent"]["solutionAlgos"] = ["dfgc"]
config["permanent"]["numRuns"] = 5


def experimentIdToTitle(experimentId):
    mapping = {
        "matMul": "Matrix multiplication",
        "permanent": "Matrix permanent"
        }
    return mapping[experimentId]

def experimentIdToXLabel(experimentId):
    mapping = {
        "matMul": "number of rows in matrix (n)",
        "permanent": "number of rows in matrix (n)",
        }
    return mapping[experimentId]

def fullAlgoName(experimentId, algoId):
    middle = "MatMul" if experimentId == "matMul" else "Permanent"
    endMap = {
        "dfgc": "DFGC",
        "base": "Base",
        "dfgc_slow": "DFGC_slow"
    }
    end = endMap[algoId]
    return "solve" + middle + end


def findDataFile(fnameStub, searchDirs):
    k = len(fnameStub)
    for sd in searchDirs:
        try:
            fileNames = os.listdir(sd)
            for fn in fileNames:
                if fn[0:k] == fnameStub:
                    return sd + "/" + fn
        except FileNotFoundError as err:
            print("Error while trying to search in directory", sd, ": ", err)
            print("skipping", sd)
    return ""

def parseDataFromFile(filename, subroutines, solutionAlgoName=None):
    # TODO: Less brittle implementation
    d = {}
    with open(filename) as f:
        contents = f.read()
        lines = contents.split("\n")
        
        # Extract solution summary
        ip = 0  # start index of summary for proof-construction
        key = "Summary:"
        while ip < len(lines) and lines[ip][0:len(key)] != key:
            ip += 1
        iv = ip + 1  # start index of summary for proof-verification
        while iv < len(lines) and lines[iv][0:len(key)] != key:
            iv += 1
        summaryStr = "\n".join(lines[ip+1 : iv-1])
        summary = json.loads(summaryStr)
        
        # Stats of the solution algorithm as a whole (without verification)
        if not (solutionAlgoName is None):
            solutionSummary = summary[solutionAlgoName]
            d["solution"] = {
                "total": solutionSummary["time"],
                "countsByParams": solutionSummary["msg_counts"],
                "timesByParams": solutionSummary["times_by_msg"]}
        
        # Stats for proof verification
        if (iv < len(lines)):  # verification summary was found (maybe)
            verificationSummaryStr = "\n".join(lines[iv+1:])
            verificationSummary = json.loads(verificationSummaryStr)
            summary["verification"] = verificationSummary["verifyProofPolynomial"]
            
        # Parse runtimes and call counts of subroutines
        for s in subroutines:
            if not (s in summary.keys()):
                continue
            d[s] = {
                "total": summary[s]["time"],
                "countsByParams": summary[s]["msg_counts"],
                "timesByParams": summary[s]["times_by_msg"]}
    
    return d
    

def subroutineToLabel(subroutineName):
    if subroutineName[0:8] == "solution":
        return "prover"
    mapping = {
        "rsDecodeXY": "decode",
        "arrayPolyInterpolate": "interpolate",
        "contractFactors": "contract",
        "evaluateProofPolynomial:permuteToEnd": "permute",
        "verification": "verifier"
        }
    return mapping.get(subroutineName, subroutineName)


def averageOverRuns(dicts):
    # TODO: less fragile implementation
    avg = dicts[0]
    N = len(dicts)
    for i in range(1, N):
        avg["total"] += dicts[i]["total"]
        counts = dicts[i]["countsByParams"]
        for p in counts.keys():
            avg["countsByParams"][p] += counts[p]
        times = dicts[i]["timesByParams"]
        for p in times.keys():
            for t in range(len(times[p])):
                avg["timesByParams"][p][t] += times[p][t]
    
    avg["total"] /= float(N)
    for p in avg["countsByParams"].keys():
        avg["countsByParams"][p] /= float(N)
    for p in avg["timesByParams"].keys():
        for t in range(len(avg["timesByParams"][p])):
            avg["timesByParams"][p][t] /= float(N)
    
    return avg


def firstOfRuns(dicts):
    return dicts[0]


def summarizeOverRuns(dicts, summary="avg"):
    if summary == "avg":
        return averageOverRuns(dicts)
    elif summary == "first":
        return firstOfRuns(dicts)
    else:
        raise NotImplementedError


def subroutineTableStr(subroutine, n, data, writeN=False):
    #TODO: Less fragile implementation
    
    counts = data["countsByParams"]
    times = data["timesByParams"]
    strToFloat = lambda s: 0 if s=="" else float(s)
    params = sorted(list(counts.keys()), key=strToFloat)
    assert params == sorted(list(times.keys()), key=strToFloat)
    assert len(params) >= 1
    
    # First line differs from later lines; construct separately
    p0 = params[0]
    nStr = str(n) if writeN else " "
    result = catenateWith(
        [nStr, subroutine, p0, int(counts[p0]),
         "{:.4f}".format(times[p0]) + " s"],
        " & ")
    
    for i in range(1, len(params)):
        result += " \\\\ \n"
        p = params[i]
        result += catenateWith(
            [" ", " ", p, int(counts[p]),
             "{:.4f}".format(times[p]) + " s"],
            " & ")
    
    return result


def plotRoutines(data, experimentName, routines, markers, colors,
                 fname="", scales=[], alpha=1.0):
    ns = np.asarray(list(data.keys()))
    if not scales:
        scales = ["linlin"]
    nSubplots = len(scales)
    fig, axes = plt.subplots(nSubplots, 1, figsize=(7,5*nSubplots))
    if not hasattr(axes, "__len__"):
        axes = [axes]
    for j in range(len(ns)):
        n = ns[j]
        runsDict = data[n]
        runIds = list(runsDict.keys())
        xs = np.repeat(n, len(runIds))
        for i in range(len(routines)):
            s = routines[i]
            if not s in runsDict[0].keys():
                continue
            ys = np.asarray(
                [runsDict[r][s]["total"] for r in runIds])
            for ax in axes:
                if j == 0:
                    ax.plot(xs, ys, markers[i], color=colors[i],
                            linestyle='', alpha=alpha,
                            label=subroutineToLabel(s))
                else:
                    ax.plot(xs, ys, markers[i], color=colors[i], 
                            linestyle='', alpha=alpha)
    axes[-1].set_xlabel(experimentIdToXLabel(experimentName))
    for axInd in range(len(axes)):
        ax = axes[axInd]
        ax.set_ylabel("wall time [s]")
        ax.set_xticks(ns)
        ax.legend(loc="upper left", frameon=True)
        if len(ns) > 20 and (scales[axInd][0:3] == "lin"):
            tickVisibility = lambda t: t%2==0
            for (t, tick) in enumerate(ax.xaxis.get_ticklabels()):
                tick.set_visible(tickVisibility(t))
        if scales[axInd] == "linlin":
            continue
        elif scales[axInd] == "linlog":
            ax.set_yscale("log", base=2)
        elif scales[axInd] == "loglog":
            ax.set_xscale("log", base=2)
            ax.set_yscale("log", base=2)
    if nSubplots == 1:
        plt.title(experimentIdToTitle(experimentName))
    else:
        fig.suptitle(experimentIdToTitle(experimentName))
    if fname == "":
        fname = "runtimes_" + timeTag() + ".png"
    plt.savefig(fname, format="png")


def main():
    
    workdir = "experiments_" + timeTag()
    try:
        os.mkdir(workdir)
    except FileExistsError:
        pass
    # where to search for existing results:
    searchDirs = ["results_data"]
    subroutines = ["arrayPolyInterpolate", "contractFactors", "rsDecodeXY",
                   "evaluateProofPolynomial:permuteToEnd"]
    executables = ["tbm1", "tbm3"] #compiled with different time-logging levels
    subroutinesMap = {
        "tbm1" : ["rsDecodeXY", "verification"],
        "tbm3" : ["arrayPolyInterpolate", "contractFactors",
                  "evaluateProofPolynomial:permuteToEnd"]
        }
    experimentNames = config.keys()
    
    # Run the experiments, gather data
    data = {}
    for e in experimentNames:
        data[e] = {}
        for a in config[e]["solutionAlgos"]:
            data[e][a] = {}
            for n in config[e]["ns"]:
                data[e][a][n] = {}
                numRuns = config[e].get("numRuns", 5)
                for r in range(numRuns):
                    p = config[e]["instanceParams"][n]
                    dataEANR = {}
                    for i in range(len(executables)):
                        ex = executables[i]
                        fnameStub = catenateWith([e,a,n,r,ex], "_")
                        fname = findDataFile(fnameStub, searchDirs)
                        if fname == "":
                            fname = workdir + "/" + fnameStub + timeTag()
                            runExperiment(ex, e, [a], [p], fname, 1,
                                          config[e].get("verbosity", 1))
                        algoName = fullAlgoName(e,a) if i==0 else None
                        dataEANR.update(parseDataFromFile(
                            fname, subroutinesMap[ex], algoName))
                    data[e][a][n][r] = dataEANR
    
    # Process/distill/display the data
    for e in experimentNames:
        for a in config[e]["solutionAlgos"]:
            
            # Plot total runtimes for prover and verifier
            fname = catenateWith(["times",e,a,"totals",timeTag()]) + ".png"
            scales = ["linlin", "loglog"] if e=="matMul" else ["linlog"]
            plotRoutines(data[e][a], e, ["solution","verification"],
                         ['x','v'], ["black",'c'], fname, scales)
            
            # Plot runtimes for prover subroutines
            fname = catenateWith(["times",e,a,"subroutines",timeTag()]) + ".png"
            scales = ["loglog"] if e=="matMul" else ["linlog"]
            plotRoutines(data[e][a], e, subroutines,
                         ['o','^','s','*','+'], ['b','g','r','m','y'],
                         fname, scales)
            
            # Write more detailed results to LaTeX table
            
            fname = workdir + "/" + e + "_" + a + "_summary_" + timeTag()
            with open(fname, "w") as f:
                
                # Write table "header"
                f.write("\\begin{table}[htb]\n"
                        + "\\begin{center}\n"
                        + "\\caption{" + e + a + "}\n"
                        + "\\label{tbl:" + e + a + "}\n"
                        + "\\begin{tabular}{c|r|r|r|r}\n"
                        + "\\hline \n"
                        + "$n$ & Subroutine & Input size & Count & Time [s]\n")
                
                # Write a table section for each n in nsToWrite
                allRoutines = ["solution"] + subroutines + ["verification"]
                ns = np.asarray(list(data[e][a].keys()))
                nFilter = lambda n: n % 4 == 0
                if e == "matMul":
                    nFilter = lambda n: n in [8,16,32]
                nsToWrite = list(filter(nFilter, ns))
                for i in range(len(nsToWrite)):
                    f.write("\\\\ \n \\hline\n")
                    n = nsToWrite[i]
                    dataN = data[e][a][n]
                    for j in range(len(allRoutines)):
                        s = allRoutines[j]
                        if not (s in dataN[0].keys()):
                            continue
                        allDataNS = [dataN[r][s] for r in dataN.keys()]
                        dataNS = summarizeOverRuns(allDataNS, "first")
                        f.write(subroutineTableStr(
                            subroutineToLabel(s), n, dataNS, j==0))
                        if j < len(allRoutines) - 1:
                            f.write("\\\\ \n\\cline{2-5}")
                
                # Write table "footer"
                f.write("\\\\ \n\hline\n\\end{tabular}\n" 
                        + "\\end{center}\n"
                        + "\\end{table}\n")

    return data


if __name__ == "__main__":
    main()

