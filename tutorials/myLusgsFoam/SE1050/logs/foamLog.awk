# Awk script for OpenFOAM log file extraction
BEGIN {
    Iteration=0
    resetCounters()
}

# Reset counters used for variable postfix
function resetCounters() {
    executionTimeCnt=0
    kCnt=0
    kAvgCnt=0
    kFinalResCnt=0
    kItersCnt=0
    kMaxCnt=0
    kMinCnt=0
    omegaCnt=0
    omegaFinalResCnt=0
    omegaItersCnt=0
    rhoCnt=0
    rhoECnt=0
    rhoEFinalResCnt=0
    rhoEItersCnt=0
    rhoFinalResCnt=0
    rhoItersCnt=0
    rhoUCnt=0
    rhoUFinalResCnt=0
    rhoUItersCnt=0
    SeparatorCnt=0
    TimeCnt=0
    # Reset counters for 'Solving for ...'
    for (varName in subIter)
    {
        subIter[varName]=0
    }
}

# Extract value after columnSel
function extract(inLine,columnSel,outVar,a,b)
{
    a=index(inLine, columnSel)
    b=length(columnSel)
    split(substr(inLine, a+b),outVar)
    gsub("[,:]","",outVar[1])
}

# Iteration separator (increments 'Iteration')
/^[ \t]*Time = / {
    Iteration++
    resetCounters()
}

# Time extraction (sets 'Time')
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    Time=val[1]
}

# Skip whole line with singularity variable
/solution singularity/ {
    next;
}

# Extract: 'Solving for ...'
/Solving for/ {
    extract($0, "Solving for ", varNameVal)

    varName=varNameVal[1]
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Initial residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "FinalRes"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Final residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "Iters"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "No Iterations ", val)
    print Time "\t" val[1] > file
}

# Extract: 'executionTime'
/ExecutionTime = / {
    extract($0, "ExecutionTime = ", val)
    file="logs/executionTime_" executionTimeCnt
    print Time "\t" val[1] > file
    executionTimeCnt++
}

# Extract: 'kAvg'
/bounding k,/ {
    extract($0, "average: ", val)
    file="logs/kAvg_" kAvgCnt
    print Time "\t" val[1] > file
    kAvgCnt++
}

# Extract: 'kMax'
/bounding k,/ {
    extract($0, "max: ", val)
    file="logs/kMax_" kMaxCnt
    print Time "\t" val[1] > file
    kMaxCnt++
}

# Extract: 'kMin'
/bounding k,/ {
    extract($0, "min: ", val)
    file="logs/kMin_" kMinCnt
    print Time "\t" val[1] > file
    kMinCnt++
}

# Extract: 'Separator'
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="logs/Separator_" SeparatorCnt
    print Time "\t" val[1] > file
    SeparatorCnt++
}

# Extract: 'Time'
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="logs/Time_" TimeCnt
    print Time "\t" val[1] > file
    TimeCnt++
}

# End
