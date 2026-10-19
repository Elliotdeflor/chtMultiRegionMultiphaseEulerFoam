    // No top-level mesh: read the outer-loop controls from system/fvSolution
    // of the case itself (same as chtMultiRegionFoam/readPIMPLEControls.H).
    fvSolution solutionDict(runTime);
    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimple.getOrDefault<int>("nOuterCorrectors", 1);
