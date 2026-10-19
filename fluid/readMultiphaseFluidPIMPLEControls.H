// Per-region algorithm switches, read from system/<region>/fvSolution PIMPLE

const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");

Switch faceMomentum
(
    pimpleDict.getOrDefault<Switch>("faceMomentum", false)
);

Switch partialElimination
(
    pimpleDict.getOrDefault<Switch>("partialElimination", false)
);

int nEnergyCorrectors
(
    pimpleDict.getOrDefault<int>("nEnergyCorrectors", 1)
);

// Same default as pimpleControl: turbulence on the final outer iteration only
const bool turbOnFinalIterOnly
(
    pimpleDict.getOrDefault<bool>("turbOnFinalIterOnly", true)
);

pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);
