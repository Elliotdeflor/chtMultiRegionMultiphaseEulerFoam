// Copy of chtMultiRegionTwoPhaseEulerFoam/solid/solveSolid.H (segregated solid
// enthalpy solve, no coupled matrix assembly). Renamed so it cannot be
// shadowed by chtMultiRegionFoam/solid/solveSolid.H.

if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rho, h)
          - thermo.heatDiffusion(betav, h)
          ==
            fvOptions(rho, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve(h.select(finalIter));

        fvOptions.correct(h);
    }

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
