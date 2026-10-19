// One outer iteration of the fluid in region i, in the order used by
// reactingMultiphaseEulerFoam. Equation files come from that solver.

if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

if (frozenFlow)
{
    #include "EEqns.H"
}
else
{
    fluid.solve();
    fluid.correct();

    #include "YEqns.H"

    if (faceMomentum)
    {
        #include "pUf/UEqns.H"
        #include "EEqns.H"
        #include "pUf/pEqn.H"
    }
    else
    {
        #include "pU/UEqns.H"
        #include "EEqns.H"
        #include "pU/pEqn.H"
    }

    fluid.correctKinematics();

    // As reactingMultiphaseEulerFoam (pimple.turbCorr()): turbulence is
    // corrected on the final outer iteration only, unless
    // turbOnFinalIterOnly is set to false in the region's PIMPLE dictionary
    if (!turbOnFinalIterOnly || finalIter)
    {
        fluid.correctTurbulence();
    }
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
