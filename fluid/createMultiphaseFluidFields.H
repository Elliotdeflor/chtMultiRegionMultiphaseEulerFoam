// Per-region fluid objects for chtMultiRegionMultiphaseEulerFoam_notGit.
// Replaces chtMultiRegionTwoPhaseEulerFoam/fluid/createFluidFields.H:
// twoPhaseSystem -> multiphaseSystem; unused fields (initialMass,
// pRefCell/pRefValue, fieldsFluid, commented debug fields) removed.

PtrList<multiphaseSystem> phaseSystemFluid(fluidRegions.size());

PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());

PtrList<volScalarField> p_rghFluid(fluidRegions.size());

List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

PtrList<pimpleControl> pimpleFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    pimpleFluid.set
    (
        i,
        new pimpleControl(fluidRegions[i])
    );

    Info<< "    Adding to phaseSystemFluid\n" << endl;
    phaseSystemFluid.set(i, multiphaseSystem::New(fluidRegions[i]).ptr());

    Info<< "    Phases in region " << fluidRegions[i].name() << ": ";
    forAll(phaseSystemFluid[i].phases(), phasei)
    {
        Info<< phaseSystemFluid[i].phases()[phasei].name() << ' ';
    }
    Info<< nl << endl;

    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p and the mixture density
    // (as in reactingMultiphaseEulerFoam/pU/pEqn.H). The two-phase CHT
    // solver used the first phase's density here; p_rgh is reset from the
    // mixture density at the first pressure solve either way.
    p_rghFluid[i] =
        phaseSystemFluid[i].phases()[0].thermo().p()
      - phaseSystemFluid[i].rho()*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");

    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    pMinFluid.set
    (
        i,
        new dimensionedScalar
        (
            "pMin",
            dimPressure,
            phaseSystemFluid[i]
        )
    );
}
