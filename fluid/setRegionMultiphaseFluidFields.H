// Binds region i's objects to the names used by the equation files of
// reactingMultiphaseEulerFoam (EEqns.H, YEqns.H, pU/, pUf/).

multiphaseSystem& fluid = phaseSystemFluid[i];

multiphaseSystem::phaseModelList& phases = fluid.phases();

forAll(phases, phasei)
{
    phases[phasei].thermoRef().validate(args.executable(), "h", "e");
}

surfaceScalarField& phi = fluid.phi();

const IOMRFZoneList& MRF = fluid.MRF();
fv::options& fvOptions = fluid.fvOptions();

volScalarField& p = phases[0].thermoRef().p();

volScalarField& p_rgh = p_rghFluid[i];

const volScalarField& gh = ghFluid[i];
const surfaceScalarField& ghf = ghfFluid[i];

bool frozenFlow = frozenFlowFluid[i];

pimpleControl& pimple = pimpleFluid[i];

const dimensionedScalar& pMin = pMinFluid[i];
