/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
    Modified 2026 E. de Florinier (University of Nottingham): multiphase
    version of chtMultiRegionTwoPhaseEulerFoam.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, licensed under the GNU General Public
    License v3 or later. See <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "phaseSystem.H"
#include "mappedPatchBase.H"
#include "solidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::compressible::
    turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::regionType
>
Foam::compressible::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::regionTypeNames_
{
    { regionType::solid, "solid" },
    { regionType::fluid, "fluid" },
};


const Foam::Enum
<
    Foam::compressible::
    turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::KMethodType
>
Foam::compressible::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::KMethodTypeNames_
{
    { KMethodType::mtSolidThermo, "solidThermo" },
    { KMethodType::mtLookup, "lookup" },
    { KMethodType::mtPhaseSystem, "phaseSystem" }
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::compressible::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::includesPhase
(
    const word& phaseName
) const
{
    if (phaseNames_.empty())
    {
        return true;
    }

    forAll(phaseNames_, namei)
    {
        if (phaseNames_[namei] == phaseName)
        {
            return true;
        }
    }

    return false;
}


Foam::tmp<Foam::scalarField> Foam::compressible::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::kappa
(
    const scalarField& Tp
) const
{
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const label patchi = patch().index();

    switch (method_)
    {
        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            return thermo.kappa(patchi);
            break;
        }

        case mtLookup:
        {
            {
                const auto* ptr =
                    mesh.cfindObject<volScalarField>(kappaName_);

                if (ptr)
                {
                    return patch().patchField(*ptr);
                }
            }
            {
                const auto* ptr =
                    mesh.cfindObject<volSymmTensorField>(kappaName_);

                if (ptr)
                {
                    const symmTensorField& KWall = patch().patchField(*ptr);

                    const vectorField n(patch().nf());

                    return n & KWall & n;
                }
            }

            FatalErrorInFunction
                << "Did not find field " << kappaName_
                << " on mesh " << mesh.name() << " patch " << patch().name()
                << nl
                << "    Please set 'kappa' to the name of a volScalarField"
                << " or volSymmTensorField."
                << exit(FatalError);

            break;
        }

        case mtPhaseSystem:
        {
            const phaseSystem& fluid =
                mesh.lookupObject<phaseSystem>("phaseProperties");

            auto tkappaEff = tmp<scalarField>::New(patch().size(), Zero);
            auto& kappaEff = tkappaEff.ref();

            forAll(fluid.phases(), phasei)
            {
                const phaseModel& phase = fluid.phases()[phasei];

                if (!includesPhase(phase.name()))
                {
                    continue;
                }

                const fvPatchScalarField& alpha = phase.boundaryField()[patchi];

                kappaEff += alpha*phase.kappaEff(patchi)();
            }

            return tkappaEff;

            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "Please set 'kappaMethod' to one of "
                << flatOutput(KMethodTypeNames_.sortedToc()) << nl
                << "and 'kappa' to the name of the volScalar"
                << exit(FatalError);
        }
    }

    // Return zero-sized (not nullptr)
    return tmp<scalarField>::New();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    regionType_(fluid),
    method_(mtLookup),
    kappaName_("none"),
    phaseNames_(),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
(
    const this_bctype& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(psf, p, iF, mapper),
    regionType_(psf.regionType_),
    method_(psf.method_),
    kappaName_(psf.kappaName_),
    phaseNames_(psf.phaseNames_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_)
{}


turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    regionType_(regionTypeNames_.get("region", dict)),
    method_(KMethodTypeNames_.get("kappaMethod", dict)),
    kappaName_(dict.getOrDefault<word>("kappa", "none")),
    phaseNames_(dict.getOrDefault<wordList>("phases", wordList())),
    TnbrName_(dict.getOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.getOrDefault<word>("qrNbr", "none")),
    qrName_(dict.getOrDefault<word>("qr", "none"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.found("otherPhase"))
    {
        WarningInFunction
            << "Entry 'otherPhase' is ignored by "
            << typeName << ": all phases (or those listed in 'phases')"
            << " take part in the wall balance." << nl
            << "    patch " << p.name()
            << " of field " << internalField().name() << endl;
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::
turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
(
    const this_bctype& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(psf, iF),
    regionType_(psf.regionType_),
    method_(psf.method_),
    kappaName_(psf.kappaName_),
    phaseNames_(psf.phaseNames_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // Get the coupling information from the mappedPatchBase
    const label patchi = patch().index();
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    scalarField& Tp = *this;

    const auto& nbrField =
        refCast<const this_bctype>
        (
            nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
        );

    scalarField qr(Tp.size(), 0.0);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField>(qrName_);
    }

    scalarField qrNbr(Tp.size(), 0.0);
    if (qrNbrName_ != "none")
    {
        qrNbr = nbrPatch.lookupPatchField<volScalarField>(qrNbrName_);
        mpp.distribute(qrNbr);
    }

    if (regionType_ == solid)
    {
        // Neighbour is the fluid: every included phase adds its conductance
        // K_k = alpha_k kappaEff_k deltaCoeffs at its near-wall temperature.
        const phaseSystem& fluid =
            nbrMesh.lookupObject<phaseSystem>("phaseProperties");

        scalarField KSum(Tp.size(), Zero);
        scalarField c(Tp.size(), Zero);

        forAll(fluid.phases(), phasei)
        {
            const phaseModel& phase = fluid.phases()[phasei];

            if (!nbrField.includesPhase(phase.name()))
            {
                continue;
            }

            const fvPatchScalarField& alphak =
                phase.boundaryField()[samplePatchi];

            scalarField Kk
            (
                alphak*phase.kappaEff(samplePatchi)*nbrPatch.deltaCoeffs()
            );
            mpp.distribute(Kk);

            scalarField Tk
            (
                phase.thermo().T().boundaryField()[samplePatchi]
               .patchInternalField()
            );
            mpp.distribute(Tk);

            KSum += Kk;
            c += Kk*Tk;
        }

        const scalarField KDelta(kappa(Tp)*patch().deltaCoeffs());

        valueFraction() = KSum/(KSum + KDelta);
        refValue() = c/max(KSum, VSMALL);
        refGrad() = (qr + qrNbr)/kappa(Tp);

        if (debug)
        {
            const scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

            Info<< "T solid (" << patch().name() << "):" << nl
                << "    heat transfer rate from solid: " << Q << nl
                << "    wall temperature min/avg/max: "
                << gMin(Tp) << " / " << gAverage(Tp) << " / " << gMax(Tp)
                << nl << endl;
        }
    }
    else if (regionType_ == fluid)
    {
        // This condition sits on the temperature of one phase ("self").
        // The solid and every other included phase form the far side of
        // the balance; self is the near side.
        const phaseSystem& fluid =
            mesh.lookupObject<phaseSystem>("phaseProperties");

        const word selfName(internalField().group());

        scalarField TcNbr(nbrField.patchInternalField());
        mpp.distribute(TcNbr);

        scalarField KDeltaNbr(nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs());
        mpp.distribute(KDeltaNbr);

        scalarField a(KDeltaNbr);
        scalarField c(KDeltaNbr*TcNbr);
        scalarField KSelf(Tp.size(), Zero);

        forAll(fluid.phases(), phasei)
        {
            const phaseModel& phase = fluid.phases()[phasei];

            const bool isSelf = (phase.name() == selfName);

            if (!isSelf && !includesPhase(phase.name()))
            {
                continue;
            }

            const fvPatchScalarField& alphak = phase.boundaryField()[patchi];

            const scalarField Kk
            (
                alphak*phase.kappaEff(patchi)*patch().deltaCoeffs()
            );

            if (isSelf)
            {
                KSelf = Kk;
            }
            else
            {
                const scalarField Tk
                (
                    phase.thermo().T().boundaryField()[patchi]
                   .patchInternalField()
                );

                a += Kk;
                c += Kk*Tk;
            }
        }

        valueFraction() = a/(a + KSelf);
        refValue() = c/a;
        refGrad() = (qr + qrNbr)/kappa(Tp);

        if (debug)
        {
            Info<< "T fluid (" << patch().name() << ", " << selfName << "):"
                << nl;

            scalar QTotal = 0;

            forAll(fluid.phases(), phasei)
            {
                const phaseModel& phase = fluid.phases()[phasei];

                if
                (
                    phase.name() != selfName
                 && !includesPhase(phase.name())
                )
                {
                    continue;
                }

                const fvPatchScalarField& alphak =
                    phase.boundaryField()[patchi];

                const scalarField Kk
                (
                    alphak*phase.kappaEff(patchi)*patch().deltaCoeffs()
                );

                const scalarField Tk
                (
                    phase.thermo().T().boundaryField()[patchi]
                   .patchInternalField()
                );

                const scalar Qk = gSum((Tp - Tk)*Kk*patch().magSf());
                QTotal += Qk;

                Info<< "    heat transfer rate to " << phase.name()
                    << ": " << Qk << nl;
            }

            Info<< "    heat transfer rate to fluid (sum): " << QTotal << nl
                << "    wall temperature min/avg/max: "
                << gMin(Tp) << " / " << gAverage(Tp) << " / " << gMax(Tp)
                << nl << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown region type. Valid types are: "
            << regionTypeNames_ << nl << exit(FatalError);
    }

    UPstream::msgType(oldTag);  // Restore tag

    this->parent_bctype::updateCoeffs();
}


void turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntry("kappaMethod", KMethodTypeNames_[method_]);
    os.writeEntryIfDifferent<word>("kappa", "none", kappaName_);

    os.writeEntry("Tnbr", TnbrName_);

    os.writeEntryIfDifferent<word>("qrNbr", "none", qrNbrName_);
    os.writeEntryIfDifferent<word>("qr", "none", qrName_);

    os.writeEntry("region", regionTypeNames_[regionType_]);

    if (!phaseNames_.empty())
    {
        os.writeEntry("phases", phaseNames_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
