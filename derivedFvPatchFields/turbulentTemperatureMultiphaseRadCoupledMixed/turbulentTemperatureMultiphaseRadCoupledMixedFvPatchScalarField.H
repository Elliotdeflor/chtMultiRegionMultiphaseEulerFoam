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

Class
    Foam::compressible::
        turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer
    between a multiphase Euler fluid region and a solid region.

    N-phase generalisation of turbulentTemperatureTwoPhaseRadCoupledMixed.
    Each fluid phase k contributes a wall conductance

        K_k = alpha_k * kappaEff_k * deltaCoeffs

    and the wall temperature is the conductance-weighted balance

        Tw = (K_s T_s + sum_k K_k T_k) / (K_s + sum_k K_k)

    where s is the solid and T_k, T_s are near-wall cell temperatures.
    Wall boiling enters through kappaEff of each phase (alphat set by the
    wall boiling condition), as in the two-phase version.

    On the fluid side the condition may be applied to the temperature of any
    phase, or of every phase: each gives the same wall temperature. The phase
    named by Tnbr on the solid side must carry this condition.

    With two phases and no 'phases' entry this reproduces
    turbulentTemperatureTwoPhaseRadCoupledMixed.

Usage
    \table
        Property    | Description                              | Required | Default
        region      | solid or fluid                           | yes      |
        Tnbr        | temperature field name in the neighbour  | no       | T
        kappaMethod | solidThermo, lookup or phaseSystem       | yes      |
        kappa       | conductivity field name (lookup only)    | no       | none
        phases      | fluid phases in the balance (read on the  | no       | all
                    | fluid side; the solid side uses Tnbr's)  |          |
        qrNbr       | radiative flux name in the neighbour     | no       | none
        qr          | radiative flux name in this region       | no       | none
    \endtable

    Fluid region (e.g. 0/water/T.liquid, optionally also T.gas, ...):
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureMultiphaseRadCoupledMixed;
        region          fluid;
        Tnbr            T;
        kappaMethod     phaseSystem;
        value           uniform 300;
    }
    \endverbatim

    Solid region (e.g. 0/solid/T):
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureMultiphaseRadCoupledMixed;
        region          solid;
        Tnbr            T.liquid;
        kappaMethod     solidThermo;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on an underlying mapped(Wall)FvPatch.

SourceFiles
    turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField_H
#define turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "scalarList.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
       Class turbulentTemperatureMultiphaseRadCoupledMixed Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    typedef
        turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
        this_bctype;

    typedef mixedFvPatchScalarField parent_bctype;

public:

    // Public Enumerations

        //- Type of supplied kappa
        enum KMethodType
        {
            mtSolidThermo,
            mtLookup,
            mtPhaseSystem
        };

        //- Region this condition belongs to
        enum regionType
        {
            solid,
            fluid
        };


private:

    // Private Data

        static const Enum<regionType> regionTypeNames_;

        static const Enum<KMethodType> KMethodTypeNames_;

        regionType regionType_;

        const KMethodType method_;

        //- Name of thermal conductivity field (lookup only)
        const word kappaName_;

        //- Phases included in the wall balance (empty = all phases)
        wordList phaseNames_;

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Name of the radiative heat flux in the neighbour region
        const word qrNbrName_;

        //- Name of the radiative heat flux in local region
        const word qrName_;


    // Private Member Functions

        //- Conductivity on this patch, by method_
        tmp<scalarField> kappa(const scalarField& Tp) const;

        //- True if the phase takes part in the wall balance
        bool includesPhase(const word& phaseName) const;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureMultiphaseRadCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
        (
            const this_bctype&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
        (
            const this_bctype&,
            const DimensionedField<scalar, volMesh>&
        );

        //- No copy without an internal field
        turbulentTemperatureMultiphaseRadCoupledMixedFvPatchScalarField
        (
            const this_bctype&
        ) = delete;

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using parent_bctype::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
