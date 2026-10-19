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

Application
    chtMultiRegionMultiphaseEulerFoam_notGit

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent multiphase fluid flow and solid
    heat conduction, with conjugate heat transfer between solid and fluid
    regions.

    Each fluid region is a multiphaseSystem (any number of phases, common
    pressure), solved with the equations of reactingMultiphaseEulerFoam.
    Solid regions, region handling and time-step control are those of
    chtMultiRegionFoam / chtMultiRegionTwoPhaseEulerFoam.

    Running without CHT: list no solid regions in constant/regionProperties
    (solid ()) and put an ordinary heat-flux condition on the heated wall.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiphaseSystem.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "loopControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent multiphase fluid flow and"
        " solid heat conduction with conjugate heat transfer"
        " between solid and fluid regions."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createMultiphaseFluidFields.H"
    #include "createSolidFields.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "multiphaseMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readOuterPIMPLEControls.H"

        #include "multiphaseMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldMultiphaseFluidFields.H"
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr-1);

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "readMultiphaseFluidPIMPLEControls.H"
                #include "setRegionMultiphaseFluidFields.H"
                #include "solveMultiphaseFluid.H"
            }

            forAll(solidRegions, i)
            {
                fvMesh& mesh = solidRegions[i];

                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "setRegionSolidFields.H"
                #include "solveSolidRegion.H"
            }

            // Additional loops for energy solution only
            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        fvMesh& mesh = fluidRegions[i];

                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        #include "readMultiphaseFluidPIMPLEControls.H"
                        #include "setRegionMultiphaseFluidFields.H"
                        frozenFlow = true;
                        #include "solveMultiphaseFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        fvMesh& mesh = solidRegions[i];

                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "setRegionSolidFields.H"
                        #include "solveSolidRegion.H"
                    }
                }
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
