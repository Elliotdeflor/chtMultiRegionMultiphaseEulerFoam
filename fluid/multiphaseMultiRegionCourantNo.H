// Courant number over all fluid regions: the larger of the mixture flux,
// every phase flux, and every pair's relative flux (generalises the
// phi1 - phi2 term of the two-phase CHT solver).

    scalar CoNum = -GREAT;

    forAll(fluidRegions, regioni)
    {
        if (fluidRegions[regioni].nInternalFaces())
        {
            multiphaseSystem& fluidCo = phaseSystemFluid[regioni];
            multiphaseSystem::phaseModelList& phasesCo = fluidCo.phases();

            scalarField sumPhi
            (
                fvc::surfaceSum(mag(fluidCo.phi()))().primitiveField()
            );

            forAll(phasesCo, phasei)
            {
                sumPhi = max
                (
                    sumPhi,
                    fvc::surfaceSum
                    (
                        mag(phasesCo[phasei].phi())
                    )().primitiveField()
                );

                for (label phasej = phasei + 1; phasej < phasesCo.size(); ++phasej)
                {
                    sumPhi = max
                    (
                        sumPhi,
                        fvc::surfaceSum
                        (
                            mag(phasesCo[phasei].phi() - phasesCo[phasej].phi())
                        )().primitiveField()
                    );
                }
            }

            const scalar regionCoNum =
                0.5*gMax
                (
                    sumPhi/fluidRegions[regioni].V().field()
                )*runTime.deltaTValue();

            CoNum = Foam::max(CoNum, regionCoNum);
        }
    }

    Info<< "Courant Number max: " << CoNum << endl;
