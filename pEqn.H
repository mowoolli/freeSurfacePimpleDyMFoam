/*---------------------------------------------------------------------------*\ 
| Modified 2010-2012 Copyright (C) Engys Ltd                                  |
\*---------------------------------------------------------------------------*/
{
//    volScalarField rAU(1.0/UEqn().A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();

    if (pimple.nCorrPISO() <= 1)
    {
        UEqn.clear();
    }

    surfaceScalarField phiHbyA("phiHbyA", (fvc::interpolate(HbyA) & mesh.Sf()) + rAUf*fvc::ddtCorr(U, phi));

    fvOptions.makeRelative(phiHbyA);
    
    adjustPhi(phiHbyA, U, p);

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p.boundaryField(),
        (
            phiHbyA.boundaryField()
          - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
