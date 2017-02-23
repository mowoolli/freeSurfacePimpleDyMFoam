    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(rho, U)
// add for fvOption support
     ==
        fvOptions(rho, U)
//
    );

    UEqn.relax();

// add for fvOption support
//    fvOptions.constrain(UEqn());
//

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

// add for fvOption support
        fvOptions.correct(U);
//
    }
