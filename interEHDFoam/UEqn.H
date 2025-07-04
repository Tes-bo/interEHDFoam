    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );


    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    //mixture.surfaceTensionForce()
                  surfForces.surfaceTensionForce()
                  //- ghf*fvc::snGrad(rho)
                  //- NOTE
                  //+ surfForces.accelerationForce(rho)
                  + surfForces.accelerationForce()
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )+Fe
        );

        fvOptions.correct(U);
    }
