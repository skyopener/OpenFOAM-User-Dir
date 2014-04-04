    perm = permeability;

    forAll (perm, celli)
    {
        if (voidFraction[celli] > 1.0 - cVoid.value())
        {
            perm[celli] = interfacePermeability.value();
        }
    }

    volScalarField SpCoeff
    (
        "SpCoeff",
        rho*(1.0-voidFraction)/perm
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(rho, U)
      + fvm::Sp(SpCoeff, U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*
                    (
                        fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                      - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                    )
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
