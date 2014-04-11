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

    volScalarField a1 = min
    (
        alpha1/(voidFraction + dimensionedScalar("SM", dimless, SMALL)),
        dimensionedScalar("1", dimless, 1)
    );
    a1 = max(a1, dimensionedScalar("0", dimless, 0));

    volScalarField a2 = 1.0 - alpha1;

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
                        fvc::interpolate(a2)*fvc::snGrad(a1)
                      - fvc::interpolate(a1)*fvc::snGrad(a2)
                    )
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
