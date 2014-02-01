{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir(phic*interface.nHatf());

    volVectorField UrS
    (
        "UrS",
        U - US
    );

    surfaceScalarField phirS
    (
        "phirS",
        fvc::interpolate(UrS) & mesh.Sf()
    );

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
          + fvc::flux
            (
                fvc::flux(phirS, 1-voidFraction, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        userMULES::explicitSolve(alpha1, voidFraction, phi, phiAlpha, 1, 0);

        alpha2 = voidFraction - alpha1;
        rhoPhi = phiAlpha*(rho1 - rho2)
               + phi*rho2
               + fvc::interpolate(1-voidFraction)*phiS*(rhoS - rho2);
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
