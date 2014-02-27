{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));

    volScalarField alphaS("alphaS", 1-voidFraction);

    dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh.V()), 1.0/3.0)
    );

    // Calculate gradient of alpha function for alpha1-alpha2 interface
    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    surfaceScalarField phir
    (
        phic*
        (
            (gradAlphaf/(mag(gradAlphaf) + deltaN))&mesh.Sf()
        )
    );

    surfaceVectorField gradAlphaSf
    (
        fvc::interpolate(alphaS)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alphaS))
    );

    surfaceScalarField phirS
    (
        phic*
        (
            (gradAlphaSf/(mag(gradAlphaSf) + deltaN))&mesh.Sf()
        )
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
            // Interface compression flux at alpha1-alphaS interface
          + fvc::flux
            (
                -fvc::flux(-phirS, alphaS, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);

        alpha2 = voidFraction - alpha1;
        rhoPhi = phiAlpha*(rho1 - rho2)
               + phi*rho2
               + fvc::interpolate(alphaS)*phiS*(rhoS - rho2);
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
