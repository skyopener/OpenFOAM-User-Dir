{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir(phic*interface.nHatf());

    surfaceScalarField phirS
    (
        phi / (fvc::interpolate(voidFraction) + 1.0 - cVoid)
    );

    volScalarField alphaS("alphaS", 1.0-voidFraction);

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
//        // Interface compression flux at alpha1-alphaS interface
//          + fvc::flux
//            (
//                -fvc::flux(-phirS, alphaS, alpharScheme),
//                alpha1,
//                alpharScheme
//            )
        );

        volScalarField aps(alpha1 + alphaS);

        MULES::explicitSolve(aps, phi, phiAlpha, 1, 0);

        alpha1 = aps - alphaS;

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
