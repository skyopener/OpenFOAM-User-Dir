{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir(phic*interface.nHatf());

    volVectorField UrS(fvc::reconstruct(phi));
    UrS -= US;

    volVectorField gradAlphaS
    (
        -fvc::grad(1-voidFraction)
    );

    surfaceVectorField gradAlphaSf(fvc::interpolate(gradAlphaS));

    dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh.V()), 1.0/3.0)
    );

    surfaceVectorField nHatSf(gradAlphaSf/(mag(gradAlphaSf) + deltaN));

    surfaceScalarField nHatS(nHatSf & mesh.Sf());

    surfaceScalarField phirS
    (
        "phirS",
        mag(fvc::interpolate(UrS))*nHatS
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
