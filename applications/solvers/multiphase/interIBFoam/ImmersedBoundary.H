{
    volVectorField Ut
    (
        "Ut",
        U
    );

    surfaceScalarField phit
    (
        "phit",
        phi
    );

    U = voidFraction*Ut + (1.0 - voidFraction)*US;

    surfaceScalarField voidFractionf(fvc::interpolate(voidFraction));

    phi = voidFractionf*phit + (1.0 - voidFractionf)*phiS;

    fp = rho * (Ut - U) / runTime.deltaT();
}
