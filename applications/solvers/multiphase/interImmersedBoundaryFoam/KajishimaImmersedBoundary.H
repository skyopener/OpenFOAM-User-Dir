{
    wordList corrUIBTypes
    (
        p_rgh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll (p_rgh.boundaryField(), i)
    {
        if (p_rgh.boundaryField()[i].fixesValue())
        {
            corrUIBTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField corrUIB
    (
        IOobject
        (
            "corrUIB",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("corrUIB", p_rgh.dimensions(), 0.0),
        corrUIBTypes
    );

    surfaceScalarField rAUf
    (
        "(1|A(U))",
        (
            upwind<scalar>(mesh, phiVoid).interpolate(voidFraction)
          + SMALL
        )*
        fvc::interpolate
        (

            dimensionedScalar("unity", dimTime, 1)/rho
        )
    );

    volVectorField Utemp
    (
        "Utemp",
        U
    );

    U = voidFraction*U + (scalar(1)-voidFraction)*US;

    U.correctBoundaryConditions();

    fp = rho*(U - Utemp) / runTime.deltaT();

    surfaceScalarField voidFractionf
    (
        fvc::interpolate(voidFraction)
    );

    phi = voidFractionf*phi
        + (scalar(1)-voidFractionf)*phiS;

    Info << "Original Kajishima model\n"
         << "Should not be divergence free\n" << endl;

    #include "continuityErrs.H"

    adjustPhi(phi, U, corrUIB);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix corrUIBEqn
        (
            fvm::laplacian(rAUf, corrUIB) == fvc::div(phi)
        );

        corrUIBEqn.setReference(pRefCell, pRefValue);
        corrUIBEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= corrUIBEqn.flux();
        }
    }

//    U = fvc::reconstruct(phi);
//    U = (U - (1-voidFraction)*US)/(voidFraction + dimensionedScalar("SM", dimless, 0.0001));

    Info << "Flux correction after IB\n"
         << "Should be divergence free now\n" << endl;

    #include "continuityErrs.H"
}
