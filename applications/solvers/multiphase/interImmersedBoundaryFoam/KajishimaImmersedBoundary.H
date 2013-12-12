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

    dimensionedScalar rAUf("(1|A(U))", dimTime/rho.dimensions(), 1.0);

    U = voidFraction*U + (scalar(1)-voidFraction)*US;

    U.correctBoundaryConditions();

    surfaceScalarField voidFractionf
    (
        fvc::interpolate(voidFraction)
    );


    for(label ibI = 0; ibI < nIBLoop; ibI++)
    {
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

        #include "continuityErrs.H"
    }
}
