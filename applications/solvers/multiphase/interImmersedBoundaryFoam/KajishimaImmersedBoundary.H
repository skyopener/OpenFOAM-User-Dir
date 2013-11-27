{
    phi = fvc::interpolate(voidFraction)*phi
        + (scalar(1)-fvc::interpolate(voidFraction))*phiS;

    U = voidFraction*U + (scalar(1)-voidFraction)*US;

    Info << "Original Kajishima model\n"
         << "Should not be divergence free\n" << endl;

    #include "continuityErrs.H"

    wordList pcorrTypes
    (
        p_rgh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll (p_rgh.boundaryField(), i)
    {
        if (p_rgh.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p_rgh.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf("(1|A(U))", dimTime/rho.dimensions(), 1.0);

    adjustPhi(phi, U, pcorr);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

    #include "continuityErrs.H"
}
