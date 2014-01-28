{
    interfaceFunc = alpha1/(voidFraction + dimensionedScalar("SMALL", dimless, SMALL));

    interfaceFunc = min(interfaceFunc, dimensionedScalar("1", dimless, 1));

    /*
    forAll(phiVoid, faceI)
    {
        label ownerI = mesh.owner()[faceI];
        label neighbourI = mesh.neighbour()[faceI];

        if
        (
            voidFraction[ownerI] < 0.01 && voidFraction[neighbourI] > 0.99
        )
        {
            interfaceFunc[ownerI] = alpha1[neighbourI];
        }
        else if
        (
            voidFraction[neighbourI] < 0.01 && voidFraction[ownerI] > 0.99
        )
        {
            interfaceFunc[neighbourI] = alpha1[ownerI];
        }
    }
    */

    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(alpha1));

    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    surfaceScalarField phiAlphaInverse
    (
        "phiAlphaInverse",
        -(gradAlphaf & mesh.Sf())
    );

    surfaceScalarField voidFractionf
    (
        linearInterpolate(voidFraction)
    );

    for(label ic = 0; ic < IBSmoothLoop; ic++)
    {
        surfaceScalarField interfaceFuncf
        (
            fvc::interpolate(interfaceFunc)
        );

        surfaceScalarField iff
        (
            mesh.magSf()*voidFractionf*interfaceFuncf
        );

        surfaceScalarField ifb
        (
            mesh.magSf()*voidFractionf
        );

        Field<scalar> ivf(voidFraction.size(), 0);
        Field<scalar> ivb(voidFraction.size(), 0);

        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        forAll(owner, facei)
        {
            ivf[owner[facei]] += iff[facei];
            ivb[owner[facei]] += ifb[facei];
            ivf[neighbour[facei]] += iff[facei];
            ivb[neighbour[facei]] += ifb[facei];
        }

        forAll(mesh.boundary(), patchi)
        {
            const labelUList& pFaceCells =
                mesh.boundary()[patchi].faceCells();

            const fvsPatchField<scalar>& piff = iff.boundaryField()[patchi];
            const fvsPatchField<scalar>& pifb = ifb.boundaryField()[patchi];

            forAll(mesh.boundary()[patchi], facei)
            {
                ivf[pFaceCells[facei]] += piff[facei];
                ivb[pFaceCells[facei]] += pifb[facei];
            }
        }

        forAll(interfaceFunc, celli)
        {
            if(voidFraction[celli] <= CVoid)
            {
                interfaceFunc[celli] = ivf[celli]/(ivb[celli]+VSMALL);

                interfaceFunc[celli] = max(interfaceFunc[celli], 0);
                interfaceFunc[celli] = min(interfaceFunc[celli], 1);
    //            if(interfaceFunc[celli] > 1-CVoid)
    //            {
    //                interfaceFunc[celli] = 1;
    //            }
    //            else
    //            {
    //                interfaceFunc[celli] = 0;
    //            }
            }
        }
    }

    interfaceFunc.correctBoundaryConditions();
}
