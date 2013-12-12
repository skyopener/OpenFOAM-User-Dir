label nAlphaCorr(readLabel(pimple.dict().lookup("nAlphaCorr")));
label nAlphaSubCycles(readLabel(pimple.dict().lookup("nAlphaSubCycles")));

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum(0.0*rhoPhi);

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;

rho = voidFraction*rho + (1-voidFraction)*rhoS;

rhoPhi = fvc::interpolate(voidFraction)*rhoPhi
       + (1-fvc::interpolate(voidFraction))*rhoS*phiS;
