    Info<< "Reading immersedBoundaryProperties\n" << endl;
    IOdictionary immersedBoundaryProperties
    (
        IOobject
        (
            "immersedBoundaryProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rhoS
    (
        immersedBoundaryProperties.lookup("rhoS")
    );

    scalar CVoid
    (
        readScalar
        (
            immersedBoundaryProperties.lookup("CVoid")
        )
    );

    label IBSmoothLoop
    (
        readLabel
        (
            immersedBoundaryProperties.lookup("IBSmoothLoop")
        )
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field voidFraction\n" << endl;
    volScalarField voidFraction
    (
        IOobject
        (
            "voidFraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField gradVoid
    (
        fvc::grad(linearInterpolate(voidFraction))
    );

    surfaceVectorField gradVoidf
    (
        fvc::interpolate(gradVoid)
    );

    // phiVoid is used for upwind interpolation of voidFraction
    surfaceScalarField phiVoid
    (
        "phiVoid",
        gradVoidf & mesh.Sf()
    );

    Info<< "Reading field US\n" << endl;
    volVectorField US
    (
        IOobject
        (
            "US",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhiS.H"

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;
    incompressibleTwoPhaseMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());
    volScalarField& alpha2(twoPhaseProperties.alpha2());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    // interfaceFunc is used for curvature calculation
    // Initialisation does not matter because interfaceFunc
    // is reset after the alpha1 solution.
    Info<< "Creating interfaceFunc\n" << endl;
    volScalarField interfaceFunc
    (
        IOobject
        (
            "interfaceFunc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
    );

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    // Construct interface from alpha1 distribution
    interfaceImmersedBoundaryProperties interface
    (
        interfaceFunc,
        voidFraction,
        phiVoid,
        U,
        twoPhaseProperties,
        immersedBoundaryProperties
    );


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    volVectorField fp
    (
        IOobject
        (
            "fp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimAcceleration*dimDensity, vector::zero)
    );

    fv::IOoptionList fvOptions(mesh);
