#include "readGravitationalAcceleration.H"

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

singlePhaseTransportModel laminarTransport(U, phi);

dimensionedScalar rhoInfValue
(
    "rhoInf",
    dimDensity,
    laminarTransport
);

volScalarField rhoInf
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoInfValue
);


volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    laminarTransport.nu()*rhoInfValue
);


Info<< "Constructing kinematicCloud "  << endl;
basicKinematicCloud kinematicCloud
(
    "kinematicCloud",
    rhoInf,
    U,
    mu,
    g
);


Info<< "Constructing kinematicCloud 2 "  << endl;

/* 
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
 */