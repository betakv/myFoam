Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

Info<< "Reading soiling properties\n" << endl;

Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field V\n" << endl;
volVectorField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field phiV\n" << endl;
surfaceScalarField phiV
(
    IOobject
    (
        "phiV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(V) & mesh.Sf()
);

Info<< "Reading soiling properties\n" << endl;
IOdictionary soilingProperties
(
    IOobject
    (
        "soilingProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar d( soilingProperties.lookup("d") );// mean droplet diameters
Info << "Mean droplet diameter " << d.value() << endl;

dimensionedScalar rhoSoil(soilingProperties.lookup("rhoSoil")); // droplet density
Info << "Soil density " << rhoSoil.value() << endl;

dimensionedScalar rhoAir(soilingProperties.lookup("rhoAir")); // air density
Info << "Air density " << rhoAir.value() << endl;

dimensionedScalar Sct(soilingProperties.lookupOrDefault("Sct",0.7)); // Schmidt turbulent number
Info << "Schmidt turbulent number " << Sct.value() << endl; 
  
dimensionedVector g( soilingProperties.lookup("g") );// gravity forces
Info << "Gravity vector " << g << endl;

dimensionedScalar twoWayCoupling(soilingProperties.lookup("twoWayCoupling")); // two way coupling 1 yes 0 no
Info << "Two way coupling " << twoWayCoupling.value() << endl;

Info<< "Reading additional difusivity\n" << endl;
volSymmTensorField delta
(
    IOobject
    (
        "delta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

