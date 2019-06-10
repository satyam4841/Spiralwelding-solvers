Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field heat source (g)\n" << endl;

volScalarField g
(
    IOobject
    (
        "g",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading q: \n" << endl;

dimensionedScalar q
(
    "q",
    dimless,
    transportProperties.lookup("q")
);

Info<< "Reading radius of gaussian heat source r:\n" << endl;

dimensionedScalar r
(
    "r",
    dimless,
    transportProperties.lookup("r")
);

Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT
(
    "DT",
    dimArea/dimTime,
    transportProperties
);

#include "createFvOptions.H"
