#include"readEleDict.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Creating electric parameters"<< nl  << endl;

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*eps1
);

surfaceScalarField epsSuf
(
    IOobject
    (
        "epsSuf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1)*eps1
);


volScalarField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*D1
);

surfaceScalarField DSuf
(
    IOobject
    (
        "DSuf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1)*D1
);

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*K1
);

surfaceScalarField KSuf
(
    IOobject
    (
        "KSuf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1)*K1
);

volScalarField sig
(
    IOobject
    (
        "sig",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*sig1
);

surfaceScalarField sigSuf
(
    IOobject
    (
        "sigSuf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1)*sig1
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading field rhoe" << nl << endl;
volScalarField rhoe
(
    IOobject
    (
        "rhoe",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating elePotential"<< nl << endl;
volScalarField Ep
(
    IOobject
    (
        "Ep",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating eleField"<< nl << endl;
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   -fvc::grad(Ep)
);

volVectorField Dis
(
    IOobject
    (
        "Dis",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   E*eps,
   E.boundaryField().types()//因为我们后期是重建的Dis场，所以Dis的边界条件应该是通过内部的场外推边界的值，通过浏览代码可以发现，grad命令得到的场其边界条件就是外推，因此这里直接进行复制。
);

surfaceScalarField elecPhi
(
    IOobject
    (
        "elecPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::interpolate(K)*mesh.magSf()*fvc::snGrad(Ep)
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Creating eleForce "<< nl << endl;

tensor uni(1,0,0,0,1,0,0,0,1);

volTensorField Maxwell
(
    IOobject
    (
        "Maxwell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    eps*(E*E-0.5*magSqr(E)*uni)
);

volVectorField Fe
(
    IOobject
    (
        "Fe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(Maxwell)
);

volVectorField columbForce
(
    IOobject
    (
        "ColumbForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoe*E
);

volVectorField dieForce
(
    IOobject
    (
        "dieForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -0.5*mag(E)*mag(E)*fvc::grad(eps)
);

Info<< "Creating eleForce subdomain"<< nl << endl;

volScalarField FeFilter
(
    IOobject
    (
        "FeFilter",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField PhiFilter
(
    IOobject
    (
        "PhiFilte",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    FeFilter
);

forAll(PhiFilter,celli)
{
    if(mesh.C()[celli][0]>jetSwitch)
    {
        PhiFilter[celli]=0;
    }
    else
    {
        PhiFilter[celli]=1;
    }
}
PhiFilter.correctBoundaryConditions();

surfaceScalarField PhiFilterSuf
(
    IOobject
    (
        "PhiFilterSuf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(PhiFilter)
);
forAll(PhiFilterSuf,facei)
{
    if(mesh.Cf()[facei][0]>jetSwitch)
    {
        PhiFilterSuf[facei]=0;
    }
    else
    {
        PhiFilterSuf[facei]=1;
    }
}
Switch cutted = false;
cutted = false;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Creating surface cut parameters"<< nl  << endl;

volPointInterpolation interpolater(mesh);

pointScalarField alphaPoints
(
    IOobject
    (
        "alphaPoints",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    interpolater.interpolate(alpha1)
);

surfaceScalarField alphaSuf
(
    IOobject
    (
        "alphaSuf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1)
);

volScalarField RDFFunction
(
    IOobject
    (
        "RDFFunction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.lookupObject<volScalarField>("reconstructedDistanceFunction")
);

volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1,
    alpha1.boundaryField().types()
);

volScalarField interfaceCellMaker
(
    IOobject
    (
        "interfaceCellMaker",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
         "interfaceCellMaker",
         dimless,
         0
    )
);
volScalarField interfaceNeighbourMaker
(
    IOobject
    (
        "interfaceNeighbourMaker",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
         "interfaceNeighbourMaker",
         dimless,
         0
    )
);

surfaceScalarField phiCorrected
(
    IOobject
    (
        "phiCorrected",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi*alphaSuf
);

surfaceScalarField phiCompression
(
    IOobject
    (
        "phiCompression",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
         "phiCompression",
         dimensionSet(0,3,-1,0,0),
         0
    )
);

surfaceScalarField phiCompressionMarker
(
    IOobject
    (
        "phiCompressionMarker",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
         "phiCompressionMarker",
         dimensionSet(0,0,0,0,0),
         0
    )
);

volScalarField phiCompressionFactor
(
    IOobject
    (
        "phiCompressionFactor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
         "phiCompressionFactor",
         dimensionSet(0,0,0,0,0),
         0
    )
);

volScalarField rhoeOriCompreesion
(
    IOobject
    (
        "rhoeOriCompreesion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(phiCompression, rhoe) 
);

dimensionedScalar deltaT
(
    "deltaT",
    dimTime,
    runTime.deltaT().value()
);

volScalarField rhoeSubGrid
(
    IOobject
    (
        "rhoeSubGrid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   rhoe,
   rhoe.boundaryField().types()
);

volScalarField alpha1withInterface
(
    IOobject
    (
        "alpha1withInterface",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   alpha1,
   alpha1.boundaryField().types()
);

volScalarField alpha1Sub
(
    IOobject
    (
        "alpha1Sub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   (1-alpha1withInterface)*alpha1,
   alpha1.boundaryField().types()
);


volVectorField GradAlpha
(
    IOobject
    (
        "GradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(alpha1)
);


/*
volScalarField rhoeConvection
(
    IOobject
    (
        "rhoeConvection",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhoe
);
volScalarField rhoeOriCompreesionAfter
(
    IOobject
    (
        "rhoeOriCompreesionAfter",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(phiCompression, rhoe) 
);
*/
