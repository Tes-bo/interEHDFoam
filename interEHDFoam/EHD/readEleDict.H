Info<< "Reading electric settings"<< nl  << endl;

word phase1Name(word("centre.")+transportProperties.get<wordList>("phases")[0]);

IOdictionary electricProperties
(
    IOobject
    (
        "electricProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

int potentialNonorthogonalCorrectTime = electricProperties.subDict("settings").lookupOrDefault<int>("potentialNonorthogonalCorrectTime", 1);
Info << "potentialNonorthogonalCorrectTime:" << potentialNonorthogonalCorrectTime << nl << endl;

int electricFieldCorrectTime = electricProperties.subDict("settings").lookupOrDefault<int>("electricFieldCorrectTime", 1);
Info << "electricFieldCorrectTime:" << electricFieldCorrectTime << nl << endl;

word equationType(electricProperties.subDict("settings").lookup("equationType"));
Info << "equationType:" << equationType << nl << endl;

Info << "Reading electric physical properties:" << endl;
dimensionedScalar K1("K1", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0);
dimensionedScalar K2("K2", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0);
dimensionedScalar sig1("sig1", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0);
dimensionedScalar sig2("sig2", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0);
if (equationType=="NP")
{

    K1 = dimensionedScalar("K1", electricProperties);
    Info << K1 << endl;
    K2 = dimensionedScalar("K2", electricProperties);
    Info << K2 << endl;
}
else if (equationType=="Herrera" || equationType=="Herrera2" || equationType=="leakyDielectric")
{
    sig1 = dimensionedScalar("sig1", electricProperties);
    Info << sig1 << endl;
    sig2 = dimensionedScalar("sig2", electricProperties);
    Info << sig2 << endl;
}
else
{
    Info << "ERROR: Unknow equationType" << endl;
    Info << "equationType can be: NP, Herrera, leakyDielectric" << nl << endl;
    exit(0);
}
dimensionedScalar eps1("eps1", electricProperties);
Info << eps1 << endl;
dimensionedScalar eps2("eps2", electricProperties);
Info << eps2 << endl;
dimensionedScalar D1("D1", electricProperties);
Info << D1 << endl;
dimensionedScalar D2("D2", electricProperties);
Info << D2 << nl << endl;

word averageMethod(electricProperties.subDict("settings").lookup("averageMethod"));
Info << "averageMethod:" << averageMethod << nl  << endl;
if(averageMethod !="linear" && averageMethod !="harmonic")
{
    Info << "ERROR: Unknow average method" << endl;
    Info << "averageMethod can be: linear, harmonic" << nl << endl;
    exit(0);  
}

word electricFieldModel(electricProperties.subDict("settings").lookup("electricFieldModel"));
Info << "electricFieldModel:" << electricFieldModel << nl  << endl;
if(electricFieldModel != "Dis" && electricFieldModel != "gradient")
{
    Info << "ERROR: Unknow electricFieldModel" << endl;
    Info << "electricFieldModel can be: Dis, gradient" << nl << endl;
    exit(0);
}

word forceType(electricProperties.subDict("settings").lookup("forceType"));
Info << "forceType:" << forceType << nl  << endl;
if(forceType != "tenser" && forceType != "DIVtenser" && forceType != "volumetric")
{
    Info << "ERROR: Unknow forceType" << endl;
    Info << "force type can be: tenser,DIVtenser,volumetric" << nl << endl;
    exit(0);
}

Switch cutfaceInEquation(electricProperties.subDict("settings").lookupOrDefault("cutfaceInEquation", true));
Info << "cutfaceInEquation:" << cutfaceInEquation << nl << endl;

Switch cutfaceInForce(electricProperties.subDict("settings").lookupOrDefault("cutfaceInForce", true));
Info << "cutfaceInForce:" << cutfaceInForce << nl << endl;

Switch cutEps(electricProperties.subDict("settings").lookupOrDefault("cutEps", true));
Info << "cutEps:" << cutEps << nl << endl;

Switch useHeavisideToAverage(electricProperties.subDict("settings").lookupOrDefault("useHeavisideToAverage", false));
Info << "useHeavisideToAverage:" << useHeavisideToAverage << nl << endl;

scalar epsilon = electricProperties.subDict("settings").lookupOrDefault("mindx", 5e-7);
epsilon = epsilon * 1.5;
if(useHeavisideToAverage)
{
    Info << "epsilon:" << epsilon << nl << endl;
}

Switch limitFlux(electricProperties.subDict("settings").lookupOrDefault("limitFlux", false));
Info << "limitFlux:" << limitFlux << nl << endl;

Switch adaptiveCompressionFactor(electricProperties.subDict("settings").lookupOrDefault("adaptiveCompressionFactor", true));
Info << "adaptiveCompressionFactor:" << adaptiveCompressionFactor << nl << endl;
scalar compressionFactor = electricProperties.subDict("settings").lookupOrDefault("compressionFactor", 3.0);
if(adaptiveCompressionFactor==false)
{
    Info << "compressionFactor:" << compressionFactor << nl << endl;
}

Switch pureConduction(electricProperties.subDict("settings").lookupOrDefault("pureConduction", false));
Info << "pureConduction:" << pureConduction << nl << endl;

scalar jetSwitch = electricProperties.subDict("settings").lookupOrDefault("jetSwitch", 1.00);
