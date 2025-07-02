H = alpha1;

if (useHeavisideToAverage)
{
    //-NOTE
    RDFFunction = mesh.lookupObject<volScalarField>("reconstructedDistanceFunction");

    forAll(mesh.cells(), celli)
    {
        if (RDFFunction[celli] == 0)
        {

            if (alpha1[celli] > 0.5)
            {
                H[celli] = scalar(1);
            }
            else
            {
                H[celli] = scalar(0);
            }
        }
        else if (RDFFunction[celli] < -epsilon)
        {
            H[celli] = double(0);
        }
        else if (epsilon < RDFFunction[celli])
        {
            H[celli] = double(1);
        }
        else
        {
            H[celli] = double(1.0) / double(2.0) * (double(1.0) + RDFFunction[celli] / epsilon + Foam::sin(M_PI * RDFFunction[celli] / epsilon) / M_PI);
        }
    }
    H.correctBoundaryConditions();
    H = min(max(H, scalar(0)), scalar(1));
}

if (equationType=="NP")
{

    if(averageMethod=="linear")
    {
        K = H * K1 + (1 - H) * K2;
        eps = H * eps1 + (1 - H) * eps2;
        D = H * D1 + (1 - H) * D2;
    }
    else if(averageMethod=="harmonic")
    {
        K = (K1 * K2) / ((H * K2) + ((1 - H) * K1));
        eps = (eps1 * eps2) / ((H * eps2) + ((1 - H) * eps1));
        D = (D1 * D2) / ((H * D2) + ((1 - H) * D1));
    }
    /*
    KSuf = fvc::interpolate(K);
    epsSuf = fvc::interpolate(eps);
    DSuf = fvc::interpolate(D);
    */
}
else
{
    if(averageMethod=="linear")
    {
        sig = H * sig1 + (1 - H) * sig2;
        eps = H * eps1 + (1 - H) * eps2;
    }
    else if(averageMethod=="harmonic")
    {
        sig = (sig1 * sig2) / ((H * sig2) + ((1 - H) * sig1));
        eps = (eps1 * eps2) / ((H * eps2) + ((1 - H) * eps1));
    }
    /*
    sigSuf = fvc::interpolate(sig);
    epsSuf = fvc::interpolate(eps);
    */
}
