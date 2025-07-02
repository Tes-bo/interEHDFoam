Maxwell = eps * (E * E - 0.5 * magSqr(E) * uni);
Maxwell.correctBoundaryConditions();

GradAlpha=fvc::grad(alpha1);

//FeFilter
forAll(FeFilter,celli)
{
    if(-1<mag(GradAlpha[celli]) && mag(GradAlpha[celli])<1)
    {
        FeFilter[celli]=0;
    }
    else
    {
        FeFilter[celli]=1;
    }
}
FeFilter.correctBoundaryConditions();

if(cutfaceInForce)
{
    #include "cutEleProperties.H"    
}
else
{
    epsSuf = fvc::interpolate(eps);
}
if(forceType=="tenser")
{

    columbForce = FeFilter*fvc::surfaceIntegrate
                    (
                    fvc::interpolate(Dis)*(fvc::interpolate(Dis)&mesh.Sf()/epsSuf)
                    );

    dieForce = FeFilter*fvc::surfaceIntegrate
                    (
                    -fvc::interpolate(Dis&Dis)/2/epsSuf*mesh.Sf()
                    );
}
else if(forceType=="DIVtenser")
{
    columbForce = FeFilter*fvc::surfaceIntegrate(epsSuf * fvc::dotInterpolate(mesh.Sf(),(E * E)));
    dieForce = FeFilter*fvc::surfaceIntegrate(epsSuf * fvc::dotInterpolate(mesh.Sf(),(-0.5 * magSqr(E) * uni)));
}
else
{
    columbForce = FeFilter * rhoe * E;
    dieForce = -0.5 * mag(E) * mag(E) * fvc::surfaceIntegrate(mesh.Sf() * epsSuf) * FeFilter;
}

columbForce = FeFilter*columbForce + (1-FeFilter)*rhoe * E;
dieForce = FeFilter*dieForce - (1-FeFilter) * 0.5 * mag(E) * mag(E) * fvc::surfaceIntegrate(mesh.Sf() * epsSuf);

columbForce.correctBoundaryConditions();
dieForce.correctBoundaryConditions(); 
Fe = columbForce + dieForce;
