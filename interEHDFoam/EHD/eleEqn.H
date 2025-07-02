if(pureConduction)
{
    phiCorrected=phiCorrected*0;
    phiCompression=phiCompression*0;
}

if(cutfaceInEquation)
{
    #include "cutEleProperties.H"
}
else
{
    if (equationType=="NP")
    {
        KSuf = fvc::interpolate(K);
        epsSuf = fvc::interpolate(eps);
        DSuf = fvc::interpolate(D);
    }
    else
    {
        sigSuf = fvc::interpolate(sig);
        epsSuf = fvc::interpolate(eps);
    }
}
for( int b = 0; b < electricFieldCorrectTime; b = b + 1 )
{
    if (equationType=="NP")
    {   
        for (int c = 0; c < potentialNonorthogonalCorrectTime; c = c + 1)
        {
            fvScalarMatrix EpEqn
            (
                fvm::laplacian(epsSuf, Ep) ==
                -rhoe
            );
            EpEqn.solve();
            Ep.correctBoundaryConditions();
        }

        elecPhi = -KSuf * mesh.magSf() * fvc::snGrad(Ep);

        fvScalarMatrix rhoeEqn
        (
            fvm::ddt(rhoe) 
            + fvm::div(elecPhi, rhoe) 
            + fvm::div(phiCorrected, rhoe) 
            + fvm::div(phiCompression, rhoe) 
            - fvm::laplacian(DSuf, rhoe)
        );
        rhoeEqn.solve();
        rhoe.correctBoundaryConditions();
    }

    else if (equationType=="Herrera")
    {   
        for (int c = 0; c < potentialNonorthogonalCorrectTime; c = c + 1)
        {
            fvScalarMatrix EpEqn
            (
                fvm::laplacian(epsSuf, Ep) ==
                -rhoe
            );
            EpEqn.solve();
            Ep.correctBoundaryConditions();
        }
        phiCorrected=phiCorrected*PhiFilterSuf+(1-PhiFilterSuf)*phi;
        fvScalarMatrix rhoeEqn
        (
            fvm::ddt(rhoe) 
            + fvm::div(phiCorrected, rhoe)       
            ==
            fvc::laplacian(sigSuf, Ep)
        );
        rhoeEqn.solve();
        rhoe.correctBoundaryConditions();
        if (limitFlux)
        {
            #include "compressRhoe.H"
            //#include "transportSubGird.H"  
        }    
    }
    
    else if (equationType=="Herrera2")
    {   
        for (int c = 0; c < potentialNonorthogonalCorrectTime; c = c + 1)
        {
            fvScalarMatrix EpEqn
            (
                fvm::laplacian(epsSuf, Ep) ==
                -rhoe
            );
            EpEqn.solve();
            Ep.correctBoundaryConditions();
        }
        phiCorrected=phi*fvc::interpolate(alpha1);
        fvScalarMatrix rhoeEqn
        (
            fvm::ddt(rhoe) 
            + fvm::div(phiCorrected, rhoe)          
            ==
            fvc::laplacian(sigSuf, Ep)
        );
        rhoeEqn.solve();
        rhoe.correctBoundaryConditions();
    }

    else if (equationType=="leakyDielectric")
    {

        for (int c = 0; c < potentialNonorthogonalCorrectTime; c = c + 1)
        {
            fvScalarMatrix EpEqn
            (
                fvm::laplacian(sigSuf, Ep)
            );
            EpEqn.solve();
            Ep.correctBoundaryConditions();
        }
        rhoe = -fvc::laplacian(epsSuf, Ep);
        rhoe.correctBoundaryConditions();
    }
}

#include "reconstructEle.H"
