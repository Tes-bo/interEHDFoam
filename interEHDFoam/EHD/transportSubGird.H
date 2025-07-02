forAll(rhoe,celli)
{
    if(rhoe[celli]>1e-5)
    {
        rhoeSubGrid[celli]=rhoe[celli]-alpha1withInterface[celli]*rhoe[celli];    
        rhoe[celli]=alpha1withInterface[celli]*rhoe[celli];    
    }
    else
    {
        rhoeSubGrid[celli]=scalar(0);       
    }
}

rhoeSubGrid.correctBoundaryConditions();
rhoeSubGrid.oldTime()=rhoeSubGrid;
fvScalarMatrix rhoeSubEqn
(
    fvm::ddt(rhoeSubGrid) 
    + fvm::div(phi, rhoeSubGrid)        
    ==
    alpha1Sub*fvc::laplacian(fvc::interpolate((sig1 * sig2) / ((alpha1Sub * sig2) + ((1 - alpha1Sub) * sig1))), Ep)
);

rhoeSubEqn.solve();
rhoeSubGrid.correctBoundaryConditions();

rhoe=rhoe+rhoeSubGrid;
rhoe.correctBoundaryConditions();
