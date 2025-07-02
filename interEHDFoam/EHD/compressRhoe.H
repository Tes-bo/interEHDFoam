//rhoeConvection=rhoe;

if(adaptiveCompressionFactor)
{
    rhoeOriCompreesion=fvc::div(phiCompression, rhoe); 
    deltaT=runTime.deltaT();
    forAll(rhoe, celli)
    {
        if(rhoeOriCompreesion[celli] <-1e-10  ||  1e-10 < rhoeOriCompreesion[celli])
        {
            phiCompressionFactor[celli]=rhoe[celli]/deltaT.value()/rhoeOriCompreesion[celli];
        }
        else
        {
            phiCompressionFactor[celli]=0;
        }
    }
    phiCompressionFactor.correctBoundaryConditions();

    forAll(phiCompression, facei)
    {
        const labelUList &owner = mesh.owner();
        const labelUList &neighbour = mesh.neighbour();

        if (interfaceNeighbourMaker[owner[facei]]>0)
        {
            phiCompression[facei]=phiCompression[facei]*phiCompressionFactor[owner[facei]];
        }
        else if (interfaceNeighbourMaker[neighbour[facei]]>0)
        {
            phiCompression[facei]=phiCompression[facei]*phiCompressionFactor[neighbour[facei]];    
        }
        
        /*
        if (phiCompressionMarker[facei]==1)
        {
            phiCompression[facei]=phiCompression[facei]*phiCompressionFactor[owner[facei]];
        }
        else if (phiCompressionMarker[facei]==-1)
        {
            phiCompression[facei]=phiCompression[facei]*phiCompressionFactor[neighbour[facei]];    
        }
        */
    }

    forAll(mesh.boundary(), patchi)
    {
        if(alphaSuf.boundaryFieldRef()[patchi].coupled())
        //处理边界条件
        {
            const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
            const scalarField phiCompressionFactorNeighbour(phiCompressionFactor.boundaryField()[patchi].patchNeighbourField());   
            const scalarField interfaceNeighbourMakerNeighbour(interfaceNeighbourMaker.boundaryField()[patchi].patchNeighbourField());          
            forAll(mesh.boundary()[patchi], facei)
            {
                if(interfaceNeighbourMaker[pFaceCells[facei]]>0)
                {
                    phiCompression.boundaryFieldRef()[patchi][facei]=phiCompression.boundaryField()[patchi][facei]*phiCompressionFactor[pFaceCells[facei]];
                }
                else if(interfaceNeighbourMakerNeighbour[facei]>0)
                {
                    phiCompression.boundaryFieldRef()[patchi][facei]=phiCompression.boundaryField()[patchi][facei]*phiCompressionFactorNeighbour[facei];
                }
                /*
                if(phiCompressionMarker.boundaryField()[patchi][facei]==1)
                {
                    phiCompression.boundaryFieldRef()[patchi][facei]=phiCompression.boundaryField()[patchi][facei]*phiCompressionFactor[pFaceCells[facei]];
                }
                else if(phiCompressionMarker.boundaryField()[patchi][facei]==-1)
                {
                    phiCompression.boundaryFieldRef()[patchi][facei]=phiCompression.boundaryField()[patchi][facei]*phiCompressionFactorNeighbour[facei];
                }
                */
            
            }
        }
    }
    //rhoeOriCompreesionAfter=fvc::div(phiCompression, rhoe);
}
else
{
   phiCompression=phiCompression*compressionFactor;
}

rhoe=rhoe-fvc::div(phiCompression, rhoe)*deltaT*PhiFilter;

rhoe.correctBoundaryConditions();
