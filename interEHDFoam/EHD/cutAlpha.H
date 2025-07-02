/*
Note:
//Info<< alphaSuf.size() << endl; //internalFaceScalar
//Info<< mesh.faces().size() << endl; //Allfaces
*/
interfaceCellMaker = scalar(0)*interfaceCellMaker;
interfaceNeighbourMaker = scalar(0)*interfaceNeighbourMaker;
if (cutfaceInEquation || cutfaceInForce || limitFlux)
{
    alphaSuf = fvc::interpolate(H);
    alphaSuf = min(max(alphaSuf, scalar(0)), scalar(1));
    alphaPoints = interpolater.interpolate(alpha1);
    const labelUList &owner = mesh.owner();
    const labelUList &neighbour = mesh.neighbour();
    
    //cut the internal face:
    forAll(alphaSuf, facei)
    {
        if (alphaSuf[facei] > 0 && alphaSuf[facei] < 1)
        {
            labelList pLabels = mesh.faces()[facei];
            int nPoints = mesh.faces()[facei].size();
            int nLarger = 0;
            int nSmall = 0;
            int nEqual = 0;            
            forAll(pLabels, pi)
            {
                if (alphaPoints[pLabels[pi]] > 0.5)
                {
                    nLarger++;
                }
                else if(alphaPoints[pLabels[pi]] < 0.5)
                {
                    nSmall++;               
                }
                else
                {
                    nEqual++;
                }
            }
            if (nLarger == nPoints)
            {
                alphaSuf[facei] = scalar(1);
            }
            else if (nSmall == nPoints)
            {
                alphaSuf[facei] = scalar(0);            
            }
            else if (nEqual == nPoints)
            {
                if(0.499 < alphaSuf[facei] && alphaSuf[facei] < 0.501 )
                {
                    alphaSuf[facei]=scalar(0);
                    if(0.99<alpha1[owner[facei]])
                    {
                        interfaceCellMaker[owner[facei]]=interfaceCellMaker[owner[facei]]+scalar(2);
                    }
                    else
                    {
                        interfaceCellMaker[neighbour[facei]]=interfaceCellMaker[neighbour[facei]]+scalar(2);
                    }     
                }
                else if(alphaSuf[facei] >= 0.501)
                {
                    alphaSuf[facei]=scalar(1);
                }
                else
                {
                    alphaSuf[facei]=scalar(0);
                }
            }
            else if ((nLarger+nEqual) == nPoints)
            {
                alphaSuf[facei] = scalar(1);
                if(0<alpha1[owner[facei]] && alpha1[owner[facei]]<1)
                {
                    interfaceCellMaker[owner[facei]]=interfaceCellMaker[owner[facei]]+scalar(1);
                }
                else
                {
                    interfaceCellMaker[neighbour[facei]]=interfaceCellMaker[neighbour[facei]]+scalar(1);
                }            
            }
            else if ((nSmall+nEqual) == nPoints)
            {
                alphaSuf[facei] = scalar(0);
                if(0<alpha1[owner[facei]] && alpha1[owner[facei]]<1)
                {
                    interfaceCellMaker[owner[facei]]=interfaceCellMaker[owner[facei]]+scalar(1);
                }
                else
                {
                    interfaceCellMaker[neighbour[facei]]=interfaceCellMaker[neighbour[facei]]+scalar(1);
                }            
            }
            else
            {
                interfaceCellMaker[owner[facei]]=interfaceCellMaker[owner[facei]]+scalar(1);
                interfaceCellMaker[neighbour[facei]]=interfaceCellMaker[neighbour[facei]]+scalar(1);       
            }
        }
    }

    
    //cut the boundary face:
    forAll(alphaSuf.boundaryFieldRef(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        forAll(alphaSuf.boundaryFieldRef()[patchi], facei)
        {
            labelList pLabels = mesh.boundaryMesh()[patchi][facei];
            int nPoints = mesh.boundaryMesh()[patchi][facei].size();
            int nLarger = 0;
            int nSmall = 0;
            int nEqual = 0;    
            forAll(pLabels, pi)
            {
                if (alphaPoints[pLabels[pi]] > 0.5)
                {
                    nLarger++;
                }                
                else if (alphaPoints[pLabels[pi]] < 0.5)
                {
                    nSmall++;               
                }
                else
                {
                    nEqual++;
                }
 
            }
            if (nLarger == nPoints)
            {
                alphaSuf.boundaryFieldRef()[patchi][facei] = scalar(1);
            }
            else if (nSmall == nPoints)
            {
                alphaSuf.boundaryFieldRef()[patchi][facei] = scalar(0);           
            }
            else if (nEqual == nPoints)
            {
                if(0.499 < alphaSuf[facei] && alphaSuf[facei] < 0.501 )
                {
                    alphaSuf.boundaryFieldRef()[patchi][facei]=scalar(0);
                    if(0.99<alpha1[pFaceCells[facei]])
                    {
                        interfaceCellMaker[pFaceCells[facei]]=interfaceCellMaker[pFaceCells[facei]]+scalar(2);
                    }  
                }
                else if(alphaSuf[facei] >= 0.501)
                {
                    alphaSuf.boundaryFieldRef()[patchi][facei]=scalar(1);
                }
                else
                {
                    alphaSuf.boundaryFieldRef()[patchi][facei]=scalar(0);
                }
            }
            else if ((nLarger+nEqual) == nPoints)
            {
                alphaSuf.boundaryFieldRef()[patchi][facei] = scalar(1);
                if(0<alpha1[pFaceCells[facei]] && alpha1[pFaceCells[facei]]<1)
                {
                    interfaceCellMaker[pFaceCells[facei]]=interfaceCellMaker[pFaceCells[facei]]+scalar(1);
                }         
            }            
            else if ((nSmall+nEqual) == nPoints)
            {
                alphaSuf.boundaryFieldRef()[patchi][facei] = scalar(0);
                if(0<alpha1[pFaceCells[facei]] && alpha1[pFaceCells[facei]]<1)
                {
                    interfaceCellMaker[pFaceCells[facei]]=interfaceCellMaker[pFaceCells[facei]]+scalar(1);
                }          
            }           
            else
            {
                interfaceCellMaker[pFaceCells[facei]]=interfaceCellMaker[pFaceCells[facei]]+scalar(1);   
            }
        }
    }
    interfaceCellMaker.correctBoundaryConditions();
    // very important when dealing with parallel runing since it need to update the patchNeighbourField() value.
    
    if (limitFlux)
    {
        forAll(alphaSuf, facei)
        {
            if (alphaSuf[facei] == 0)
            {
                if (interfaceCellMaker[owner[facei]] > 1 && interfaceCellMaker[neighbour[facei]] < 2)
                {
                        interfaceNeighbourMaker[neighbour[facei]]=interfaceNeighbourMaker[neighbour[facei]]+1;
                }
                else if (interfaceCellMaker[neighbour[facei]] > 1 && interfaceCellMaker[owner[facei]] < 2)
                {
                        interfaceNeighbourMaker[owner[facei]]=interfaceNeighbourMaker[owner[facei]]+1;
                }      
            }
        }
        forAll(mesh.boundary(), patchi)
        {
            if(alphaSuf.boundaryFieldRef()[patchi].coupled())
            {
                const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
                const scalarField interfaceCellMakerNeighbour(interfaceCellMaker.boundaryField()[patchi].patchNeighbourField());        
                forAll(mesh.boundary()[patchi], facei)
                {
                    scalar alphaSuf_face( alphaSuf.boundaryField()[patchi][facei] );        
                    if (alphaSuf_face == 0 && interfaceCellMakerNeighbour[facei]  > 1 && interfaceCellMaker[pFaceCells[facei]] < 2)
                    {
                        interfaceNeighbourMaker[pFaceCells[facei]]=interfaceNeighbourMaker[pFaceCells[facei]]+1;                   
                    }
                }
            }    
        }
        interfaceNeighbourMaker.correctBoundaryConditions();    
    }  
}

if (limitFlux)
{
    forAll(alpha1, celli)
    {
        if (alpha1[celli]>0.5 || interfaceCellMaker[celli] > 1 || interfaceNeighbourMaker[celli] > 0 )
        {
            alpha1withInterface[celli]=scalar(1);
        }
        else
        {
            alpha1withInterface[celli]=scalar(0);
        }
    }
    alpha1withInterface.correctBoundaryConditions();

    alpha1Sub=(1-alpha1withInterface)*alpha1;
    forAll(alpha1Sub,celli)
    {
        if(alpha1Sub[celli]<1e-3)
        {
            alpha1Sub[celli]=0;
        }
    }
    alpha1Sub.correctBoundaryConditions();
}
