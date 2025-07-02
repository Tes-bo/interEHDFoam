//速度可以从0的网格流向0.5的网格（由小到大,但不能从1的网格流向0.5的网格（由大到小）
phiCorrected = phi;
phiCompression = 0*phiCompression;
//phiCompressionMarker=0*phiCompressionMarker;
if (limitFlux && (! pureConduction)) 
{
    phiCorrected = phiCorrected*alphaSuf;
    const labelUList &owner = mesh.owner();
    const labelUList &neighbour = mesh.neighbour(); 
    forAll(alphaSuf, facei)
    //修正内部通量
    {
        if (alphaSuf[facei] == 0)
        //当前面相分数值为0
        {
            if (interfaceCellMaker[owner[facei]] > 1 && interfaceCellMaker[neighbour[facei]] < 2 && phi[facei] < 0)
            //当前面为0，此时如果这个面的ownerCell是界面Cell(即alpha=0.5的cell)，则此时的面法向矢量Sf是从0.5指向0,如果此时的通量小于0，说明速度和Sf相反，即U从0指向0.5（由小到大），需要加倍
            {
                    phiCompression[facei] = phi[facei];
                    //phiCompressionMarker[facei]=scalar(-1);
                    //interfaceNeighbourMaker[neighbour[facei]]=interfaceNeighbourMaker[neighbour[facei]]+1;
            }
            else if (interfaceCellMaker[neighbour[facei]] > 1 && interfaceCellMaker[owner[facei]] < 2 && phi[facei] > 0)
            //当前面为0，此时如果这个面的neighbourCell是界面Cell(即alpha=0.5的cell)，则此时的面法向矢量Sf是从0指向0.5
            {
                    phiCompression[facei] =  phi[facei];
                    //phiCompressionMarker[facei]=scalar(1);
                    //interfaceNeighbourMaker[owner[facei]]=interfaceNeighbourMaker[owner[facei]]+1;
            }      
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        if(alphaSuf.boundaryFieldRef()[patchi].coupled())
        //处理边界条件
        {
            const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
            const scalarField interfaceCellMakerNeighbour(interfaceCellMaker.boundaryField()[patchi].patchNeighbourField());
            const scalarField alpha1Neighbour(alpha1.boundaryField()[patchi].patchNeighbourField());
                     
            forAll(mesh.boundary()[patchi], facei)
            {
                scalar alphaSuf_face( alphaSuf.boundaryField()[patchi][facei] );
                scalar phi_face( phi.boundaryField()[patchi][facei] );
                
                if (alphaSuf_face == 0)
                //当前面为0
                {
                    if (interfaceCellMaker[pFaceCells[facei]] > 1 && interfaceCellMakerNeighbour[facei] < 2 && phi_face < 0)
                    //当前面为0，此时如果这个面的ownerCell是界面Cell(即alpha=0.5的cell)，则此时的面法向矢量Sf是从0.5指向0,如果此时的通量小于0，说明速度和Sf相反，即U从0指向0.5（由小到大），需要加倍
                    {
                        phiCompression.boundaryFieldRef()[patchi][facei] = phi_face;
                        //phiCompressionMarker.boundaryFieldRef()[patchi][facei]=scalar(-1);
                    }
                    else if (interfaceCellMakerNeighbour[facei]  > 1 && interfaceCellMaker[pFaceCells[facei]] < 2 && phi_face > 0)
                    //当前面为0，此时如果这个面的neighbourCell是界面Cell(即alpha=0.5的cell)，则此时的面法向矢量Sf是从0指向0.5
                    {
                        phiCompression.boundaryFieldRef()[patchi][facei] = phi_face;
                        //phiCompressionMarker.boundaryFieldRef()[patchi][facei]=scalar(1); 
                        //interfaceNeighbourMaker[pFaceCells[facei]]=interfaceNeighbourMaker[pFaceCells[facei]]+1;                   
                    }  
                }
            }
        }    
    }
   
}
