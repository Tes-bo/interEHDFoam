if(electricFieldModel=="Dis")
{

    const labelUList &owners = mesh.owner();
    const labelUList &neighbours = mesh.neighbour();
    surfaceScalarField DisSuf =-fvc::snGrad(Ep)*mesh.magSf()*epsSuf;
    Dis=0*Dis;
    forAll(owners, facei)
    {
        Dis[owners[facei]] += ( DisSuf[facei]*(mesh.Cf()[facei]-mesh.C()[owners[facei]]) )/mesh.V()[owners[facei]];
        Dis[neighbours[facei]] -= ( DisSuf[facei]*(mesh.Cf()[facei]-mesh.C()[neighbours[facei]]) )/mesh.V()[owners[facei]];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], facei)
        {
            Dis[pFaceCells[facei]] +=  ( DisSuf.boundaryField()[patchi][facei]*(mesh.Cf().boundaryField()[patchi][facei]-mesh.C()[pFaceCells[facei]]) )/mesh.V()[pFaceCells[facei]];
        }
    }
    Dis.correctBoundaryConditions();
    E=Dis/eps;
    E.correctBoundaryConditions();

}
else
{
    E=-fvc::grad(Ep);
}

