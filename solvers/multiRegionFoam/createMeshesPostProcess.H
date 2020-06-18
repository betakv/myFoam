#include "createMeshes.H"

if (!slowFluidRegions.size() && !solidRegions.size())
{
    FatalErrorIn(args.executable())
        << "No region meshes present" << exit(FatalError);
}

fvMesh& mesh = slowFluidRegions.size() ? slowFluidRegions[0] : solidRegions[0];
