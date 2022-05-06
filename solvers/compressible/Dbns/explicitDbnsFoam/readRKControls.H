dictionary RKDict = mesh.solutionDict().subDict("RK");

List<scalar> alpha( RKDict.lookup("alpha") ); 

bool LTS = 
    RKDict.lookupOrDefault("localTimestepping", false);

IRSSmooth IRS(mesh);
