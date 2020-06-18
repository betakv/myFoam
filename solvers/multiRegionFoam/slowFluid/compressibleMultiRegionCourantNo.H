    scalar CoNum = -GREAT;

    forAll(slowFluidRegions, regioni)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                slowFluidRegions[regioni],
                runTime,
                rhoSlowFluid[regioni],
                phiSlowFluid[regioni]
            ),
            CoNum
        );
    }
/*
    forAll(porousFluidRegions, porousi)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                porousFluidRegions[porousi],
                runTime,
                rhoPorous[porousi],
                phiPorous[porousi]
            ),
            CoNum
        );
    }
*/
