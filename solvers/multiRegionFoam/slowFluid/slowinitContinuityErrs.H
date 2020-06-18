PtrList<uniformDimensionedScalarField> cumulativeContErrIO(slowFluidRegions.size());
forAll(cumulativeContErrIO, i)
{
    #include "setRegionSlowFluidFields.H"
    cumulativeContErrIO.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "cumulativeContErr",
                runTime.timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar(dimless, Zero)
        )
    );
}

UPtrList<scalar> cumulativeContErr(cumulativeContErrIO.size());
forAll(cumulativeContErrIO, i)
{
    cumulativeContErr.set(i, &cumulativeContErrIO[i].value());
}
