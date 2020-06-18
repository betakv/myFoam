// Initialise fluid field pointer lists
PtrList<rhoReactionThermo> slowThermoFluid(slowFluidRegions.size());
PtrList<volScalarField> rhoSlowFluid(slowFluidRegions.size());
PtrList<volVectorField> USlowFluid(slowFluidRegions.size());
PtrList<surfaceScalarField> phiSlowFluid(slowFluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefSlowFluid(slowFluidRegions.size());
PtrList<volScalarField> ghSlowFluid(slowFluidRegions.size());
PtrList<surfaceScalarField> ghfSlowFluid(slowFluidRegions.size());
PtrList<compressible::turbulenceModel> turbulenceSlowFluid(slowFluidRegions.size());
PtrList<CombustionModel<rhoReactionThermo>> reactionSlowFluid(slowFluidRegions.size());
PtrList<volScalarField> p_rghSlowFluid(slowFluidRegions.size());
PtrList<radiation::radiationModel> radiationSlow(slowFluidRegions.size());
PtrList<volScalarField> KSlowFluid(slowFluidRegions.size());
PtrList<volScalarField> dpdtSlowFluid(slowFluidRegions.size());
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsSlowFluid(slowFluidRegions.size());
PtrList<volScalarField> QdotSlowFluid(slowFluidRegions.size());

List<scalar> initialMassSlowFluid(slowFluidRegions.size());
List<bool> frozenFlowSlowFluid(slowFluidRegions.size(), false);

PtrList<IOMRFZoneList> MRFSlowFluid(slowFluidRegions.size());
PtrList<fv::options> slowFluidFvOptions(slowFluidRegions.size());

List<label> pRefCellSlowFluid(slowFluidRegions.size());
List<scalar> pRefValueSlowFluid(slowFluidRegions.size());

PtrList<dimensionedScalar> rhoMinSlowFluid(slowFluidRegions.size());
PtrList<dimensionedScalar> rhoMaxSlowFluid(slowFluidRegions.size());

PtrList<pressureControl> pressureControlsSlow(slowFluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(slowFluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << slowFluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;
    slowThermoFluid.set(i, rhoReactionThermo::New(slowFluidRegions[i]).ptr());

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoSlowFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                slowFluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            slowThermoFluid[i].rho()
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    USlowFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                slowFluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            slowFluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiSlowFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                slowFluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoSlowFluid[i]*USlowFluid[i])
                & slowFluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefSlowFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                slowFluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero) // uses name
        )
    );

    dimensionedScalar ghSlowRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefSlowFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghSlowFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & slowFluidRegions[i].C()) - ghSlowRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfSlowFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & slowFluidRegions[i].Cf()) - ghSlowRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceSlowFluid.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoSlowFluid[i],
            USlowFluid[i],
            phiSlowFluid[i],
            slowThermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to reactionFluid\n" << endl;
    reactionSlowFluid.set
    (
        i,
        CombustionModel<rhoReactionThermo>::New
        (
            slowThermoFluid[i],
            turbulenceSlowFluid[i]
        )
    );

    p_rghSlowFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                slowFluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            slowFluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghSlowFluid[i] = slowThermoFluid[i].p() - rhoSlowFluid[i]*ghSlowFluid[i];

    slowFluidRegions[i].setFluxRequired(p_rghSlowFluid[i].name());

    Info<< "    Adding to radiationFluid\n" << endl;
    radiationSlow.set
    (
        i,
        radiation::radiationModel::New(slowThermoFluid[i].T())
    );

    initialMassSlowFluid[i] = fvc::domainIntegrate(rhoSlowFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KSlowFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(USlowFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtSlowFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                slowFluidRegions[i]
            ),
            slowFluidRegions[i],
            dimensionedScalar(slowThermoFluid[i].p().dimensions()/dimTime, Zero)
        )
    );

    Info<< "    Adding to fieldsFluid\n" << endl;
    fieldsSlowFluid.set
    (
        i,
        new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
    );
    forAll(slowThermoFluid[i].composition().Y(), j)
    {
        fieldsSlowFluid[i].add(slowThermoFluid[i].composition().Y()[j]);
    }
    fieldsSlowFluid[i].add(slowThermoFluid[i].he());

    Info<< "    Adding to QdotFluid\n" << endl;
    QdotSlowFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                runTime.timeName(),
                slowFluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            slowFluidRegions[i],
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    const dictionary& pimpleDict =
        slowFluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowSlowFluid[i]);

    rhoMaxSlowFluid.set
    (
        i,
        new dimensionedScalar("rhoMax", dimDensity, GREAT, pimpleDict)
    );

    rhoMinSlowFluid.set
    (
        i,
        new dimensionedScalar("rhoMin", dimDensity, Zero, pimpleDict)
    );

    pressureControlsSlow.set
    (
        i,
        new pressureControl(slowThermoFluid[i].p(), rhoSlowFluid[i], pimpleDict, false)
    );

    Info<< "    Adding MRF\n" << endl;
    MRFSlowFluid.set
    (
        i,
        new IOMRFZoneList(slowFluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    slowFluidFvOptions.set
    (
        i,
        new fv::options(slowFluidRegions[i])
    );

    turbulenceSlowFluid[i].validate();

    pRefCellSlowFluid[i] = -1;
    pRefValueSlowFluid[i] = 0.0;

    if (p_rghSlowFluid[i].needReference())
    {
        setRefCell
        (
            slowThermoFluid[i].p(),
            p_rghSlowFluid[i],
            pimpleDict,
            pRefCellSlowFluid[i],
            pRefValueSlowFluid[i]
        );
    }
}
