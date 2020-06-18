    fvMesh& mesh = slowFluidRegions[i];

    CombustionModel<rhoReactionThermo>& reaction = reactionSlowFluid[i];

    rhoReactionThermo& thermo = reaction.thermo();
    thermo.validate(args.executable(), "h", "e");

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    label inertIndex = -1;
    if (Y.size())
    {
        const word inertSpecie(thermo.get<word>("inertSpecie"));
        if (!composition.species().found(inertSpecie))
        {
            FatalIOErrorIn(args.executable().c_str(), thermo)
                << "Inert specie " << inertSpecie
                << " not found in available species "
                << composition.species()
                << exit(FatalIOError);
        }
        inertIndex = composition.species()[inertSpecie];
    }

    volScalarField& rho = rhoSlowFluid[i];
    volVectorField& U = USlowFluid[i];
    surfaceScalarField& phi = phiSlowFluid[i];

    compressible::turbulenceModel& turbulence = turbulenceSlowFluid[i];
    volScalarField& K = KSlowFluid[i];
    volScalarField& dpdt = dpdtSlowFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghSlowFluid[i];

    const volScalarField& gh = ghSlowFluid[i];
    const surfaceScalarField& ghf = ghfSlowFluid[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
        fieldsSlowFluid[i];

    volScalarField& Qdot = QdotSlowFluid[i];

    radiation::radiationModel& rad = radiationSlow[i];

    IOMRFZoneList& MRF = MRFSlowFluid[i];
    fv::options& fvOptions = slowFluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassSlowFluid[i]
    );

    bool frozenFlow = frozenFlowSlowFluid[i];

    const label pRefCell = pRefCellSlowFluid[i];
    const scalar pRefValue = pRefValueSlowFluid[i];

    const dimensionedScalar rhoMax = rhoMaxSlowFluid[i];
    const dimensionedScalar rhoMin = rhoMinSlowFluid[i];

    const pressureControl& pressureControl = pressureControlsSlow[i];

