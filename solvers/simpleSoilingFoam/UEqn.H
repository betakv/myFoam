    // Momentum predictor

    MRF.correctBoundaryVelocity(U);
    const scalarField& volume = mesh.V();
    
    volScalarField massFraction(0*C);
    
    forAll(C, celli)
    {
        scalar massSoil = C[celli] * volume[celli];
        
        scalar soilVolume = massSoil / rhoSoil.value();
        
        soilVolume = min(soilVolume, volume[celli]);
        
        scalar airVolume = volume[celli] - soilVolume;
        
        airVolume = max(airVolume, SMALL);
        
        scalar massAir = airVolume * rhoAir.value();
        
        massFraction[celli] = massSoil / massAir;
        
    }

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
      -  twoWayCoupling * massFraction * liftForce
     // -  twoWayCoupling * massFraction * dragForce * VsubU
      +  twoWayCoupling * massFraction * fvm::SuSp( dragForce, U )
      -  twoWayCoupling * massFraction * dragForce * V
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
