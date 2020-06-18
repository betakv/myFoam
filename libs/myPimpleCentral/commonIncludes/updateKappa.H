if( fluxScheme == "Kurganov" || fluxScheme == "Tadmor" )
{
    aphiv_own = phi_own / rho_own;
    aphiv_nei = phi_nei / rho_nei;
    dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
    amaxSf = max(mag(aphiv_own), mag(aphiv_nei)) + amaxSmall;
    amaxSf.setOriented(true);
    
    kappaFuncPtr().update();
    kappa = kappaFuncPtr().kappa()();

    Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
    phi_own += (1.0 - kappa) * phi_nei;
    phi_nei = kappa * phi_nei;
}
else if ( fluxScheme == "AUSM" )
{
    aphiv_own = phi_own / rho_own;
    aphiv_nei = phi_nei / rho_nei;
    dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
    amaxSf = max(mag(aphiv_own), mag(aphiv_nei)) + amaxSmall;
    amaxSf.setOriented(true);
    
    kappaFuncPtr().update();
    kappa = kappaFuncPtr().kappa()();

    Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
    
    const surfaceScalarField rho12("rho12", 0.5 * ( rho_own + rho_nei ) );
        
    const surfaceScalarField c12("c12", min( cf_own, cf_nei) );
    
    surfaceScalarField Un_own ( phiv_own/mesh.magSf()); Un_own.setOriented(false);
    surfaceScalarField Un_nei ( phiv_nei/mesh.magSf()); Un_nei.setOriented(false);
    const surfaceScalarField M_pos("M_pos", Un_own / c12); 
    const surfaceScalarField M_neg("M_neg", Un_nei / c12);
        
    const surfaceScalarField M_avg("M_avg", sqrt( 0.5*( sqr(M_pos) + sqr(M_neg)+SMALL))); 
          
    surfaceScalarField M_o ("M_o",  min( 1.0, max( sqr(M_avg), sqr(MInf)))); 
    M_o.max(0.0);
    M_o.min(1.0);
    M_o=sqrt(M_o+SMALL);
        
    surfaceScalarField fa  ("fa", M_o*(2.0-M_o));//( "fa", sqrt( sqr(1.0-sqr(M_o))*sqr(M_avg)+4*sqr(M_o) ) / (1.0+sqr(M_o)));
    fa.max(0.0);
    fa.min(1.0); 
        
    surfaceScalarField M_pos_one ("M_pos_one", sign (mag(M_pos) - 1.0)); 
    M_pos_one.max(0.0);
    surfaceScalarField M_pos_one_low("M_pos_one_low", 1.0-M_pos_one); 
        
    surfaceScalarField M_neg_one ("M_neg_one", sign(mag(M_neg) - 1.0 )); 
    M_neg_one.max(0.0);
    surfaceScalarField M_neg_one_low("M_pos_one_low", 1.0-M_neg_one); 

    const surfaceScalarField mag_M_pos("mag_M_pos", mag(M_pos)); 
    const surfaceScalarField mag_M_neg("mag_M_pos", mag(M_neg)); 
    
    const surfaceScalarField M1_pos_p("M1_pos_p", 0.5 * (M_pos + mag_M_pos)); 
    const surfaceScalarField M1_pos_n("M1_pos_n", 0.5 * (M_pos - mag_M_pos)); 
    
    const surfaceScalarField M1_neg_p("M1_neg_p", 0.5 * (M_neg + mag_M_neg)); 
    const surfaceScalarField M1_neg_n("M1_neg_n", 0.5 * (M_neg - mag_M_neg)); 
    
    const surfaceScalarField M2_pos_p("M2_pos_p",  0.25 * sqr(M_pos + 1.0)); 
    const surfaceScalarField M2_pos_n("M2_pos_n", -0.25 * sqr(M_pos - 1.0)); 
        
    const surfaceScalarField M2_neg_p("M2_neg_p",  0.25 * sqr(M_neg + 1.0)); 
    const surfaceScalarField M2_neg_n("M2_neg_n", -0.25 * sqr(M_neg - 1.0));
    
    const surfaceScalarField M4_pos
    (
        "M4_pos", 
        M_pos_one * M1_pos_p 
      + M_pos_one_low * M2_pos_p * (1.0 - 2 * M2_pos_n)
    );
        
    const surfaceScalarField M4_neg
    (
        "M4_neg", 
        M_neg_one * M1_neg_n 
      + M_neg_one_low * M2_neg_n * (1.0 + 2 * M2_neg_p) 
    );
    
    phi_own = rho_own * cf_own * M4_pos * mesh.magSf();
    phi_nei = rho_nei * cf_nei * M4_neg * mesh.magSf();
}
