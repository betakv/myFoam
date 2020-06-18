    if( fluxScheme == "Kurganov" || fluxScheme == "Tadmor" )
    {
        if (noFluxWeighting)
        {
            phiv_own    = (fvc::interpolate(U, own, "reconstruct(U)")) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(U, nei, "reconstruct(U)")) & mesh.Sf();
        }
        else
        {
            rho_own     = fvc::interpolate(rho, own, "reconstruct(rho)");
            rho_nei     = fvc::interpolate(rho, nei, "reconstruct(rho)");
            correctCentralACMIInterpolation(rho_nei);
            phiv_own    = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        }
        
        if (mesh.changing())
        {
            phiv_own   -= mesh.phi();
            phiv_nei   -= mesh.phi();
        }
        
        c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
        cf_own      = fvc::interpolate(c, own, "reconstruct(psi)");
        cf_nei      = fvc::interpolate(c, nei, "reconstruct(psi)");
        correctCentralACMIInterpolation(cf_nei);
        cSf_own     = cf_own * mesh.magSf();
        cSf_own.setOriented(true);
        cSf_nei     = cf_nei * mesh.magSf();
        cSf_nei.setOriented(true);
        
        ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
        am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);
        
        if (fluxScheme == "Kurganov")
        {
            alpha_own   = ap/(ap - am);
            aSf         = am*alpha_own;
            alpha_nei   = 1.0 - alpha_own;
        }
        else if (fluxScheme == "Tadmor")
        {
            aSf         = -0.5*max(mag(ap), mag(am));
        }
        else
        {
            FatalErrorIn
            (
                    "updateCentralWeights.H"
            )   << "Unknown central flux scheme " << fluxScheme << nl << nl
                << "Valid flux schemes are:" << nl
                << "1) Kurganov" << nl
                << "2) Tamodr" << endl;
        }
    }
    else if ( fluxScheme == "AUSM" )
    {
        rho_own     = fvc::interpolate(rho, own, "reconstruct(rho)");
        rho_nei     = fvc::interpolate(rho, nei, "reconstruct(rho)");
        correctCentralACMIInterpolation(rho_nei);
        phiv_own    = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        
        if (mesh.changing())
        {
            phiv_own   -= mesh.phi();
            phiv_nei   -= mesh.phi();
        }
        
        volScalarField H ( 0.5 * magSqr(U) + h );
    
        H.max(SMALL);
        
        const volScalarField gamma( "gamma",  thermo.Cp()/thermo.Cv());
    
        const surfaceScalarField cCrit
        (  
            "cCrit", 
             fvc::interpolate
             (
                sqrt( 2 * H * ( gamma - 1 )/( gamma + 1) )
             )
        );
        
        
        surfaceScalarField Un_own ( phiv_own/mesh.magSf()); Un_own.setOriented(false);
        surfaceScalarField Un_nei ( phiv_nei/mesh.magSf()); Un_nei.setOriented(false);
        
        cf_own = sqr(cCrit) / max(cCrit,  Un_own);
        cf_nei = sqr(cCrit) / max(cCrit, -Un_nei );
        correctCentralACMIInterpolation(cf_nei);
        
        cSf_own = cf_own * mesh.magSf();
        cSf_own.setOriented(true);
        cSf_nei = cf_nei * mesh.magSf();
        cSf_nei.setOriented(true);
        
        ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
        am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);
        
        const surfaceScalarField rho12("rho12", 0.5 * ( rho_own + rho_nei ) );
        
        const surfaceScalarField c12("c12", min( cf_own, cf_nei) );
    
        
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
        
        const surfaceScalarField M12
        (
            "M12", 
            M4_pos + M4_neg 
        );
        
        const surfaceScalarField phiV( c12 * M12 * mesh.magSf());
        surfaceScalarField magPhiV = mag(phiV); magPhiV.max(SMALL);
        
        aSf = 0.5*(phiV *(mag(ap) + mag(am)) - magPhiV*(mag(ap) - mag(am)))/magPhiV;
        aSf.setOriented(true);
        
        
        surfaceScalarField alpha ("alha", 0.1875 * (-4 + 5 * sqr(fa))); 
        alpha.max(-0.75); 
        alpha.min(0.1875);
        
        alpha_own   = M_pos_one * 0.5*(1+sign(M_pos)) + (1.0 - M_pos_one) * M2_pos_p * (( 2.0-M_pos) - 16.0 * alpha * M_pos * M2_pos_n);
        alpha_nei   = M_neg_one * 0.5*(1-sign(M_neg)) + (1.0 - M_neg_one) * M2_neg_n * ((-2.0-M_neg) + 16.0 * alpha * M_neg * M2_neg_p);
        
    }
