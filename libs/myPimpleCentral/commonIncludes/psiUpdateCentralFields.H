//if( fluxScheme == "Kurganov" || fluxScheme == "Tadmor" )
{
    psif = linearInterpolate(psi);
    psi_own     = fvc::interpolate(psi, own, "reconstruct(psi)")*kappa
                  + (1.0 - kappa)*psif;
    psi_nei     = fvc::interpolate(psi, nei, "reconstruct(psi)")*kappa
                  + (1.0 - kappa)*psif;
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);
    if (noFluxWeighting)
    {
        phiv_own    = (fvc::interpolate(HbyA, own, "reconstruct(U)")) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(HbyA, nei, "reconstruct(U)")) & mesh.Sf();
    }
    else
    {
        psiU_own    = fvc::interpolate(psi*HbyA, own, "reconstruct(U)");
        psiU_nei    = fvc::interpolate(psi*HbyA, nei, "reconstruct(U)");
        phiv_own    = (psiU_own / psi_own) & mesh.Sf();
        phiv_nei    = (psiU_nei / psi_nei) & mesh.Sf();
    }
    if (mesh.changing())
    {
        phiv_own   -= mesh.phi();
        phiv_nei   -= mesh.phi();
    }
    
    MRF.makeRelative(phiv_own);
    MRF.makeRelative(phiv_nei);
    
    aphiv_own   = alpha_own*phiv_own - aSf;
    aphiv_nei   = alpha_nei*phiv_nei + aSf;
    
    phid_own    = aphiv_own * psi_own;
    phid_nei    = aphiv_nei * psi_nei;
    
    Dp_own = alpha_own*fvc::interpolate(rho*rAU, own, "reconstruct(Dp)");
    Dp_nei = alpha_nei*fvc::interpolate(rho*rAU, nei, "reconstruct(Dp)");
}
/*else if ( fluxScheme == "AUSM" )
{
    psif = linearInterpolate(psi);
    psi_own     = fvc::interpolate(psi, own, "reconstruct(psi)")*kappa
                  + (1.0 - kappa)*psif;
    psi_nei     = fvc::interpolate(psi, nei, "reconstruct(psi)")*kappa
                  + (1.0 - kappa)*psif;
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);
    
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
    
    surfaceScalarField Un_own ( phiv_own/mesh.magSf()); Un_own.setOriented(false);
    surfaceScalarField Un_nei ( phiv_nei/mesh.magSf()); Un_nei.setOriented(false);
    
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
    
    
}*/

