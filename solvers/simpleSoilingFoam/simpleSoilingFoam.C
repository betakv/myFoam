/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        const volVectorField VsubU (V-U);
        
        // computing turbulent viscosity for droples concentration
        
        const volScalarField gamma
        ( 
            mag(VsubU)
          / sqrt(2./3.*turbulence->k())
        ); 
        
        const volScalarField dropletTurbVis
        (
            "dropletTurbVis",
            turbulence->nuEff()/Sct
          * sqrt(gamma*gamma+1.0)/(gamma*gamma+1.0)
        );
        
        const volSymmTensorField additViscosity
        (
            "additViscosity",
            delta
          - 0.5 *  turbulence->nuEff()/Sct * gamma/(gamma*gamma+1.0) * I
        );
        
        // solving concentration equation
        
        fvScalarMatrix Ceqn
        (
            fvm::div( phiV, C)
          - fvm::laplacian(dropletTurbVis, C)
          ==
            fvc::div(additViscosity & fvc::grad(C))
        );
        
        Ceqn.relax();
        Ceqn.solve();
        
        // bounding negative values
        C.max(0.0);
        
        // computing particle forces
        
        //volScalarField Re(mag(VsubU) * d / turbulence->nuEff());
        volScalarField Re(mag(VsubU) * d / laminarTransport.nu());
        
        Re.max(.1);
        Re.min(1e6);
        
        volScalarField cd
        (
            24/Re
          + 2.6 * (Re / 5) / (1 +pow(Re/5,1.52))
          + 0.411 * pow (Re / 2.63e5,-7.94) / (1 + pow (Re / 2.63e5,-8.0))
          + 0.25 *  (Re / 1e6) / (1 + (Re / 1e6))
        );
       /* 
        volScalarField cd = Re;
        
        forAll(cd,i)
        {
            scalar phi_1 = pow( 0.4 , 10) + pow(24.0 / Re[i], 10)
                                + pow(21.0 / ::pow(Re[i],0.67), 10)
                                + pow( 4.0 / ::pow(Re[i],0.33), 10);
             scalar phi_2 = 1.0 / 
        						(  ::pow(0.148*::pow(Re[i],0.11),-10.0)
        						   + ::pow(0.5,-10)
        						);
        scalar phi_3 = pow(1.57e8/::pow(Re[i],1.625),10);
        scalar phi_4 = 1.0 / 
        						(  ::pow(1.67e-17*::pow(Re[i],2.63),-10.0)
        						   + ::pow(0.2,-10)
        						); 
        cd[i]=::pow(1.0/(1.0/(phi_1+phi_2)+1.0/phi_3)
                   +phi_4,0.01);						                    
        
    }
    */
        
        Info << "Cd min/max " << (min(cd)).value() << " " << (max(cd)).value() << endl;

        const dimensionedScalar rho = rhoAir / rhoSoil;
        const dimensionedScalar dragConst = -0.75 * rho / d;
        const dimensionedVector liftForce = ( 1 - rho ) * g;
        
        const volScalarField dragForce = dragConst * cd * mag(VsubU);

        fvVectorMatrix Veqn
        (
            fvm::div(phiV, V)
          ==
            liftForce
          + fvm::SuSp( dragForce, V )
          - dragForce * U
          +  fvOptions(V)
        );

        Veqn.relax();
        Veqn.solve();
        
        fvOptions.correct(V);

        // correction of boundary conditions of droplets fields and drolets flux 
        
        V.correctBoundaryConditions();
        C.correctBoundaryConditions();
        
        phiV =  linearInterpolate(V) & mesh.Sf();

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }
        
        /*dimensionedScalar magU(max(mag(U)));
        
        forAlll(V,celli)
        {
            if(V[celli] > magU.value())
            {
            
            }
        }
*/
        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
