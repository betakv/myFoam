/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "IRSSmooth.H"
#include "fvCFD.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(IRSSmooth, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IRSSmooth::IRSSmooth( const fvMesh& mesh )
:
    mesh_(mesh),
    errMax_(mesh_.solutionDict().subDict("RK").lookupOrDefault<scalar>("errMax", 0.01)),
    epsilon_(mesh_.solutionDict().subDict("RK").lookupOrDefault<scalar>("epsilon", 0.75)),
    iterations_(mesh_.solutionDict().subDict("RK").lookupOrDefault<label>("iterations", 100)),
    IRSSwitch_(mesh_.solutionDict().subDict("RK").lookupOrDefault("IRS", false))

   {
    
   }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IRSSmooth::smoothScalar(volScalarField& u)
{
    if(IRSSwitch_)
    {
        volScalarField u_o  = u;
        volScalarField u_oo = u;
    
        volScalarField magGradU( mag(fvc::grad(u) ) );

        surfaceScalarField psi( fvc::interpolate( magGradU) );
        volScalarField sumPsi(fvc::surfaceSum(psi)); 
        dimensionedScalar sumPsiMin = min(sumPsi); if( sumPsiMin.value() == 0 ) sumPsiMin.value() = 1e-10;

        label k       = 0;
        scalar err     = 1.0;
    
        // Loop on k-th Jacobi sub-iteration 
        while ( ( k < iterations_ ) && ( err > errMax_ ) )
        {    
            surfaceScalarField u_o_s( fvc::interpolate(u_o) );
            u_oo = u + epsilon_ * fvc::surfaceSum (u_o_s * psi )/( sumPsi + sumPsiMin );
        
            forAll( mesh_.V(), i )
            {
                u_oo[i] = u_oo[i]/( 1 + epsilon_ * mesh_.cells()[i].nFaces() );   
            }
            u_oo.correctBoundaryConditions();
            
            // Update error
            scalarField cit = mag( u_oo - u_o );
            scalarField jmen = mag( u_o );
            forAll(jmen, i)
            {
                if(jmen[i] == 0 ) jmen[i] = 1e-10;
            }
        
            err = gSum( cit/jmen*mesh_.V() )/gSum( mesh_.V() );   
        
            // Update arrays the next (k + 1)-th Jacobi sub-iteration
            u_o = u_oo; 
            k   = k + 1;
        }
        u = u_oo;
    }    
}        
void Foam::IRSSmooth::smoothVector(volVectorField& u)
{
    if(IRSSwitch_)
    {
        volVectorField u_o  = u;
        volVectorField u_oo = u;
    
        volScalarField magGradU( mag(fvc::grad(u) ) );

        surfaceScalarField psi( fvc::interpolate( magGradU) );
        volScalarField sumPsi(fvc::surfaceSum(psi)); 
        dimensionedScalar sumPsiMin = min(sumPsi); if( sumPsiMin.value() == 0 ) sumPsiMin.value() = 1e-10;

        label k       = 0;
        scalar err     = 1.0;
    
        // Loop on k-th Jacobi sub-iteration 
        while ( ( k < iterations_ ) && ( err > errMax_ ) )
        {    
            surfaceVectorField u_o_s( fvc::interpolate(u_o) );
            u_oo = u + epsilon_ * fvc::surfaceSum (u_o_s * psi )/( sumPsi + sumPsiMin );
        
            forAll( mesh_.V(), i )
            {
                u_oo[i] = u_oo[i]/( 1 + epsilon_ * mesh_.cells()[i].nFaces() );   
            }
            u_oo.correctBoundaryConditions();
            
            // Update error
            scalarField cit = mag( u_oo - u_o );
            scalarField jmen = mag( u_o );
            forAll(jmen, i)
            {
                if(jmen[i] == 0 ) jmen[i] = 1e-10;
            }
        
            err = gSum( cit/jmen*mesh_.V() )/gSum( mesh_.V() );   
        
            // Update arrays the next (k + 1)-th Jacobi sub-iteration
            u_o = u_oo; 
            k   = k + 1;
        }
        u = u_oo;
    }   
}


// ************************************************************************* //
