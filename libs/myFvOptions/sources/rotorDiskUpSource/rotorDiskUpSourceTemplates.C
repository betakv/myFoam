/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "rotorDiskUpSource.H"
#include "volFields.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::rotorDiskUpSource::calculate
(
    const RhoFieldType& rho,
    const volScalarField& p,
    const vectorField& U,
    const volScalarField& nu,
    vectorField& force,
    const bool divideVolume,
    const bool output
) const
{
    const scalarField& V = mesh_.V();

    // Logging info
    scalar torque = 0.0;
    scalar thrust = 0.0;
    
    scalar massFlow = 0.0;
    
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;

    // Cached position-dependent rotations available?
    const bool hasCache = Rcyl_.valid();

    List<vector> UcBEM(nR_*nTheta_,Zero);
    List<scalar> rhoBEM(nR_*nTheta_,Zero);
    List<scalar> pBEM(nR_*nTheta_,Zero);
    List<scalar> nuBEM(nR_*nTheta_,Zero);
    List<vector> forceBEM(nR_*nTheta_,Zero);

    forAll(cells_, i)
    {
        if (x_[i].x() >= rMin_ && x_[i].x() <= rMax_)
        {
            const label celli = cells_[i];

            const scalar radius = x_[i].x();
            const tensor Rcyl =
            (
                hasCache
              ? (*Rcyl_)[i]
              : coordSys_.R(mesh_.C()[celli])
            );
            
            // Transform velocity into local cylindrical reference frame
            vector Uc = invTransform(Rcyl, U[celli]);

            // Transform velocity into local coning system
            Uc = transform(Rcone_[i], Uc);

            // Set radial component of velocity to zero
            Uc.x() = 0.0;

            // Set blade normal component of velocity
            Uc.y() = radius*omega_ - Uc.y();
            
            // Interpolate velocity on BEM disk
            for( label j = 0; j < nTheta_; j++)
            {
                for( label k = 0; k < nR_; k++)
                {
                    label interpolationCelli = i*nR_*nTheta_ + j*nR_ + k;
                    label BEMcelli = j*nR_ + k;
                    
                    scalar weight = interpolationWeights_[interpolationCelli];
                    
                    UcBEM[BEMcelli] += weight * Uc;
                    rhoBEM[BEMcelli] += weight * rho[celli];
                    pBEM[BEMcelli] += weight * p[celli];
                    nuBEM[BEMcelli] += weight * nu[celli];
                }
            }
        }
    }
    
    for( label i = 0; i < nTheta_; i++)
    {
        for( label j = 0; j < nR_; j++)
        {
            label BEMcelli = i*nR_ + j;
                    
            UcBEM[BEMcelli] *= 1.0/sumInterpolationWeights_[BEMcelli];
            rhoBEM[BEMcelli] *= 1.0/sumInterpolationWeights_[BEMcelli];
            pBEM[BEMcelli] *= 1.0/sumInterpolationWeights_[BEMcelli];
            nuBEM[BEMcelli] *= 1.0/sumInterpolationWeights_[BEMcelli];
        }
    }
    
    reduce(UcBEM, sumOp<List<vector>>());
    reduce(rhoBEM, sumOp<List<scalar>>());
    reduce(nuBEM, sumOp<List<scalar>>());
    
        
    // Calculate forces on BEM disk
     
    for( label i = 0; i < nTheta_; i++)
    {
        for( label j = 0; j < nR_; j++)
        {
            label BEMcelli = i*nR_ + j;
                
            scalar radius = BEMcoordinates_[BEMcelli].x();
            vector Uc = UcBEM[BEMcelli];
            
            massFlow += UcBEM[BEMcelli].z() * rhoBEM[BEMcelli] * area2_[BEMcelli];
            
            //Info << radius << endl;
                
            scalar twist = 0.0;
            scalar chord = 0.0;
            label i1 = -1;
            label i2 = -1;
            scalar invDr = 0.0;
            blade_.interpolate(radius, twist, chord, i1, i2, invDr);
            
            //Info << twist << " " << chord << " " << i1 << " " << i2 << endl;
                
            // Flip geometric angle if blade is spinning in reverse (clockwise)
            scalar alphaGeom = + twist;
            if (omega_ < 0)
            {
                alphaGeom = mathematical::pi - alphaGeom;
            }
                
            // Effective angle of attack
            scalar alphaEff = alphaGeom - atan2(Uc.z(), Uc.y());
            if (alphaEff > mathematical::pi)
            {
                alphaEff -= mathematical::twoPi;
            }
            if (alphaEff < -mathematical::pi)
            {
                alphaEff += mathematical::twoPi;
            }
                
            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            // Determine profile data for this radius and angle of attack
            const label profile1 = blade_.profileID()[i1];
            const label profile2 = blade_.profileID()[i2];
            
            scalar Re = mag(Uc)*chord/nuBEM[BEMcelli];
            
            scalar Ma = 0;
            
            if( gamma_ == 1 )
            {
                Ma = mag(Uc) / aRef_;
            }
            else
            {
                Ma = mag(Uc) / sqrt (gamma_*rhoBEM[BEMcelli]*pBEM[BEMcelli]);
            }
            
                 
            scalar Cd1 = 0.0;
            scalar Cl1 = 0.0;
           
            if (ReMa_ == "Re" )
            {
                profiles_[profile1].Cdl(Re, alphaEff, Cd1, Cl1);
            }
            else if (ReMa_ == "Ma" )
            {
                profiles_[profile1].Cdl(Ma, alphaEff, Cd1, Cl1);
            }
            else
            {
                Info << "Wrong interpolation type of characteristic curve" << endl;
            }
           
           
            
           //  Info << Re << " " << Cd1 << " " << Cl1 << endl;   
            scalar Cd2 = 0.0;
            scalar Cl2 = 0.0;
            if (ReMa_ == "Re" )
            {
                profiles_[profile2].Cdl(Re, alphaEff, Cd1, Cl1);
            }
            else if (ReMa_ == "Ma" )
            {
                profiles_[profile2].Cdl(Ma, alphaEff, Cd1, Cl1);
            }
            else
            {
                Info << "Wrong interpolation type of characteristic curve" << endl;
            }
                
            scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
            scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

            // Calculate forces perpendicular to blade
            scalar pDyn = 0.5*rho[BEMcelli]*magSqr(Uc);
                
            scalar f = pDyn*chord*nBlades_*dR_*area_[BEMcelli];
            
            if( mode_ == "propeller" )
            {
                forceBEM[BEMcelli] = vector(0.0, -f*(Cd*cos(alphaEff)+Cl*sin(alphaEff)), -f*(Cl*cos(alphaEff)-Cd*sin(alphaEff)));
            }
            else if( mode_ == "turbine" )
            {
                forceBEM[BEMcelli] = vector(0.0, f*(Cd*cos(alphaEff)+Cl*sin(alphaEff)), f*(Cl*cos(alphaEff)-Cd*sin(alphaEff)));
            }
            else
            {
                forceBEM[BEMcelli] = vector(0.0, 0.0, 0.0);
            }    
            
        }
    }
        
    vector totalForce = sum(forceBEM);
    
    Info << totalForce << endl;
        
    forAll(cells_, i)
    {
        const label celli = cells_[i];
            
        vector localForce(Zero);
            
        if (x_[i].x() >= rMin_ && x_[i].x() <= tipFactor_*rMax_)
        {
            for( label j = 0; j < nTheta_; j++)
            {
                for( label k = 0; k < nR_; k++)
                {
                    label interpolationCelli = i*nR_*nTheta_ + j*nR_ + k;
                    label BEMcelli = j*nR_ + k;
                    
                    if(backInterpolationWeights_[interpolationCelli] != 0)
                    {
                       localForce += backInterpolationWeights_[interpolationCelli] / sumBackInterpolationWeights_[BEMcelli]
                                  * forceBEM[BEMcelli];
                    }
               }
            }
        }    
            
        // Accumulate forces
        torque += rhoRef_*localForce.y()*x_[i].x();
        thrust += rhoRef_*localForce.z();
            
        const tensor Rcyl =
        (
            hasCache
         ? (*Rcyl_)[i]
         : coordSys_.R(mesh_.C()[celli])
        );


        // Transform force from local coning system into rotor cylindrical
        localForce = invTransform(Rcone_[i], localForce);

        // Transform force into global Cartesian coordinate system
        force[celli] = transform(Rcyl, localForce);
        
        force[celli] /= V[celli];

         //   if (divideVolume)
         //   {
         //       
         //   }
            
            
        }
          
    


        reduce(AOAmin, minOp<scalar>());
        reduce(AOAmax, maxOp<scalar>());
        reduce(torque, sumOp<scalar>());
        reduce(thrust, sumOp<scalar>());

        Info<< type() << " output:" << nl
            << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
            << radToDeg(AOAmax) << nl
            << "    Torque =    " << -1*torque << " Nm " << nl
            << "    Thrust =    " << -1*thrust << " N " << nl
            << "    Mass flow = " << massFlow << " kg/s" << endl;
    
}


template<class Type>
void Foam::fv::rotorDiskUpSource::writeField
(
    const word& name,
    const List<Type>& values,
    const bool writeNow
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    if (mesh_.time().writeTime() || writeNow)
    {
        auto tfield = tmp<FieldType>::New
        (
            IOobject
            (
                name,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<Type>(dimless, Zero)
        );

        auto& field = tfield.ref().primitiveFieldRef();

        if (cells_.size() != values.size())
        {
            FatalErrorInFunction
                << abort(FatalError);
        }

        forAll(cells_, i)
        {
            const label celli = cells_[i];
            field[celli] = values[i];
        }

        tfield().write();
    }
}


// ************************************************************************* //
