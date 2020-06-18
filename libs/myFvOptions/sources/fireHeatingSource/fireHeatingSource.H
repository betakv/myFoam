/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Class
    Foam::fv::fireHeatingSource

Group
    grpFvOptionsSources

Description
    


Usage
    

SourceFiles
    fireHeatingSource.C

SeeAlso
    Foam::Function1

\*---------------------------------------------------------------------------*/

#ifndef fv_fireHeatingSource_H
#define fv_fireHeatingSource_H

#include "fvOption.H"
#include "Function1.H"
#include "coordinateSystem.H"
#include "cellSetOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class fireHeatingSource Declaration
\*---------------------------------------------------------------------------*/

class fireHeatingSource
:
    public cellSetOption
    //public option
{
    // Private data

        //- Name of electrical conductivity field
        static const word sigmaName;

        //- Name of temperature field - default = "T" (optional)
        word TName_;

        // Total heat flux in volume
        scalar Q_; 
        
        // Distribution of fire heat source in time Q_ =c1 * sum(q[i]*t^i)
        List<scalar> q_;
        
        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- No copy construct
        fireHeatingSource(const fireHeatingSource&) = delete;

        //- No copy assignment
        void operator=(const fireHeatingSource&) = delete;

        


public:

    //- Runtime type information
    TypeName("fireHeatingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        fireHeatingSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~fireHeatingSource();


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fireHeatingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
