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

Class
    Foam::IRSSmooth

Description
    MRF zone definition based on cell zone and parameters
    obtained from a control dictionary constructed from the given stream.

    The rotation of the MRF region is defined by an origin and axis of
    rotation and an angular speed.

SourceFiles
    IRSSmooth.C

\*---------------------------------------------------------------------------*/

#ifndef IRSSmooth_H
#define IRSSmooth_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
//#include "fvMatricesFwd.H"
//#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class IRSSmooth Declaration
\*---------------------------------------------------------------------------*/

class IRSSmooth
{
    // Private data

        const fvMesh& mesh_;
        const scalar errMax_;
        const scalar epsilon_;
        const label iterations_;
        const bool IRSSwitch_;
        
        // Private Member Functions
        
        //- Disallow default bitwise copy construct
        IRSSmooth(const IRSSmooth&);

        //- Disallow default bitwise assignment
        void operator=(const IRSSmooth&);


public:

    // Declare name of the class and its debug switch
    ClassName("IRSSmooth");


    // Constructors

        //- Construct from fvMesh
        IRSSmooth( const fvMesh& mesh );

    // Member Functions

        //- Smooth scalar field
        void smoothScalar(volScalarField& u);
        
        //- Smooth vector field
        void smoothVector(volVectorField& u);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
