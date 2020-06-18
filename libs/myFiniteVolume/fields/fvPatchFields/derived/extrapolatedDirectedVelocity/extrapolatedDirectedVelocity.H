/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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
    Foam::extrapolatedDirectedVelocity

Description
    extrapolatedDirectedVelocity implements a mixed type boundary condition.
    It uses zero gradient condition for normal component of the velocity
    whereas the tangential components are calculated with given direction.

    i.e.

    \partial u_n / \partial n = 0
    U = inletDirection * (U.n) / (inletDirection.n)

\*---------------------------------------------------------------------------*/

#ifndef extrapolatedDirectedVelocityFvPatchVectorFields_H
#define extrapolatedDirectedVelocityFvPatchVectorFields_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class extrapolatedDirectedVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class extrapolatedDirectedVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // Private data

    //- Inlet direction
    vectorField inletDir_;

public:

    //- Runtime type information
    TypeName("extrapolatedDirectedVelocity");


    // Constructors

        //- Construct from patch and internal field
        extrapolatedDirectedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedDirectedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given extrapolatedDirectedVelocityFvPatchVectorField
        //  onto a new patch
        extrapolatedDirectedVelocityFvPatchVectorField
        (
            const extrapolatedDirectedVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        extrapolatedDirectedVelocityFvPatchVectorField
        (
            const extrapolatedDirectedVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new extrapolatedDirectedVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extrapolatedDirectedVelocityFvPatchVectorField
        (
            const extrapolatedDirectedVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new extrapolatedDirectedVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& inletDir() 
            {
                return inletDir_;
            }


            // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


            //- Evaluate the patch field
            virtual void evaluate
            (
	     
	     );

  

        //- Write
        virtual void write(Ostream&) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
