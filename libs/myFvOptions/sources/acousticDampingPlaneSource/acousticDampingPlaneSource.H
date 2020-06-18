/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fv::acousticDampingPlaneSource

Group
    grpFvOptionsSources

Description
    Acoustic damping source

    \heading Source usage

    Example usage:
    \verbatim
    acousticDampingPlaneSourceCoeffs
    {
    }
    \endverbatim

SourceFiles
    acousticDampingPlaneSource.C

\*---------------------------------------------------------------------------*/

#ifndef acousticDampingPlaneSource_H
#define acousticDampingPlaneSource_H

#include "cellSetOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
                    Class acousticDampingPlaneSource Declaration
\*---------------------------------------------------------------------------*/

class acousticDampingPlaneSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Frequency [Hz]
        dimensionedScalar frequency_;

        //- Blending factor []
        volScalarField blendFactor_;

        //- Name of reference velocity field
        word URefName_;

        // Centre location of damping planes
        List<vector> x0_;
        
        // Axis of damping planes

        List<vector> axis_;
        
        // Width of damping planes
        
        List<scalar> h_;

        //- Stencil width, default = 20
        label w_;


    // Protected Member Functions

        //- Helper function to set the blending factor
        void setBlendingFactor();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        acousticDampingPlaneSource(const acousticDampingPlaneSource&);

        //- Disallow default bitwise assignment
        void operator=(const acousticDampingPlaneSource&) = delete;


public:

    //- Runtime type information
    TypeName("acousticDampingPlaneSource");


    // Constructors

        //- Construct from components
        acousticDampingPlaneSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~acousticDampingPlaneSource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
