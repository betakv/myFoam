/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "lookupUpProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "vector.H"
#include "symmTensor.H"
#include "unitConversion.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lookupUpProfile, 0);
    addToRunTimeSelectionTable(profileModelUp, lookupUpProfile, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::lookupUpProfile::interpolateWeights
(
    const scalar& xIn,
    const List<scalar>& values,
    label& i1,
    label& i2,
    scalar& ddx
) const
{
   
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookupUpProfile::lookupUpProfile
(
    const dictionary& dict,
    const word& modelName
)
:
    profileModelUp(dict, modelName),
    Re_(),
    AOA_(),
    Cd_(),
    Cl_(),
    interpolationRe_(""),
    interpolationAlpha_("")
{
    List<symmTensor> data;
    if (readFromFile())
    {
        IFstream is(fName_);
        is  >> data;
    }
    else
    {
        dict.readEntry("data", data);
    }

    if (data.size())
    {
        Re_.setSize(data.size());
        AOA_.setSize(data.size());
        Cd_.setSize(data.size());
        Cl_.setSize(data.size());

        forAll(data, i)
        {
            Re_[i] = data[i][0];
            AOA_[i] = degToRad(data[i][1]);
            Cd_[i] = data[i][2];
            Cl_[i] = data[i][3];
        }
    }
    else
    {
        FatalErrorInFunction
            << "No profile data specified" << exit(FatalError);
    }
    
    interpolationRe_ = dict.get<word>("interoplationRe");
    interpolationAlpha_ = dict.get<word>("interoplationAlpha");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lookupUpProfile::Cdl(const scalar Re, const scalar alpha, scalar& Cd, scalar& Cl) const
{
    scalar ReL = -1;
    scalar ReH = -1;
    
    for( label i = 0; i < Re_.size(); i++ )
    {
        if(Re >= Re_[i]) 
        {
            ReL = Re_[i];
        }
        else
        {
            break;
        }
    }
    
    for( label i = Re_.size() -1; i >= 0 ; i++ )
    {
        if(Re <= Re_[i]) 
        {
            ReH = Re_[i];
        }
        else
        {
            break;
        }
    }
    
    if( ReL == -1 ) ReL = ReH;
    if( ReH == -1 ) ReH = ReL;
    
    scalar clL = 0, cdL = 0;
    
    for( label i = 0; i < Re_.size() - 1; i++ )
    {
        if( Re_[i] == ReL )
        {
            if(alpha >= AOA_[i] && alpha <=AOA_[i+1])
            {
                if( interpolationAlpha_ == "weighted" )
                {
                    scalar weightL = 1.0 / max(sqr(AOA_[i] - alpha), SMALL);
                    scalar weightR = 1.0 / max(sqr(AOA_[i+1] - alpha), SMALL);
                
                    cdL = ( weightL*Cd_[i] + weightR*Cd_[i+1] )/( weightL + weightR );
                    clL = ( weightL*Cl_[i] + weightR*Cl_[i+1] )/( weightL + weightR );
                }
                else if(  interpolationAlpha_ == "linear" )
                {
                    cdL = Cd_[i] + (Cd_[i+1]-Cd_[i])/(AOA_[i+1] - AOA_[i]) * alpha;
                    clL = Cl_[i] + (Cl_[i+1]-Cl_[i])/(AOA_[i+1] - AOA_[i]) * alpha;
                }
                else
                {
                    cdL = 0.0;
                    clL = 0.0;
                }    
            }
        }
    }
    
    scalar clH = 0, cdH = 0;
    
    for( label i = 0; i < Re_.size() - 1; i++ )
    {
        if( Re_[i] == ReH )
        {
            if(alpha >= AOA_[i] && alpha <=AOA_[i+1])
            {
                if( interpolationAlpha_ == "weighted" )
                {
                    scalar weightL = 1.0 / max(sqr(AOA_[i]- alpha), SMALL);
                    scalar weightR = 1.0 / max(sqr(AOA_[i+1]- alpha), SMALL);
                
                    cdH = ( weightL*Cd_[i] + weightR*Cd_[i+1] )/( weightL + weightR );
                    clH = ( weightL*Cl_[i] + weightR*Cl_[i+1] )/( weightL + weightR );
                }
                else if(  interpolationAlpha_ == "linear" )
                {
                    cdH = Cd_[i] + (Cd_[i+1]-Cd_[i])/(AOA_[i+1] - AOA_[i]) * alpha;
                    clH = Cl_[i] + (Cl_[i+1]-Cl_[i])/(AOA_[i+1] - AOA_[i]) * alpha;
                }
                else
                {
                    cdH = 0.0;
                    clH = 0.0;
                }    
            }
        }
    }
    
    if( interpolationRe_ == "weighted" )
    {
        scalar weightReL = 1.0 / max(sqr(ReL - Re), SMALL);
        scalar weightReH = 1.0 / max(sqr(ReH - Re), SMALL);
    
        Cd = ( weightReL*cdL + weightReH*cdH)/( weightReL + weightReH );
        Cl = ( weightReL*clL + weightReH*clH)/( weightReL + weightReH );
    }
    else if(  interpolationRe_ == "linear" )
    {
        Cd = cdL + (cdH-cdL)/(ReH-ReL) * Re;
        Cl = clL + (clH-clL)/(ReH-ReL) * Re;
    }
    else
    {
        Cd = 0.0;
        Cl = 0.0;
    }    
}


// ************************************************************************* //
