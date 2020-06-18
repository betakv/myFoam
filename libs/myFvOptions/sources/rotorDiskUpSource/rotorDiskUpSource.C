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
#include "addToRunTimeSelectionTable.H"
//#include "trimModelUp.H"
#include "pimpleControl.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "fvmSup.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(rotorDiskUpSource, 0);
        addToRunTimeSelectionTable(option, rotorDiskUpSource, dictionary);
    }
}


const Foam::Enum
<
    Foam::fv::rotorDiskUpSource::geometryModeType
>
Foam::fv::rotorDiskUpSource::geometryModeTypeNames_
({
    { geometryModeType::gmAuto, "auto" },
    { geometryModeType::gmSpecified, "specified" },
});


const Foam::Enum
<
    Foam::fv::rotorDiskUpSource::inletFlowType
>
Foam::fv::rotorDiskUpSource::inletFlowTypeNames_
({
    { inletFlowType::ifFixed, "fixed" },
    { inletFlowType::ifSurfaceNormal, "surfaceNormal" },
    { inletFlowType::ifLocal, "local" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::rotorDiskUpSource::checkData()
{

}


void Foam::fv::rotorDiskUpSource::setFaceArea(vector& axis, const bool correct)
{

}


void Foam::fv::rotorDiskUpSource::createCoordinateSystem()
{
    // Construct the local rotor coordinate system
    vector origin(Zero);
    vector axis(Zero);
    vector refDir(Zero);
    
    coeffs_.readEntry("origin", origin_);
    coeffs_.readEntry("axis", axis_);
    coeffs_.readEntry("refDirection", refDir_);

    coordSys_ = coordSystem::cylindrical(origin_, axis_, refDir_);

    Info<< "    Rotor gometry:" << nl
        << "    - origin        = " << coordSys_.origin() << nl
        << "    - r-axis        = " << coordSys_.e1() << nl
        << "    - psi-axis      = " << coordSys_.e2() << nl
        << "    - z-axis        = " << coordSys_.e3() << endl;
}


void Foam::fv::rotorDiskUpSource::constructGeometry()
{
    const pointUIndList cc(mesh_.C(), cells_);

    // Optional: for later transform(), invTransform()
    /// Rcyl_.reset(coordSys_.R(cc).ptr());

    forAll(cells_, i)
    {
        //if (area_[i] > ROOTVSMALL)
        {
            // Position in (planar) rotor coordinate system
            x_[i] = coordSys_.localPosition(cc[i]);

            // Swept angle relative to rDir axis [radians] in range 0 -> 2*pi
            scalar psi = x_[i].y();

            // Blade flap angle [radians]
            //scalar beta =
            //    flap_.beta0 - flap_.beta1c*cos(psi) - flap_.beta2s*sin(psi);

            // Determine rotation tensor to convert from planar system into the
            // rotor cone system
            //scalar c = cos(beta);
            //scalar s = sin(beta);
            //Rcone_[i] = tensor(c, 0, -s, 0, 1, 0, s, 0, c);
        }
    }
}

void Foam::fv::rotorDiskUpSource::updateGeometry()
{
    const pointField& points = mesh_.points();
    const pointField& oldPoints = mesh_.oldPoints();
    
    scalarField distance = mag(points-oldPoints);
    vector point1(Zero), point2(Zero), point3(Zero), point4(Zero);
    vector oldPoint1(Zero), oldPoint2(Zero), oldPoint3(Zero), oldPoint4(Zero);
    List<vector> centers; centers.setSize(4,Zero);
    List<vector> oldCenters; oldCenters.setSize(4,Zero);
    
    for( label i = 0; i < 4; i++)
    {    
        const label celli = cells_[i];
        labelList cellPoints = mesh_.cellPoints()[celli];
        forAll( cellPoints, cellPointi)
        {
            centers[i] += points[cellPoints[cellPointi]];
            oldCenters[i] += oldPoints[cellPoints[cellPointi]];
        }
        
        centers[i] *= 1./cellPoints.size();
        oldCenters[i] *= 1./cellPoints.size();
    
    }
    
    scalar magDistance(SMALL);
    
    label nChanges = 0;
  
    tensor rotation(Zero);
    vector translation(Zero);
    
    if(nChanges >= 4)
    {
        scalarSquareMatrix ATA(4, scalar(1));
        List<scalar> ATb(4, scalar(0));
        
        // create matrix
        
        ATA(0,1) = oldCenters[0].x();
        ATA(0,2) = oldCenters[0].y();
        ATA(0,3) = oldCenters[0].z();
        
        ATA(1,1) = oldCenters[1].x();
        ATA(1,2) = oldCenters[1].y();
        ATA(1,3) = oldCenters[1].z();
        
        ATA(2,1) = oldCenters[2].x();
        ATA(2,2) = oldCenters[2].y();
        ATA(2,3) = oldCenters[2].z();
        
        ATA(3,1) = oldCenters[3].x();
        ATA(3,2) = oldCenters[3].y();
        ATA(3,3) = oldCenters[3].z();
        
        //create right hand side for x
        
        ATb[0] = centers[0].x();
        ATb[1] = centers[1].x();
        ATb[2] = centers[2].x();
        ATb[3] = centers[3].x();
        
        LUsolve(ATA, ATb);
        
        translation.x() = ATb[0];
        rotation.xx() = ATb[1];
        rotation.xy() = ATb[2];
        rotation.xz() = ATb[3];
        
        //create right hand side for y
        
        ATb[0] = centers[0].y();
        ATb[1] = centers[1].y();
        ATb[2] = centers[2].y();
        ATb[3] = centers[3].y();
        
        LUsolve(ATA, ATb);
        
        translation.y() = ATb[0];
        rotation.yx() = ATb[1];
        rotation.yy() = ATb[2];
        rotation.yz() = ATb[3];
        
        //create right hand side for z
        
        ATb[0] = centers[0].z();
        ATb[1] = centers[1].z();
        ATb[2] = centers[2].z();
        ATb[3] = centers[3].z();
        
        LUsolve(ATA, ATb);
        
        translation.z() = ATb[0];
        rotation.zx() = ATb[1];
        rotation.zy() = ATb[2];
        rotation.zz() = ATb[3];
    }
    
    reduce(rotation, maxMagSqrOp<tensor>());
    reduce(translation, maxMagSqrOp<vector>());
    
    axis_ = rotation & axis_;
    refDir_ = rotation & refDir_;
    
    origin_ =  rotation & origin_;
    origin_ += translation;
    
    // actualization of local coordinate system
    
    coordSys_ = coordSystem::cylindrical(origin_, axis_, refDir_);
    
     for( label i = 0; i < nTheta_; i++)
        {
            scalar theta = i * dTheta_;
            
            for( label j = 0; j < nR_; j++)
            {
                label celli = i*nR_ + j;
                
                BEMcoordinates_[celli].x() = rMin_ + (j + 0.5) * dR_;
                BEMcoordinates_[celli].y() = theta;
                BEMcoordinates_[celli].z() = 0.0;
               
                area_[celli] = mathematical::pi* dTheta_ / mathematical::twoPi;
                
                area2_[celli] = dTheta_ / 2 * (sqr( BEMcoordinates_[celli].x() + 0.5 * dR_ )-(sqr( BEMcoordinates_[celli].x() - 0.5 * dR_ )));
                             
            }
        }

        label totCellSize = cells_.size();
        
        reduce( totCellSize, sumOp<label>());
        
        interpolationWeights_.setSize(cells_.size()*nR_*nTheta_,Zero);
        backInterpolationWeights_.setSize(cells_.size()*nR_*nTheta_,Zero);
        sumInterpolationWeights_.setSize(nR_*nTheta_,Zero);
        sumBackInterpolationWeights_.setSize(nR_*nTheta_,Zero);
        
        const scalarField& V = mesh_.V();
        
        for(label i = 0; i < cells_.size(); i++ )
        {
            for( label j = 0; j < nTheta_; j++ )
            {
                for( label k = 0; k < nR_; k++)
                {
                    label BEMcelli = j*nR_ + k;
                    label interpolationCelli = i*nR_*nTheta_ + j*nR_ + k;
                    
                    scalar limitL = dTheta_ * BEMcoordinates_[BEMcelli].x();
                    
                    scalar dT = mag (BEMcoordinates_[BEMcelli].y() - x_[i].y());
                    if ( dT >= mathematical::twoPi ) dT -= mathematical::twoPi;
                    
                    scalar L = dT * BEMcoordinates_[BEMcelli].x();
                    scalar Ltest = ( mathematical::twoPi - dT ) * BEMcoordinates_[BEMcelli].x();
                    
                    scalar dL = 0.0;
                    
                    if(L <= limitL )
                    {
                        dL = sqr(L);
                    }
                    else if( Ltest <= limitL)
                    {
                        dL = sqr(Ltest);
                    }
                    else
                    {
                        continue;
                    }
                    
                    scalar dR = mag(x_[i].x() - BEMcoordinates_[BEMcelli].x());
                    
                    if( dR > dR_) continue;
                    
                    interpolationWeights_[interpolationCelli] = sqr(dL)+sqr(dR);
                    interpolationWeights_[interpolationCelli] = max(interpolationWeights_[interpolationCelli], SMALL);
                    interpolationWeights_[interpolationCelli] = 1.0 / interpolationWeights_[interpolationCelli];
                    sumInterpolationWeights_[BEMcelli] += interpolationWeights_[interpolationCelli];
                    
                    backInterpolationWeights_[interpolationCelli] = V[cells_[i]];// * interpolationWeights_[interpolationCelli];
                    sumBackInterpolationWeights_[BEMcelli] += backInterpolationWeights_[interpolationCelli];
                    
                    
                    
                }
            }
        }
        
        reduce( sumInterpolationWeights_, sumOp<List<scalar>>());
        reduce( sumBackInterpolationWeights_, sumOp<List<scalar>>());
}


Foam::tmp<Foam::vectorField> Foam::fv::rotorDiskUpSource::inflowVelocity
(
    const volVectorField& U
) const
{
    return U.primitiveField();
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::rotorDiskUpSource::rotorDiskUpSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh

)
:
    cellSetOption(name, modelType, dict, mesh),
    rhoRef_(1.0),
    //muRef_(SMALL),
    omega_(0.0),
    omegaMax_(0.0),
    nOmega_(0),
    mode_(""),
    nBlades_(0),
    inletFlow_(ifLocal),
    //inletVelocity_(Zero),
    //tipEffect_(1.0),
    //flap_(),
    x_(cells_.size(), Zero),
    Rcone_(cells_.size(), I),
    area_(),
    coordSys_(),
    rMax_(0.0),
    rMin_(0.0),
    nR_(0),
    nTheta_(0),
    //trim_(trimModelUp::New(*this, coeffs_)),
    blade_(coeffs_.subDict("blade")),
    profiles_(coeffs_.subDict("profiles"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::rotorDiskUpSource::~rotorDiskUpSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::rotorDiskUpSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    // Read the reference density for incompressible flow
    coeffs_.readEntry("rhoRef", rhoRef_);
    changeOmega();

    const vectorField Uin(inflowVelocity(eqn.psi()));
    
    const auto* turbPtr =
            mesh_.findObject<incompressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );
    const volScalarField& nu = turbPtr->nu(); 
    
    const volScalarField& p =  mesh_.lookupObject<volScalarField>("p"); 
    
    if(mesh_.moving())
    {
        //const dictionary& pimple = mesh_.solutionDict().subDict("PIMPLE");
        
        const auto* pimple=mesh_.findObject<pimpleControl>("pimple");
        
        // updating geometry of virtual disk
        
        Info << "Updating geometry of virtual disk" << endl;
        if (pimple->firstIter())
        updateGeometry();
    }
    
    calculate(geometricOneField(), p, Uin, nu, force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().writeTime())
    {
        force.write();
    }
}


void Foam::fv::rotorDiskUpSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );
    
    changeOmega();

    const vectorField Uin(inflowVelocity(eqn.psi()));
    
    const auto* turbPtr =
            mesh_.findObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );
    
    const volScalarField& nu = turbPtr->nu();
    
    const volScalarField& p =  mesh_.lookupObject<volScalarField>("p"); 
    
    if(mesh_.moving())
    {
        //const dictionary& pimple = mesh_.solutionDict().subDict("PIMPLE");
        
        const auto* pimple=mesh_.findObject<pimpleControl>("pimple");
        
        // updating geometry of virtual disk
        
        Info << "Updating geometry of virtual disk" << endl;
        if (pimple->firstIter())
        updateGeometry();
    }
           
    calculate(rho, p, Uin, nu, force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().writeTime())
    {
        force.write();
    }
}


bool Foam::fv::rotorDiskUpSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readEntry("fields", fieldNames_);
        applied_.setSize(fieldNames_.size(), false);

        // Read coordinate system/geometry invariant properties
        
        //muRef_ = coeffs_.get<scalar>("mu");
        mode_ = coeffs_.get<word>("mode");
        
        omega_ = 0.0;
        omegaMax_ = rpmToRads(coeffs_.get<scalar>("rpm"));
        
        
        
        coeffs_.readEntry("nOmega", nOmega_);

        coeffs_.readEntry("nBlades", nBlades_);

        // Create coordinate system
        createCoordinateSystem();

        // Read coordinate system dependent properties
        checkData();
        
        // Set the profile ID for each blade section
        profiles_.connectBlades(blade_.profileName(), blade_.profileID());

        constructGeometry();
        
        rMax_ = coeffs_.get<scalar>("rMax");
        rMin_ = coeffs_.get<scalar>("rMin");
        
        tipFactor_ = coeffs_.get<scalar>("tipFactor");
        gamma_ = coeffs_.get<scalar>("gamma");
        aRef_ = coeffs_.get<scalar>("aRef");
        ReMa_ = coeffs_.get<word>("ReMa");
        
        coeffs_.readEntry("nR", nR_);
        coeffs_.readEntry("nTheta", nTheta_);
        
        BEMcoordinates_.setSize(nR_*nTheta_,Zero);
        area_.setSize(nR_*nTheta_,Zero);
        area2_.setSize(nR_*nTheta_,Zero);
        
        dR_ = ( rMax_ - rMin_ ) / nR_;
        dTheta_ = mathematical::twoPi / nTheta_;
        
        for( label i = 0; i < nTheta_; i++)
        {
            scalar theta = i * dTheta_;
            
            for( label j = 0; j < nR_; j++)
            {
                label celli = i*nR_ + j;
                
                BEMcoordinates_[celli].x() = rMin_ + (j + 0.5) * dR_;
                BEMcoordinates_[celli].y() = theta;
                BEMcoordinates_[celli].z() = 0.0;
               
                area_[celli] = mathematical::pi* dTheta_ / mathematical::twoPi;
                
                area2_[celli] = dTheta_ / 2 * (sqr( BEMcoordinates_[celli].x() + 0.5 * dR_ )-(sqr( BEMcoordinates_[celli].x() - 0.5 * dR_ )));
                             
            }
        }

        label totCellSize = cells_.size();
        
        reduce( totCellSize, sumOp<label>());
        
        interpolationWeights_.setSize(cells_.size()*nR_*nTheta_,Zero);
        backInterpolationWeights_.setSize(cells_.size()*nR_*nTheta_,Zero);
        sumInterpolationWeights_.setSize(nR_*nTheta_,Zero);
        sumBackInterpolationWeights_.setSize(nR_*nTheta_,Zero);
        
        const scalarField& V = mesh_.V();
        
        for(label i = 0; i < cells_.size(); i++ )
        {
            for( label j = 0; j < nTheta_; j++ )
            {
                for( label k = 0; k < nR_; k++)
                {
                    label BEMcelli = j*nR_ + k;
                    label interpolationCelli = i*nR_*nTheta_ + j*nR_ + k;
                    
                    scalar limitL = dTheta_ * BEMcoordinates_[BEMcelli].x();
                    
                    scalar dT = mag (BEMcoordinates_[BEMcelli].y() - x_[i].y());
                    if ( dT >= mathematical::twoPi ) dT -= mathematical::twoPi;
                    
                    scalar L = dT * BEMcoordinates_[BEMcelli].x();
                    scalar Ltest = ( mathematical::twoPi - dT ) * BEMcoordinates_[BEMcelli].x();
                    
                    scalar dL = 0.0;
                    
                    if(L <= limitL )
                    {
                        dL = sqr(L);
                    }
                    else if( Ltest <= limitL)
                    {
                        dL = sqr(Ltest);
                    }
                    else
                    {
                        continue;
                    }
                    
                    scalar dR = mag(x_[i].x() - BEMcoordinates_[BEMcelli].x());
                    
                    if( dR > dR_) continue;
                    
                    interpolationWeights_[interpolationCelli] = sqr(dL)+sqr(dR);
                    interpolationWeights_[interpolationCelli] = max(interpolationWeights_[interpolationCelli], SMALL);
                    interpolationWeights_[interpolationCelli] = 1.0 / interpolationWeights_[interpolationCelli];
                    sumInterpolationWeights_[BEMcelli] += interpolationWeights_[interpolationCelli];
                    
                    backInterpolationWeights_[interpolationCelli] = V[cells_[i]];// * interpolationWeights_[interpolationCelli];
                    sumBackInterpolationWeights_[BEMcelli] += backInterpolationWeights_[interpolationCelli];
                    
                    
                    
                }
            }
        }
        
        reduce( sumInterpolationWeights_, sumOp<List<scalar>>());
        reduce( sumBackInterpolationWeights_, sumOp<List<scalar>>());

        if (debug)
        {
            //writeField("thetag", trim_->thetag()(), true);
            writeField("faceArea", area_, true);
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
