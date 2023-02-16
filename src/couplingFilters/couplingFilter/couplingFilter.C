/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "couplingFilter.H"
#include "fvcLaplacianPlus.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::couplingFilter::couplingFilter
(
    const fvMesh& mesh,
    Time& diffusionRunTime,
    const fvMesh& diffusionmesh,
    simpleControl& simple
)
:   IOdictionary
    (
        IOobject
        (
            "couplingDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),    
    mesh_(mesh),
    diffusionRunTime_(diffusionRunTime),
    diffusionMesh_(diffusionmesh),
    simple_(simple),
    implicitFvm_(*this, "useImplicitLaplacian", true),
    smoothDirection_
    (
        this->lookupOrDefault
        (
            "smoothDirection",
            tensor(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0)
        )
    ),
    DT("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_),
    startTime(diffusionRunTime_.startTime()),
    startTimeIndex(diffusionRunTime_.startTimeIndex()),
    diffusionBandWidth_(this->lookupOrDefault<scalar>("diffusionBandWidth", 0.003)),
    diffusionSteps_(this->lookupOrDefault("diffusionSteps", 5)),
    diffusionTime_(0),
    diffusionDeltaT_(0),
    adjustDeltaT_(*this, "adjustDiffusionSteps", true),
    tempDiffScalar_
    (
        IOobject
        (
            "tempDiffScalar",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedScalar
        (
            "zero",
            dimless,
            scalar(0.0)
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    tempDiffVector_
    (
        IOobject
        (
            "tempDiffVector",
            diffusionRunTime_.timeName(),
            diffusionMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        diffusionMesh_,
        dimensionedVector
        (
            "zero",
            dimVelocity,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
        
    ),
    tempDiffScalarInterFeildRef_(tempDiffScalar_.internalFieldRef()),
    tempDiffVectorInterFeildRef_(tempDiffVector_.internalFieldRef())
{
    //- initialization information output
    if (!this->found("diffusionBandWidth"))
    {
        Info<< "Diffusion band width for phase set to 0.003 by default" << endl;
    }
    else
    {
        Info<< "Diffusion band width set to " << diffusionBandWidth_ << endl;
    }
    
    // determine the time and time step in diffusion procedure
    diffusionTime_ = pow(diffusionBandWidth_, 2)/4;
    diffusionDeltaT_ = diffusionTime_/diffusionSteps_;

}
    
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// const access
inline Foam::scalar
Foam::couplingFilter::diffusionBandWidth() const
{
    return diffusionBandWidth_;
}


inline Foam::label
Foam::couplingFilter::diffusionSteps() const
{
    return diffusionSteps_;
}


inline Foam::scalar
Foam::couplingFilter::diffusionDeltaT() const
{
    return diffusionDeltaT_;
}


inline Foam::scalar
Foam::couplingFilter::diffusionTime() const
{
    return diffusionTime_;
}


// access to diffusion settings
inline Foam::scalar&
Foam::couplingFilter::diffusionBandWidth()
{
    return diffusionBandWidth_;
}


inline Foam::label&
Foam::couplingFilter::diffusionSteps()
{
    return diffusionSteps_;
}


inline Foam::scalar&
Foam::couplingFilter::diffusionDeltaT()
{
    return diffusionDeltaT_;
}


inline Foam::scalar&
Foam::couplingFilter::diffusionTime()
{
    return diffusionTime_;
}


void Foam::couplingFilter::filter
(
    volScalarField::Internal& s
)
{    


    if(tempDiffScalar_.dimensions() != s.dimensions())
        tempDiffScalar_.dimensions().reset(s.dimensions());

    diffusionRunTime_.setTime(startTime, startTimeIndex);
    diffusionRunTime_.setEndTime(diffusionTime_);
    diffusionRunTime_.setDeltaT(diffusionDeltaT_);
    
    tempDiffScalarInterFeildRef_ = s;
    tempDiffScalar_.correctBoundaryConditions();
    
    if (implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(tempDiffScalar_) - fvm::laplacian(DT, tempDiffScalar_,"couplingFilterDiffusion"));
                    tempDiffScalar_.correctBoundaryConditions();
                }
            }
            else
            {
                solve(fvm::ddt(tempDiffScalar_) - fvm::laplacian(DT, tempDiffScalar_,"couplingFilterDiffusion"));
                tempDiffScalar_.correctBoundaryConditions();
            }
        }

    }
    else
    {

        label stepIndex = 1;// 1-10 diffusionDeltaT_/10; 11-20, diffusionDeltaT_/5, 
                            //>21 diffusionDeltaT_/2 >51 diffusionDeltaT_

        while (diffusionRunTime_.loop())
        {

            if (adjustDeltaT_.value())
            {
                scalar adjustDiffusionDeltaT = diffusionDeltaT_;
            
                if (stepIndex <= 4 )
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/20;
                }
                else if (stepIndex <= 10)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/7.5;
                }
                else if (stepIndex <= 15)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/5;
                }
                else if (stepIndex <= 21)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/2;
                }
                else
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_;
                }
                
                diffusionRunTime_.setDeltaT(adjustDiffusionDeltaT);
            }

            tempDiffScalar_ += fvc::laplacian(DT, tempDiffScalar_,"couplingFilterDiffusion");
            tempDiffScalar_.correctBoundaryConditions();

            stepIndex++;
        
        }
    }
    
    scalarField& sRef = s;
    sRef = tempDiffScalarInterFeildRef_;
    
}

void Foam::couplingFilter::filter
(
    volVectorField::Internal& s
)
{

    if(tempDiffVector_.dimensions() != s.dimensions())
        tempDiffVector_.dimensions().reset(s.dimensions());

    diffusionRunTime_.setTime(startTime, startTimeIndex);
    diffusionRunTime_.setEndTime(diffusionTime_);
    diffusionRunTime_.setDeltaT(diffusionDeltaT_);
    
    tempDiffVectorInterFeildRef_ = s;
    tempDiffVector_.correctBoundaryConditions();

    if (implicitFvm_.value())
    {
        while (diffusionRunTime_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(tempDiffVector_) - fvm::laplacian(DT, tempDiffVector_,"couplingFilterDiffusion"));
                    tempDiffVector_.correctBoundaryConditions();
                }
            }
            else
            {
                solve(fvm::ddt(tempDiffVector_) - fvm::laplacian(DT, tempDiffVector_,"couplingFilterDiffusion"));
                tempDiffVector_.correctBoundaryConditions();
            }

        }

    }
    else
    {
        label stepIndex = 1;// 1-10 diffusionDeltaT_/10; 11-20, diffusionDeltaT_/5, 
                            //>21 diffusionDeltaT_/2 >51 diffusionDeltaT_

        while (diffusionRunTime_.loop())
        {

            if (adjustDeltaT_.value())
            {
                scalar adjustDiffusionDeltaT = diffusionDeltaT_;
            
                if (stepIndex <= 4 )
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/20;
                }
                else if (stepIndex <= 10)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/7.5;
                }
                else if (stepIndex <= 15)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/5;
                }
                else if (stepIndex <= 21)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/2;
                }
                else
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_;
                }
                
                diffusionRunTime_.setDeltaT(adjustDiffusionDeltaT);
            }

            tempDiffVector_ += fvc::laplacian(DT, tempDiffVector_,"couplingFilterDiffusion");
            tempDiffVector_.correctBoundaryConditions();

            stepIndex++;
        
        }
    }
    
    vectorField& sRef = s;
    sRef = tempDiffVectorInterFeildRef_;

}



// Return the filtered field from the given volScalarField F

void Foam::couplingFilter::filter
(
    volScalarField& F
)
{
    filter(F.internalFieldRef()); 
    F.correctBoundaryConditions();
}


void Foam::couplingFilter::filter
(
    volVectorField& F
)
{
    filter(F.internalFieldRef()); 
    F.correctBoundaryConditions();
}

// Return the filtered field from the given volScalarField F
tmp<volScalarField>
Foam::couplingFilter::filteredField(const volScalarField& F)
{

    tmp<volScalarField> tF
    (
        volScalarField::New
        (
            "tF",
            F
        )
    );

    volScalarField& S = tF.ref();

    filter(S);

    S.correctBoundaryConditions();

    return tF;

}

// Return the filtered field from the given volScalarField F
tmp<volVectorField>
Foam::couplingFilter::filteredField(const volVectorField& F)
{

    tmp<volVectorField> tF
    (
        volVectorField::New
        (
            "tF",
            F
        )
    );

    volVectorField& S = tF.ref();

    filter(S);

    S.correctBoundaryConditions();

    return tF;
    
}


tmp<volScalarField>
Foam::couplingFilter::filteredField
(
    const volScalarField& alpha,
    const volScalarField& F
)
{

    if(tempDiffScalar_.dimensions() != F.dimensions())
        tempDiffScalar_.dimensions().reset(F.dimensions());

    diffusionRunTime_.setTime(startTime, startTimeIndex);   
    diffusionRunTime_.setEndTime(diffusionTime_);
    diffusionRunTime_.setDeltaT(diffusionDeltaT_);
    
    tmp<volScalarField> tF
    (
        volScalarField::New
        (
            "tF",
            F
        )
    );

    volScalarField& S = tF.ref();
    
    scalarField& iS = S;
    
    tempDiffScalarInterFeildRef_ = iS;
    
    if (implicitFvm_.value())
    {
        while (simple_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(tempDiffScalar_) - fvm::laplacian(fvc::interpolate(alpha)*DT, tempDiffScalar_,"couplingFilterDiffusion"));
                    tempDiffScalar_.correctBoundaryConditions();
                }
            }
            else
            {
                solve(fvm::ddt(tempDiffScalar_) - fvm::laplacian(fvc::interpolate(alpha)*DT, tempDiffScalar_,"couplingFilterDiffusion"));
                tempDiffScalar_.correctBoundaryConditions();
            }
        }
    }
    else
    {
        label stepIndex = 1;// 1-10 diffusionDeltaT_/10; 11-20, diffusionDeltaT_/5, 
                            //>21 diffusionDeltaT_/2 >51 diffusionDeltaT_

        while (diffusionRunTime_.loop())
        {

            if (adjustDeltaT_.value())
            {
                scalar adjustDiffusionDeltaT = diffusionDeltaT_;
            
                if (stepIndex <= 4 )
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/20;
                }
                else if (stepIndex <= 10)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/7.5;
                }
                else if (stepIndex <= 15)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/5;
                }
                else if (stepIndex <= 21)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/2;
                }
                else
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_;
                }
                
                diffusionRunTime_.setDeltaT(adjustDiffusionDeltaT);
            }

            tempDiffScalar_ += fvc::laplacian(fvc::interpolate(alpha)*DT, tempDiffScalar_,"couplingFilterDiffusion");
            tempDiffScalar_.correctBoundaryConditions();

            stepIndex++;
        
        }

    }

    iS = tempDiffScalarInterFeildRef_;
    
    S.correctBoundaryConditions();
    
    return tF;


}


tmp<volVectorField>
Foam::couplingFilter::filteredField
(
    const volScalarField& alpha,
    const volVectorField& F
)
{

    if(tempDiffVector_.dimensions() != F.dimensions())
        tempDiffVector_.dimensions().reset(F.dimensions());

    diffusionRunTime_.setTime(startTime, startTimeIndex);     
    diffusionRunTime_.setEndTime(diffusionTime_);
    diffusionRunTime_.setDeltaT(diffusionDeltaT_);
    
    tmp<volVectorField> tF
    (
        volVectorField::New
        (
            "tF",
            F
        )
    );

    volVectorField& S = tF.ref();
    
    vectorField& iS = S;
    
    tempDiffVectorInterFeildRef_ = iS;
    
    if (implicitFvm_.value())
    {
        while (simple_.loop())
        {
            if (diffusionRunTime_.timeIndex() == 1)
            {
                while (simple_.correctNonOrthogonal())
                {
                    solve(fvm::ddt(tempDiffVector_) - fvm::laplacian(fvc::interpolate(alpha)*DT, tempDiffVector_,"couplingFilterDiffusion"));
                    tempDiffVector_.correctBoundaryConditions();
                }
            }
            else
            {
                solve(fvm::ddt(tempDiffVector_) - fvm::laplacian(fvc::interpolate(alpha)*DT, tempDiffVector_,"couplingFilterDiffusion"));
                tempDiffVector_.correctBoundaryConditions();
            }
        }
    }
    else
    {

        label stepIndex = 1;// 1-10 diffusionDeltaT_/10; 11-20, diffusionDeltaT_/5, 
                            //>21 diffusionDeltaT_/2 >51 diffusionDeltaT_

        while (diffusionRunTime_.loop())
        {

            if (adjustDeltaT_.value())
            {
                scalar adjustDiffusionDeltaT = diffusionDeltaT_;
            
                if (stepIndex <= 4 )
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/20;
                }
                else if (stepIndex <= 10)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/7.5;
                }
                else if (stepIndex <= 15)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/5;
                }
                else if (stepIndex <= 21)
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_/2;
                }
                else
                {
                    adjustDiffusionDeltaT = diffusionDeltaT_;
                }
                
                diffusionRunTime_.setDeltaT(adjustDiffusionDeltaT);
            }

            tempDiffVector_ += fvc::laplacian(fvc::interpolate(alpha)*DT, tempDiffVector_,"couplingFilterDiffusion");
            tempDiffVector_.correctBoundaryConditions();

            stepIndex++;
        }
        
    }

    iS = tempDiffVectorInterFeildRef_;
    
    S.correctBoundaryConditions();
    
    return tF;

}

// ************************************************************************* //
