/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    reconstructPar

Description
    Reconstructs fields of a case that is decomposed for parallel
    execution of OpenFOAM.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "IOobjectList.H"
#include "processorMeshes.H"

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "coordSetWriter.H"

#include "couplingFilter.H"

#include "basicKinematicCloud.H"
#define basicKinematicTypeCloud basicKinematicCloud


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



int main(int argc, char *argv[])
{

    // Enable -constant ... if someone really wants it
    // Enable -withZero to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "filter.H"
    #include "readGravitationalAcceleration.H"

    volScalarField rhoc
    (
        IOobject
        (
            "rhoc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimDensity, 0.0)
    );

    volScalarField muc
    (
        IOobject
        (
            "muc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimViscosity, 0.0)
    );


    volScalarField alphap
    (
        IOobject
        (
            "alphap",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );


    volVectorField Up
    (
        IOobject
        (
            "Up",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

    volScalarField alphapMean
    (
        IOobject
        (
            "alphap.Mean",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );


    volVectorField UpMean
    (
        IOobject
        (
            "Up.mean",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

    volVectorField UpP2M
    (
        IOobject
        (
            "Up.p2m",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
    );

    volVectorField UlMean
    (
        IOobject
        (
            "Ul.mean",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

    volVectorField UlP2M
    (
        IOobject
        (
            "Ul.p2m",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity*dimVelocity, vector::zero)
    );


    // Loop over all times
    forAll(timeDirs, timeI)
    {
        // Set time for global database
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl << endl;

        mesh.readUpdate();

        basicKinematicCloud kinematicCloud
        (
            "kinematicCloud",
            rhoc,
            Up,
            muc,
            g
        );

        Info<< "    Read " << returnReduce(kinematicCloud.size(), sumOp<label>())
            << " particles" << endl;

        alphap ==  dimensionedScalar("zero", dimless, 0.0);
        Up == dimensionedVector("zero", dimVelocity, vector::zero);

        for (const auto& p : kinematicCloud)
        {
            alphap[p.cell()] += p.nParticle()*p.volume();
            Up[p.cell()] += p.nParticle()*p.volume()*p.U();
        }

        filterModel.filter(alphap);
        filterModel.filter(Up);

        alphapMean += alphap;

        volScalarField alphaLimited(max(alphap,SMALL));
        Up /= alphaLimited;
        UpMean += Up;
        UpP2M.component(vector::X) = UpP2M.component(vector::X) + sqr(Up.component(vector::X));
        UpP2M.component(vector::Y) = UpP2M.component(vector::Y) + sqr(Up.component(vector::Y));
        UpP2M.component(vector::Z) = UpP2M.component(vector::Z) + sqr(Up.component(vector::Z));

        IOobject UlHeader
        (
            "U.water",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        volVectorField Ul(UlHeader,mesh);   //Ul is already the phase-averaged quantity. 
        UlMean += Ul;

        UlP2M.component(vector::X) = UlP2M.component(vector::X) + sqr(Ul.component(vector::X));
        UlP2M.component(vector::Y) = UlP2M.component(vector::Y) + sqr(Ul.component(vector::Y));
        UlP2M.component(vector::Z) = UlP2M.component(vector::Z) + sqr(Ul.component(vector::Z));

    }

	Info << "Dividing the results" << endl ;

    UpMean  /= timeDirs.size();;
    UpP2M  /= timeDirs.size();;
    alphapMean /= timeDirs.size();

    UpP2M.component(vector::X) = UpP2M.component(vector::X) - sqr(UpMean.component(vector::X));
    UpP2M.component(vector::Y) = UpP2M.component(vector::Y) -  sqr(UpMean.component(vector::Y));
    UpP2M.component(vector::Z) = UpP2M.component(vector::Z) -  sqr(UpMean.component(vector::Z));

    UlMean  /= timeDirs.size();
    UlP2M.component(vector::X) = UlP2M.component(vector::X) -  sqr(UlMean.component(vector::X));
    UlP2M.component(vector::Y) = UlP2M.component(vector::Y) - sqr(UlMean.component(vector::Y));
    UlP2M.component(vector::Z) = UlP2M.component(vector::Z)  - sqr(UlMean.component(vector::Z));

    alphapMean.write();
    UpMean.write();
    UpP2M.write();
    UlMean.write();
    UlP2M.write();


    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
