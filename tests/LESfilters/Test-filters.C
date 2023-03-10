/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Test-volField

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "LESfilter.H"

#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "basicKinematicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("zip", "run zip/unzip tests");

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"


    IOdictionary filterDict
    (
        IOobject
        (
            "filterDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "dict" << filterDict << nl;

    autoPtr<LESfilter>  filterPtr(LESfilter::New(mesh, filterDict));

    const LESfilter& filter(filterPtr());

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "Evolving " << kinematicCloud.name() << endl;

        kinematicCloud.evolve();

        kinematicCloud.write();
        
        volScalarField alphap("theta", kinematicCloud.theta());

        Info<< "theta! " << endl;

        alphap.write();

        volScalarField  fAlphap("fAlphap", filter(alphap));

        fAlphap.write(); 
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
