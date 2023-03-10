/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    
InNamespace
    Foam::fvc

Description
    Calculate the laplacian gradient of the given field.

SourceFiles
    fvcLaplacianPlus.C

\*---------------------------------------------------------------------------*/


#ifndef fvcLaplacianPlus_H
#define fvcLaplacianPlus_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    
    template<class Type, class GType>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    laplacianScheme2Grad
    (
        const GeometricField<GType, fvsPatchField, volMesh>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );


    template<class Type, class GType>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    laplacianScheme2Grad
    (
        const GeometricField<GType, fvsPatchField, surfaceMesh>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );


    template<class Type, class GType>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    laplacianScheme2Grad
    (
        const dimensioned<GType>& gamma,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );


    template<class Type, class GType>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    laplacianScheme2Grad
    (
        const tmp<GeometricField<GType, fvsPatchField, surfaceMesh>>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcLaplacianPlus.C"
#endif

// ************************************************************************* //

#endif

// ************************************************************************* //
