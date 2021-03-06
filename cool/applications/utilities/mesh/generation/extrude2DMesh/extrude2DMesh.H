/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::extrude2DMesh

Description
    Given 2D mesh insert all the topology changes to extrude. Does not work
    in parallel

SourceFiles
    extrude2DMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extrude2DMesh_H
#define extrude2DMesh_H

#include "label.H"
#include "scalar.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class directTopoChange;
class mapPolyMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class extrude2DMesh Declaration
\*---------------------------------------------------------------------------*/

class extrude2DMesh
{
    // Private data

        //- Reference to 2D mesh
        const polyMesh& mesh_;


public:

    //- Runtime type information
    ClassName("extrude2DMesh");


    // Constructors

        //- Construct from mesh
        extrude2DMesh(const polyMesh&);


    // Member Functions

        //- Play commands into directTopoChange to extrude mesh.
        void setRefinement
        (
            const direction extrudeDir,
            const scalar thickness,
            const label frontPatchI,
            directTopoChange&
        ) const;

        //- Force recalculation of locally stored data on topological change
        void updateMesh(const mapPolyMesh&)
        {}

        //- Force recalculation of locally stored data for mesh distribution
        void distribute(const mapDistributePolyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
