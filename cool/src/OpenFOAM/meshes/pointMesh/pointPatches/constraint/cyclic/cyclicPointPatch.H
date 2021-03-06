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
    Foam::cyclicPointPatch

Description
    Cyclic patch for post-processing.

SourceFiles
    cyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPointPatch_H
#define cyclicPointPatch_H

#include "coupledFacePointPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicPolyPatch& cyclicPolyPatch_;


    // Demand driven private data

        //- The set of pairs of points that require transformation
        //  and/or mapping
        edgeList transformPairs_;


    // Private Member Functions

        //- Disallow default construct as copy
        cyclicPointPatch(const cyclicPointPatch&);

        //- Disallow default assignment
        void operator=(const cyclicPointPatch&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry()
        {}

        //- Correct patches after moving points
        virtual void movePoints()
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


public:

    typedef coupledPointPatch CoupledPointPatch;
    typedef pointBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName(cyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~cyclicPointPatch();


    // Member Functions

        // Access

            //- Return the underlying cyclicPolyPatch
            const cyclicPolyPatch& cyclicPatch() const
            {
                return cyclicPolyPatch_;
            }

            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return cyclicPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                return cyclicPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                return cyclicPolyPatch_.reverseT();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
