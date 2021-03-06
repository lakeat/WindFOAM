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
    Foam::pyrMatcher

Description
    A cellMatcher for pyr cells

See Also
    cellMatcher

SourceFiles
    pyrMatcher.C

\*---------------------------------------------------------------------------*/

#ifndef pyrMatcher_H
#define pyrMatcher_H

#include "cellMatcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pyrMatcher Declaration
\*---------------------------------------------------------------------------*/

class pyrMatcher
:
    public cellMatcher
{
    // Static data members

        //- constants for this shape
        static const label vertPerCell;
        static const label facePerCell;
        static const label maxVertPerFace;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pyrMatcher(const pyrMatcher&);

        //- Disallow default bitwise assignment
        void operator=(const pyrMatcher&);


public:

    // Constructors

        //- Construct null
        pyrMatcher();

    // Destructor

        ~pyrMatcher();


    // Member Functions

        virtual label nVertPerCell() const
        {
            return vertPerCell;
        }

        virtual label nFacePerCell() const
        {
            return facePerCell;
        }

        virtual label nMaxVertPerFace() const
        {
            return maxVertPerFace;
        }

        virtual label faceHashValue() const;

        virtual bool faceSizeMatch(const faceList&, const labelList&) const;

        virtual bool matchShape
        (
            const bool checkOnly,
            const faceList& faces,
            const labelList& faceOwner,
            const label cellI,
            const labelList& myFaces
        );

        virtual bool isA(const primitiveMesh& mesh, const label cellI);

        virtual bool isA(const faceList&);

        virtual bool matches
        (
            const primitiveMesh& mesh,
            const label cellI,
            cellShape& shape
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
