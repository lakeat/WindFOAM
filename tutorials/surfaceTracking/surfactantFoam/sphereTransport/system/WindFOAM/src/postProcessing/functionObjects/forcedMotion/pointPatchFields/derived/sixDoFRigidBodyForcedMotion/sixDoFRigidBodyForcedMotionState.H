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
    Foam::sixDoFRigidBodyForcedMotionState

Description
    Holds the motion state of sixDoF object.  Wrapped up together
    to allow rapid scatter to other processors.  The processors must all
    maintain exactly the same state data to avoid any drift or inconsistency.

SourceFiles
    sixDoFRigidBodyForcedMotionStateI.H
    sixDoFRigidBodyForcedMotionState.C
    sixDoFRigidBodyForcedMotionStateIO.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyForcedMotionState_H
#define sixDoFRigidBodyForcedMotionState_H

#include "vector.H"
#include "point.H"
#include "diagTensor.H"
#include "tensor.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class sixDoFRigidBodyForcedMotionState;
Istream& operator>>(Istream&, sixDoFRigidBodyForcedMotionState&);
Ostream& operator<<(Ostream&, const sixDoFRigidBodyForcedMotionState&);


/*---------------------------------------------------------------------------*\
                  Class sixDoFRigidBodyForcedMotionState Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyForcedMotionState
{
    // Private data

        //- Current position of the centre of mass of the body
        point centreOfMass_;

        //- Orientation, stored as the rotation tensor to transform
        //  from the body to the global reference frame, i.e.:
        //  globalVector = Q_ & bodyLocalVector
        //  bodyLocalVector = Q_.T() & globalVector
        tensor Q_;


public:

    // Constructors

        //- Construct null
        sixDoFRigidBodyForcedMotionState();

        //- Construct from components
        sixDoFRigidBodyForcedMotionState
        (
            const point& centreOfMass,
            const tensor& Q
        );

        //- Construct from dictionary
        sixDoFRigidBodyForcedMotionState(const dictionary& dict);

        //- Construct as copy
        sixDoFRigidBodyForcedMotionState(const sixDoFRigidBodyForcedMotionState&);


    //- Destructor
    ~sixDoFRigidBodyForcedMotionState();


    // Member Functions

        // Access

            //- Return access to the centre of mass
            inline const point& centreOfMass() const;

            //- Return access to the orientation
            inline const tensor& Q() const;


        // Edit

            //- Return non-const access to the centre of mass
            inline point& centreOfMass();

            //- Return non-const access to the orientation
            inline tensor& Q();


        //- Write
        void write(Ostream&) const;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, sixDoFRigidBodyForcedMotionState&);
        friend Ostream& operator<<(Ostream&, const sixDoFRigidBodyForcedMotionState&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFRigidBodyForcedMotionStateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
