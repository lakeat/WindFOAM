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
    Foam::sixDoFRigidBodyForcedMotion

Description
    Six degree of freedom motion for a rigid body.  Angular momentum
    stored in body fixed reference frame.  Reference orientation of
    the body (where Q = I) must align with the cartesian axes such
    that the Inertia tensor is in principle component form.

    Symplectic motion as per:

    title = {Symplectic splitting methods for rigid body molecular dynamics},
    publisher = {AIP},
    year = {1997},
    journal = {The Journal of Chemical Physics},
    volume = {107},
    number = {15},
    pages = {5840-5851},
    url = {http://link.aip.org/link/?JCP/107/5840/1},
    doi = {10.1063/1.474310}

    Can add restraints (i.e. a spring) and constraints (i.e. motion
    may only be on a plane).

SourceFiles
    sixDoFRigidBodyForcedMotionI.H
    sixDoFRigidBodyForcedMotion.C
    sixDoFRigidBodyForcedMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyForcedMotion_H
#define sixDoFRigidBodyForcedMotion_H

#include "sixDoFRigidBodyForcedMotionState.H"
#include "pointField.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class sixDoFRigidBodyForcedMotion;
Istream& operator>>(Istream&, sixDoFRigidBodyForcedMotion&);
Ostream& operator<<(Ostream&, const sixDoFRigidBodyForcedMotion&);


/*---------------------------------------------------------------------------*\
                      Class sixDoFRigidBodyForcedMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyForcedMotion
{
    // Private data

        //- Motion state data object
        sixDoFRigidBodyForcedMotionState motionState_;

        //- Centre of mass of initial state
        point initialCentreOfMass_;

        //- Orientation of initial state
        tensor initialQ_;

        //- Translation Amplitude
        vector translationAmplitude_;

        //- Translation Frequency
        vector translationFrequency_;

        //- Rotation Amplitude, unit in radian
        scalar rotationAmplitude_;

        //- Rotation Frequency, unit in radian
        scalar rotationFrequency_;


    // Private Member Functions

        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return access to the motion state
            inline const sixDoFRigidBodyForcedMotionState& motionState() const;

            //- Return access to the initial centre of mass
            inline const point& initialCentreOfMass() const;

            //- Return access to the initial orientation
            inline const tensor& initialQ() const;

            //- Return access to the orientation
            inline const tensor& Q() const;

            //- Return access to the Translation Amplitude
            inline const vector& translationAmplitude() const;

            //- Return access to the Translation Frequency
            inline const vector& translationFrequency() const;

            //- Return access to the Rotation Amplitude
            inline const scalar& rotationAmplitude() const;

            //- Return access to the Rotation Frequency
            inline const scalar& rotationFrequency() const;


        // Edit

            //- Return non-const access to the centre of mass
            inline point& initialCentreOfMass();

            //- Return non-const access to the centre of mass
            inline tensor& initialQ();

            //- Return non-const access to the orientation
            inline tensor& Q();

            //- Return non-const access to the Translation Amplitude
            inline vector& translationAmplitude();

            //- Return non-const access to the Translation Frequency
            inline vector& translationFrequency();

            //- Return non-const access to the Rotation Amplitude
            inline scalar& rotationAmplitude();

            //- Return non-const access to the Rotation Frequency
            inline scalar& rotationFrequency();


public:

    // Constructors

        //- Construct null
        sixDoFRigidBodyForcedMotion();

        //- Construct from components
        sixDoFRigidBodyForcedMotion
        (
            const point& centreOfMass,
            const tensor& Q,
            const point& initialCentreOfMass,
            const tensor& initialQ,
            const vector translationAmplitude,
            const vector translationFrequency,
            const scalar rotationAmplitude,
            const scalar rotationFrequency
        );

        //- Construct from dictionary
        sixDoFRigidBodyForcedMotion
        (
            const dictionary& dict
        );

        //- Construct as copy
        sixDoFRigidBodyForcedMotion(const sixDoFRigidBodyForcedMotion&);


    //- Destructor
    ~sixDoFRigidBodyForcedMotion();


    // Member Functions

        //- First leapfrog velocity adjust and motion part, required
        // before force calculation
        void updatePosition
        (
            scalar t,
            scalar deltaT
        );

        //- Transform the given initial state pointField by the current
        //  motion state
        inline tmp<pointField> currentPosition
        (
            const pointField& pInitial
        ) const;

        //- Transform the given initial state point by the current motion
        //  state
        inline point currentPosition(const point& pInitial) const;

        //- Access the orientation tensor, Q.
        //  globalVector = Q & bodyLocalVector
        //  bodyLocalVector = Q.T() & globalVector
        inline const tensor& orientation() const;
        
        inline tensor& orientation();

        //- Report the status of the motion
        void status() const;


        // Access

            //- Return const access to the centre of mass
            inline const point& centreOfMass() const;


        // Edit

            //- Return non-const access to the centre of mass
            inline point& centreOfMass();


        //- Write
        void write(Ostream&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, sixDoFRigidBodyForcedMotion&);
        friend Ostream& operator<<(Ostream&, const sixDoFRigidBodyForcedMotion&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFRigidBodyForcedMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
