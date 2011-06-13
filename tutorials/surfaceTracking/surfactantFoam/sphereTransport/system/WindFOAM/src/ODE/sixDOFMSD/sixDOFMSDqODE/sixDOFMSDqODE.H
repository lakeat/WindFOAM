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
    sixDOFMSDqODE

Description
    6-DOF solver using quaternions

Author
    Daniel Wei

SourceFiles
    sixDOFMSDqODEI.H
    sixDOFMSDqODE.C

\*---------------------------------------------------------------------------*/

#ifndef sixDOFMSDqODE_H
#define sixDOFMSDqODE_H

#include "ODE.H"
#include "IOdictionary.H"
#include "dimensionedTypes.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sixDOFMSDqODE Declaration
\*---------------------------------------------------------------------------*/

class sixDOFMSDqODE
:
    public IOdictionary,
    public ODE
{
    // Private data

        // Body data

            //- Mass
            dimensionedScalar mass_;

            //- Rotational moment of inertia around centre of mass
            //  in body (relative) coordinates - aligned with main axes
            dimensionedDiagTensor momentOfInertia_;


        // Platform variables

            //- Translation spring equilibrium position
            const dimensionedVector Xequilibrium_;

            //- Rotation spring equilibrium position
            const dimensionedVector Tequilibrium_;

            //- Translation spring coeffs
            const dimensionedDiagTensor linSpringCoeffs_;

            //- Translation damping coeffs
            const dimensionedDiagTensor linDampingCoeffs_;

            //- Rotation spring coeffs
            const dimensionedDiagTensor rotSpringCoeffs_;

            //- Rotation damping coeffs
            const dimensionedDiagTensor rotDampingCoeffs_;


        // Initial body state variables

            //- Translational displacement relative to spring equilibrium
            dimensionedVector Xrel_;

            //- Rotation displacement relative to spring equilibrium
            dimensionedVector Trel_;

            //- Translational displacement relative to spring equilibrium
            dimensionedVector Xdiff_;

            //- Rotation displacement relative to spring equilibrium
            dimensionedVector Tdiff_;

            //- Translational velocity of mass centroid
            dimensionedVector U_;

            //- Rotational velocity of mass centroid
            dimensionedVector omega_;


        // External forces

            //- Force driving the motion in inertial coord. sys.
            dimensionedVector force_;

            //- Moment driving the motion in inertial coord. sys.
            dimensionedVector moment_;


        //- ODE coefficients
        scalarField coeffs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDOFMSDqODE(const sixDOFMSDqODE&);

        //- Disallow default bitwise assignment
        void operator=(const sixDOFMSDqODE&);


        //- Set ODE coefficients from position and rotation
        void setCoeffs();


        // Variables in relative coordinate system (solved for)

            //- Return acceleration in relative coordinate system
            //  given current values of relative translational displacement and
            //  translational velocity
            dimensionedVector A
            (
                const dimensionedVector& xR,
                const dimensionedVector& uR
            ) const;


            //- Return rotational acceleration in relative coordinate system
            //  given current values of relative rotational displacement and
            //  rotational velocity
            dimensionedVector omegaDot
            (
                const dimensionedVector& thetaR,
                const dimensionedVector& omegaR
            ) const;


public:


//     //- Runtime type information
//       Not possible because of I/O error: incorrect type, expecting dictionary
//       HJ, 11/Feb/2008
//     TypeName("sixDOFMSDqODE");


    // Constructors

        //- Construct from dictionary
        sixDOFMSDqODE(const IOobject& io);


    // Destructor

        virtual ~sixDOFMSDqODE();


    // Member Functions

        //- Return mass
        inline const dimensionedScalar& mass() const;

        //- Return moment of inertia
        inline const dimensionedDiagTensor& momentOfInertia() const;

        //- Return equilibrium position of origin
        inline const dimensionedVector& Xequilibrium() const;

        //- Return equilibrium position of origin
        inline const dimensionedVector& Tequilibrium() const;


        // Variables in relative coordinate system (solved for)

            //- Return displacement in translated coordinate system
            //  relative to spring equilibrium
            inline const dimensionedVector& Xrel() const;

            //- Return displacement in translated coordinate system
            //  relative to spring equilibrium
            inline const dimensionedVector& Trel() const;

            //- Return rotational velocity in relative coordinate system
            inline const dimensionedVector& U() const;

            //- Return rotational acceleration in relative coordinate system
            inline dimensionedVector A() const;

            //- Return rotational velocity in relative coordinate system
            inline const dimensionedVector& omega() const;

            //- Return rotational acceleration in relative coordinate system
            inline dimensionedVector omegaDot() const;


        // Displacement and rotation in the absolute coordinate system

            //- Return position of origin in absolute coordinate system
            inline dimensionedVector X() const;

            //- Return position of origin in absolute coordinate system
            inline dimensionedVector T() const;


        // Displacement and rotation increase in deltaT

            //- Return position of origin in absolute coordinate system
            inline dimensionedVector Xdiff() const;

            //- Return position of origin in absolute coordinate system
            inline dimensionedVector Tdiff() const;


        // Non-constant access to velocities

            //- Return access to velocity of origin
            inline dimensionedVector& U();

            //- Return access to rotational velocity in relative
            //  coordinate system
            inline dimensionedVector& omega();


        // Force and moment

            //- Return force in inertial coord. sys.
            inline const dimensionedVector& force() const;

            //- Return access to force in inertial coord. sys.
            inline dimensionedVector& force();

            //- Return access to force in inertial coord. sys.
            void editForce(dimensionedVector fTmp)
            {
                force_ = fTmp;
            }

            //- Return moment in inertial coord. sys.
            inline const dimensionedVector& moment() const;

            //- Return access to moment in inertial coord. sys.
            inline dimensionedVector& moment();

            //- Return access to moment in inertial coord. sys.
            void editMoment(dimensionedVector mTmp)
            {
                moment_ = mTmp;
            }


        // ODE parameters

            //- Return number of equations
            virtual label nEqns() const
            {
                return 12;
            }

            //- Return access to coefficients
            virtual scalarField& coeffs()
            {
                return coeffs_;
            }

            //- Return reference to coefficients
            virtual const scalarField& coeffs() const
            {
                return coeffs_;
            }

            //- Return derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Return Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta);


        //- WriteData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const sixDOFMSDqODE&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDOFMSDqODEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

