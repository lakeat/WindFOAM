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
    sixDOFMSDbodies

Description
    6-DOF solver for multiple bodies

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "sixDOFMSDbodies.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDOFMSDbodies::setBodies()
{
    // Find if duplicate name existes
    forAll (names_, bodyI)
    {
        for
        (
            label otherBody = bodyI + 1;
            otherBody < names_.size();
            otherBody++
        )
        {
            if (names_[bodyI] == names_[otherBody])
            {
                FatalErrorIn("sixDOFMSDbodies::setBodies()")
                    << "Duplicate names of bodies: this is not allowed"
                    << exit(FatalError);
            }
        }
    }

    odes_.setSize(names_.size());
    solvers_.setSize(names_.size());

    forAll (names_, bodyI)
    {
        odes_.set
        (
            bodyI,
            new sixDOFMSDqODE
            (
                IOobject
                (
                    names_[bodyI],
                    runTime_.timeName(),
                    runTime_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            )
        );

        solvers_.set
        (
            bodyI,
            ODESolver::New
            (
                lookup("solver"),
                odes_[bodyI]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDOFMSDbodies::sixDOFMSDbodies
(
    const Time& runTime
)
:
    IOdictionary
    (
        IOobject
        (
            "sixDOFMSDsolverDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime),
    odes_(),
    solvers_(),
    names_(lookup("bodies"))
{
    setBodies();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sixDOFMSDbodies::solve()
{
    const scalar tol = readScalar(lookup("eps"));

    forAll (odes_, bodyI)
    {
        Info << "Solving 6-DOF MSD for " << names_[bodyI] << " in time"
         << tab << "T = " << runTime_.value() << " s" << endl;

        solvers_[bodyI].solve
        (
            runTime_.value(),
            runTime_.value() + runTime_.deltaT().value(),
            tol,
            runTime_.deltaT().value()
        );
    }
}


const Foam::wordList& Foam::sixDOFMSDbodies::names() const
{
    return names_;
}


const Foam::PtrList<Foam::sixDOFMSDqODE>& Foam::sixDOFMSDbodies::operator()() const
{
    return odes_;
}


Foam::PtrList<Foam::sixDOFMSDqODE>& Foam::sixDOFMSDbodies::getOde()
{
    return odes_;
}


// ************************************************************************* //

