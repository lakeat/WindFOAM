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
    Foam::lduScheduleEntry

Description
    Struct to hold the patch index and the initialisation flag for
    the patch schedule

\*---------------------------------------------------------------------------*/

#ifndef lduSchedule_H
#define lduSchedule_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class lduScheduleEntry;
Ostream& operator<<(Ostream& os, const lduScheduleEntry& lb);


/*---------------------------------------------------------------------------*\
                      Class lduSchedule Declaration
\*---------------------------------------------------------------------------*/

struct lduScheduleEntry
{
    label patch;
    bool init;

    friend bool operator!=
    (
        const lduScheduleEntry& lb1,
        const lduScheduleEntry& lb2
    )
    {
        return true;
    }

    friend Ostream& operator<<(Ostream& os, const lduScheduleEntry& lb)
    {
        os  << lb.patch << token::SPACE << lb.init;
        return os;
    }
};


typedef List<lduScheduleEntry> lduSchedule;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
