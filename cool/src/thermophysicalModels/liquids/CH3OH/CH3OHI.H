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

\*---------------------------------------------------------------------------*/

inline Foam::scalar Foam::CH3OH::rho(scalar p, scalar T) const
{
    return rho_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::pv(scalar p, scalar T) const
{
    return pv_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::hl(scalar p, scalar T) const
{
    return hl_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::cp(scalar p, scalar T) const
{
    return cp_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::h(scalar p, scalar T) const
{
    return h_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::cpg(scalar p, scalar T) const
{
    return cpg_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::B(scalar p, scalar T) const
{
    return B_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::mu(scalar p, scalar T) const
{
    return mu_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::mug(scalar p, scalar T) const
{
    return mug_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::K(scalar p, scalar T) const
{
    return K_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::Kg(scalar p, scalar T) const
{
    return Kg_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::sigma(scalar p, scalar T) const
{
    return sigma_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::D(scalar p, scalar T) const
{
    return D_.f(p, T);
}


inline Foam::scalar Foam::CH3OH::D(scalar p, scalar T, scalar Wb) const
{
    return D_.f(p, T, Wb);
}


// ************************************************************************* //
