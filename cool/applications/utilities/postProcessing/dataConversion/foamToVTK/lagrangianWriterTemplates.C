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

#include "lagrangianWriter.H"
#include "writeFuns.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::lagrangianWriter::writeIOField(const wordList& objects)
{
    forAll(objects, i)
    {
        const word& object = objects[i];

        IOobject header
        (
            object,
            vMesh_.mesh().time().timeName(),
            cloud::prefix/cloudName_,
            vMesh_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        IOField<Type> fld(header);

        os_ << object << ' ' << pTraits<Type>::nComponents << ' '
            << fld.size() << " float" << std::endl;

        DynamicList<floatScalar> fField(pTraits<Type>::nComponents*fld.size());

        writeFuns::insert(fld, fField);

        writeFuns::write(os_, binary_, fField);
    }
}


// ************************************************************************* //
