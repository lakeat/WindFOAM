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

#include "processorFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void processorFvPatchField<scalar>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    procPatch_.compressedSend
    (
        commsType,
        patch().patchInternalField(psiInternal)()
    );
}


template<>
void processorFvPatchField<scalar>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    scalarField pnf
    (
        procPatch_.compressedReceive<scalar>(commsType, this->size())()
    );

    const unallocLabelList& faceCells = patch().faceCells();

    forAll(faceCells, facei)
    {
        result[faceCells[facei]] -= coeffs[facei]*pnf[facei];
    }
}


template<>
void processorFvPatchField<scalar>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const BlockLduMatrix<scalar>&,
    const CoeffField<scalar>&,
    const Pstream::commsTypes commsType
) const
{
    procPatch_.compressedSend
    (
        commsType,
        patch().patchInternalField(psiInternal)()
    );
}


template<>
void processorFvPatchField<scalar>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const BlockLduMatrix<scalar>&,
    const CoeffField<scalar>& coeffs,
    const Pstream::commsTypes commsType
) const
{  
    scalarField pnf
    (
        procPatch_.compressedReceive<scalar>(commsType, this->size())()
    );

    const unallocLabelList& faceCells = patch().faceCells();
    const scalarField& scalarCoeffs = coeffs.asScalar();
    
    forAll(faceCells, facei)
    {
        result[faceCells[facei]] -= scalarCoeffs[facei]*pnf[facei];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
