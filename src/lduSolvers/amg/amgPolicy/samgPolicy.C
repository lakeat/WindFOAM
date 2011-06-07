/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    samgPolicy

Description
    Selective AMG policy

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "samgPolicy.H"
#include "crMatrix.H"
#include "amgMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGInterfaceField.H"
#include "PriorityArray.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(samgPolicy, 0);

    addToRunTimeSelectionTable(amgPolicy, samgPolicy, matrix);

} // End namespace Foam


const Foam::scalar Foam::samgPolicy::epsilon_ = 0.25;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::crMatrix> Foam::samgPolicy::filterMatrix() const
{
    // Get addressing
    const label nRows = matrix_.lduAddr().size();
    const unallocLabelList& col = matrix_.lduAddr().upperAddr();
    const unallocLabelList& row = matrix_.lduAddr().ownerStartAddr();

    // Get off-diagonal matrix coefficients
    const scalarField& diag = matrix_.diag();
    const scalarField& upper = matrix_.upper();
    const scalarField& lower = matrix_.upper();

    // Find largest negative coefficient in each row
    // Here, "negative" is defined as having sign opposite of that of diagonal
    scalarField strongCoeff(nRows, 0);

    // Filter for strong coefficients

    for (label i = 0; i < nRows; i++)
    {
        scalar signDiag = sign(diag[i]);

        for (label ip = row[i]; ip < row[i + 1]; ip++)
        {
            label j = col[ip];

            // Do row coefficient
            scalar magAij = mag(min(signDiag*upper[ip], 0));

            if (magAij > strongCoeff[i])
            {
                strongCoeff[i] = magAij;
            }

            // Do col coefficient
            scalar magAji = mag(min(signDiag*lower[ip], 0));

            if (magAji > strongCoeff[j])
            {
                strongCoeff[j] = magAji;
            }
        }
    }

    // Count coefficients in row stronger than fraction of strongest
    labelList count(nRows, 0);

    for (label i = 0; i < nRows; i++)
    {
        scalar signDiag = sign(diag[i]);

        for (label ip = row[i]; ip < row[i + 1]; ip++)
        {
            label j = col[ip];

            // Do row coefficient
            scalar magAij = mag(min(signDiag*upper[ip], 0));

            if (magAij > epsilon_*strongCoeff[i])
            {
                count[i]++;
            }

            scalar magAji = mag(min(signDiag*lower[ip], 0));

            // Do col coefficient
            if (magAji > epsilon_*strongCoeff[j])
            {
                count[j]++;
            }
        }
    }

    // Create filtered matrix and addressing

    tmp<crMatrix> tS(new crMatrix(nRows, nRows, count));
    crMatrix& S = tS();

    // Set column index and coefficients
    const labelList& sRow = S.crAddr().row();
    labelList& sCol = S.col();
    scalarField& sCoeffs = S.coeffs();

    // Reset count for re-use
    count = 0;

    for (label i = 0; i < nRows; i++)
    {
        scalar signDiag = sign(diag[i]);

        for (label ip = row[i]; ip < row[i + 1]; ip++)
        {
            label j = col[ip];

            scalar magAij = mag(min(signDiag*upper[ip], 0));

            if (magAij > epsilon_*strongCoeff[i])
            {
                sCol[sRow[i] + count[i]] = j;
                sCoeffs[sRow[i] + count[i]] = upper[ip];
                count[i]++;
            }

            scalar magAji = mag(min(signDiag*lower[ip], 0));

            if (magAji > epsilon_*strongCoeff[j])
            {
                sCol[sRow[j] + count[j]] = i;
                sCoeffs[sRow[j] + count[j]] = lower[ip];
                count[j]++;
            }
        }
    }

    return tS;
}


Foam::tmp<Foam::labelField> Foam::samgPolicy::markEquations(const crMatrix& S)
{
    // Get addressing
    const label sNRows = S.crAddr().nRows();
    const labelList& sRow = S.crAddr().row();
    const labelList& sCol = S.crAddr().col();

    // Get transpose addressing of S
    crAddressing Taddr = S.crAddr().T();
    const labelList& tRow = Taddr.row();
    const labelList& tCol = Taddr.col();

    // Create mark array
    tmp<labelField> tmark(new labelField(sNRows, UNDECIDED));
    labelField& mark = tmark();

    // Set initial cardinality
    PriorityArray<label, std::greater<label> > cardinality(sNRows);

    for (label i = 0; i < sNRows; i++)
    {
        cardinality.set(i, sRow[i + 1] - sRow[i]);
    }

    // Mark highly diagonally dominant rows as fine
    for (label i = 0; i < sNRows; i++)
    {
        if (sRow[i + 1] == sRow[i])
        {
            mark[i] = FINE;
        }
    }

    // Start counting coarse equations
    nCoarseEqns_ = 0;

    while (!cardinality.empty())
    {
        label i = cardinality.top();
        cardinality.pop();

        if (mark[i] == UNDECIDED)
        {
            // Make highest cardinality equation coarse, decrement cardinality
            mark[i] = nCoarseEqns_;
            nCoarseEqns_++;

            for (label ip = sRow[i]; ip < sRow[i + 1]; ip++)
            {
                label j = sCol[ip];

                if (mark[j] == UNDECIDED)
                {
                    cardinality.increment(j, -1);
                }
            }

            // Make all neighbours fine and increment cardinality
            for (label ip = tRow[i]; ip < tRow[i + 1]; ip++)
            {
                label j = tCol[ip];

                if (mark[j] == UNDECIDED)
                {
                    mark[j] = FINE;

                    for (label jp = sRow[j]; jp < sRow[j + 1]; jp++)
                    {
                        label k = sCol[jp];

                        if (mark[k] == UNDECIDED)
                        {
                            cardinality.increment(k, 2);
                        }
                    }
                }
            }
        }
    }

    return tmark;
}


void Foam::samgPolicy::calcRP()
{
    // Get addressing
    const label nRows = matrix_.lduAddr().size();
    const unallocLabelList& col = matrix_.lduAddr().upperAddr();
    const unallocLabelList& row = matrix_.lduAddr().ownerStartAddr();

    // Get matrix coefficients
    const scalarField& diag = matrix_.diag();
    const scalarField& upper = matrix_.upper();
    const scalarField& lower = matrix_.upper();

    // Specify natural ordering of equations
    forAll (eqnOrder_, i)
    {
        eqnOrder_[i] = i;
    }

    // Grab filter matrix
    const crMatrix S(filterMatrix());

    // Mark equations and calculate coarse size
    const labelField mark(markEquations(S));

    // Assemble prolongation addressing
    const labelList& sRow = S.crAddr().row();
    const labelList& sCol = S.crAddr().col();
    const scalarField& sCoeffs = S.coeffs();

    // Cannot create P addressing in first pass.  Count and resize arrays
    label maxPCount = 0;

    forAll (mark, i)
    {
        if (mark[i] == FINE)
        {
            // Interpolation involves neighbourhood
            maxPCount += sRow[i + 1] - sRow[i];
        }
        else
        {
            // Coarse point, injection
            maxPCount++;
        }
    }

    // Assemble row, column and coefficients

    // Note: numerator and denominator will be assembled for the matrix
    // on the fly.  Each coefficient will add its contributions to the
    // neighbours.  By the time the row is reached, all lower triangular
    // contributions will be added.  It remains to do the upper triangle
    // and complete the work for the row.
    // This increases the storage for num and den to complete fields

    // Memory management
    {
        labelList pRow(nRows + 1);
        labelList pCol(maxPCount, 0);
        scalarField pCoeffs(maxPCount, 0);

        scalarField num(nRows, 0);
        scalarField Dii(sign(diag)*diag);

        scalar Dij, Dji, den, signDii;

        // Start row assembly
        pRow[0] = 0;

        for (label i = 0; i < nRows; i++)
        {
            label np = pRow[i];

            // Handle the unknown sign of diagonal: multiplying
            // row coeffs by the sign
            signDii = sign(diag[i]);

            for (label ip = row[i]; ip < row[i + 1]; ip++)
            {
                // Add to local coeff
                Dij = signDii*upper[ip];

                num[i] += Foam::min(Dij, 0);
                Dii[i] += Foam::max(Dij, 0);

                // Scatter lower triangular contribution
                label j = col[ip];

                Dji = sign(diag[j])*lower[ip];
                num[j] += Foam::min(Dji, 0);
                Dii[j] += Foam::max(Dji, 0);
            }

            // Row i completed.  Calculate weights
            if (mark[i] == FINE)
            {
                // Fine equation
                den = 0;

                for (label sip = sRow[i]; sip < sRow[i + 1]; sip++)
                {
                    label js = sCol[sip];

                    if (mark[js] != FINE)
                    {
                        den += sCoeffs[sip];
                    }
                }

                for (label sip = sRow[i]; sip < sRow[i + 1]; sip++)
                {
                    label js = sCol[sip];

                    if (mark[js] != FINE)
                    {
                        pCoeffs[np] = -(num[i]/den)*sCoeffs[sip]/Dii[i];
                        pCol[np] = mark[js];
                        np++;
                    }
                }
            }
            else
            {
                // Coarse equation
                pCoeffs[np] = 1;
                pCol[np] = mark[i];
                np++;
            }

            // Grab row start
            pRow[i + 1] = np;
        }

        // Resize column and coeffs
        pCoeffs.setSize(pRow[nRows]);
        pCol.setSize(pRow[nRows]);

        // Create prolongation matrix
        Pptr_ = new crMatrix(nRows, nCoarseEqns_, pRow, pCol);
        crMatrix& P = *Pptr_;

        // Copy coefficients
        P.coeffs().transfer(pCoeffs);
    } // End memory management

    // RUBBISH TO HERE!!!

    // Create ordering for smoothing sweeps
    label nOrder = 0;

    // All coarse first
    for (label i = 0; i < nRows; i++)
    {
        if (mark[i] != FINE)
        {
            eqnOrder_[nOrder] = i;
            nOrder++;
        }
    }

    // Followed by all fine
    for (label i = 0; i < nRows; i++)
    {
        if (mark[i] == FINE)
        {
            eqnOrder_[nOrder] = i;
            nOrder++;
        }
    }

    // Create restriction matrix

    if (nCoarseEqns_ > minCoarseEqns() && 3*nCoarseEqns_ <= 2*nRows)
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (coarsen_)
    {
        // Coarsening OK, make restriction matrix
        Rptr_ = new crMatrix(Pptr_->T());
    }
    else
    {
        // Coarsening did not succeed.  Delete Pptr
        deleteDemandDrivenData(Pptr_);
    }        
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::samgPolicy::samgPolicy
(
    const lduMatrix& matrix,
    const label groupSize,
    const label minCoarseEqns
)
:
    amgPolicy(groupSize, minCoarseEqns),
    matrix_(matrix),
    groupSize_(groupSize),
    nCoarseEqns_(0),
    coarsen_(false),
    Pptr_(NULL),
    Rptr_(NULL),
    eqnOrder_(matrix_.lduAddr().size())
{
    calcRP();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::samgPolicy::~samgPolicy()
{
    deleteDemandDrivenData(Rptr_);
    deleteDemandDrivenData(Pptr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::amgMatrix> Foam::samgPolicy::restrictMatrix
(
    const FieldField<Field, scalar>& bouCoeffs,
    const FieldField<Field, scalar>& intCoeffs,
    const lduInterfaceFieldPtrsList& interfaceFields
) const
{
    if (!coarsen_)
    {
        FatalErrorIn("autoPtr<amgMatrix> samgPolicy::restrictMatrix() const")
            << "Requesting coarse matrix when it cannot be created"
            << abort(FatalError);
    }

    // Get references to restriction and prolongation matrix

    const crMatrix& R = *Rptr_;
    const crMatrix& P = *Pptr_;

    // Get connectivity
    const crAddressing& crR = R.crAddr();
    const crAddressing& crP = P.crAddr();

#   ifdef FULLDEBUG
    // Check sized chain rule
    const label nEqns = matrix_.lduAddr().size();

    if
    (
        crR.nCols() != nEqns
     || crP.nRows() != nEqns
     || crR.nRows() != nCoarseEqns_
    )
    {
        FatalErrorIn("")
            << "Incompatible matrices for triple product: "
            << "R( " << crR.nRows() << " ," << crR.nCols() << ") "
            << "A( " << nEqns << " ," << nEqns << ") "
            << "P( " << crP.nRows() << " ," << crP.nCols() << ") "
            << abort(FatalError);
    }
#   endif

    const labelList& rRow = crR.row();
    const labelList& rCol = crR.col();

    const unallocLabelList& aRow = matrix_.lduAddr().ownerStartAddr();
    const unallocLabelList& aCol = matrix_.lduAddr().upperAddr();

    // Addressing for lower triangle loop
    const unallocLabelList& aOwn = matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losort = matrix_.lduAddr().losortAddr();
    const unallocLabelList& aLsrt = matrix_.lduAddr().losortStartAddr();

    const labelList& pRow = crP.row();
    const labelList& pCol = crP.col();


    // Calculate the addressing for the coarse matrix

    // Notes:
    //   - Addressing is always symmetrical and will be calculated
    //     using only the upper triangle
    //   - First loop counts coefficients

    labelList mark(nCoarseEqns_, -1);

    label nCoarseCoeffs = 0;

    for (label ia = 0; ia < crR.nRows(); ia++)
    {
        for (label ja = rRow[ia]; ja < rRow[ia + 1]; ja++)
        {
            const label jpa = rCol[ja];

            // Off-diagonal from fine off-diagonal, lower triangle
            for (label jbl = aLsrt[jpa]; jbl < aLsrt[jpa + 1]; jbl++)
            {
                // Get coefficient index
                const label jb = losort[jbl];

                // Get column index
                const label jpb = aOwn[jb];

                for (label jc = pRow[jpb]; jc < pRow[jpb + 1]; jc++)
                {
                    const label jpc = pCol[jc];

                    if (ia > jpc)
                    {
                        // Lower triangle, not stored
                        continue;
                    }
                    else if (ia != jpc)
                    {
                        // Found coarse off-diagonal contribution.
                        // Check if it exists and if not add the coefficient
                        if (mark[jpc] != ia)
                        {
                            mark[jpc] = ia;
                            nCoarseCoeffs++;
                        }
                    }
                }
            }

            // Off-diagonal from fine off-diagonal, upper triangle
            for (label jb = aRow[jpa]; jb < aRow[jpa + 1]; jb++)
            {
                const label jpb = aCol[jb];

                for (label jc = pRow[jpb]; jc < pRow[jpb + 1]; jc++)
                {
                    const label jpc = pCol[jc];

                    if (ia > jpc)
                    {
                        // Lower triangle, not stored
                        continue;
                    }
                    else if (ia != jpc)
                    {
                        // Found coarse off-diagonal contribution.
                        // Check if it exists and if not add the coefficient
                        if (mark[jpc] != ia)
                        {
                            mark[jpc] = ia;
                            nCoarseCoeffs++;
                        }
                    }
                }
            }

            // Off-diagonal from fine diagonal
            for (label jc = pRow[jpa]; jc < pRow[jpa + 1]; jc++)
            {
                const label jpc = pCol[jc];

                if (ia > jpc)
                {
                    // Lower triangle, not stored
                    continue;
                }
                else if (ia != jpc)
                {
                    // Found coarse off-diagonal contribution.
                    // Check if it exists and if not add the coefficient
                    if (mark[jpc] != ia)
                    {
                        mark[jpc] = ia;
                        nCoarseCoeffs++;
                    }
                }
            }
        }
    }

    // Make coarse matrix
    labelList coarseOwner(nCoarseCoeffs);
    labelList coarseNeighbour(nCoarseCoeffs);
    scalarField coarseUpper(nCoarseCoeffs, 0);
    scalarField coarseDiag(nCoarseEqns_, 0);

    // Insert upper coefficients
    const scalarField& rCoeffs = R.coeffs();
    const scalarField& aDiag = matrix_.diag();
    const scalarField& aUpper = matrix_.upper();
    const scalarField& pCoeffs = P.coeffs();

    // Re-initialise mark vector
    mark = -1;
    nCoarseCoeffs = 0;

    for (label ia = 0; ia < crR.nRows(); ia++)
    {
        label rsCoarseCoeffs = nCoarseCoeffs;

        for (label ja = rRow[ia]; ja < rRow[ia + 1]; ja++)
        {
            const label jpa = rCol[ja];

            // Off-diagonal from fine off-diagonal, lower triangle
            for (label jbl = aLsrt[jpa]; jbl < aLsrt[jpa + 1]; jbl++)
            {
                // Get coefficient index
                const label jb = losort[jbl];

                // Get column index
                const label jpb = aOwn[jb];

                const scalar ab = rCoeffs[ja]*aUpper[jb];

                for (label jc = pRow[jpb]; jc < pRow[jpb + 1]; jc++)
                {
                    const label jpc = pCol[jc];

                    if (ia > jpc)
                    {
                        // Lower triangle, not stored
                        continue;
                    }
                    else if (ia == jpc)
                    {
                        // Found diagonal contribution from off-diagonal
                        coarseDiag[ia] += ab*pCoeffs[jc];
                    }
                    else
                    {
                        // Found coarse off-diagonal contribution.
                        // Check if it exists and if not add the coefficient
                        label jp = mark[jpc];

                        if (jp == -1)
                        {
                            jp = nCoarseCoeffs;
                            mark[jpc] = nCoarseCoeffs;

                            // Grab addressing
                            coarseOwner[nCoarseCoeffs] = ia;
                            coarseNeighbour[nCoarseCoeffs] = jpc;

                            nCoarseCoeffs++;
                        }

                        coarseUpper[jp] += ab*pCoeffs[jc];
                    }
                }
            }

            // Off-diagonal from fine off-diagonal, upper triangle
            for (label jb = aRow[jpa]; jb < aRow[jpa + 1]; jb++)
            {
                const label jpb = aCol[jb];

                const scalar ab = rCoeffs[ja]*aUpper[jb];

                for (label jc = pRow[jpb]; jc < pRow[jpb + 1]; jc++)
                {
                    const label jpc = pCol[jc];

                    if (ia > jpc)
                    {
                        // Lower triangle, not stored
                        continue;
                    }
                    else if (ia == jpc)
                    {
                        // Found diagonal contribution from off-diagonal
                        coarseDiag[ia] += ab*pCoeffs[jc];
                    }
                    else
                    {
                        // Found coarse off-diagonal contribution.
                        // Check if it exists and if not add the coefficient
                        label jp = mark[jpc];

                        if (jp == -1)
                        {
                            jp = nCoarseCoeffs;
                            mark[jpc] = nCoarseCoeffs;

                            // Grab addressing
                            coarseOwner[nCoarseCoeffs] = ia;
                            coarseNeighbour[nCoarseCoeffs] = jpc;

                            nCoarseCoeffs++;
                        }

                        coarseUpper[jp] += ab*pCoeffs[jc];
                    }
                }
            }

            // Off-diagonal from fine diagonal
            const scalar ab = rCoeffs[ja]*aDiag[jpa];

            for (label jc = pRow[jpa]; jc < pRow[jpa + 1]; jc++)
            {
                const label jpc = pCol[jc];

                if (ia > jpc)
                {
                    // Lower triangle, not stored
                    continue;
                }
                else if (ia == jpc)
                {
                    coarseDiag[ia] += ab*pCoeffs[jc];
                }
                else
                {
                    // Found coarse off-diagonal contribution.
                    // Check if it exists and if not add the coefficient
                    label jp = mark[jpc];

                    if (jp == -1)
                    {
                        jp = nCoarseCoeffs;
                        mark[jpc] = nCoarseCoeffs;

                        // Grab addressing
                        coarseOwner[nCoarseCoeffs] = ia;
                        coarseNeighbour[nCoarseCoeffs] = jpc;

                        nCoarseCoeffs++;
                    }

                    coarseUpper[jp] += ab*pCoeffs[jc];
                }
            }
        }

        // Reset mark of the use in next row
        for (label j = rsCoarseCoeffs; j < nCoarseCoeffs; j++)
        {
            mark[coarseNeighbour[j]] = -1;
        }
    }


    // Create coarse-level coupled interface fields.  HJ, hacked!

    // Set the coarse interfaces and coefficients
    lduInterfacePtrsList* coarseInterfacesPtr =
        new lduInterfacePtrsList(interfaceFields.size());
    lduInterfacePtrsList& coarseInterfaces = *coarseInterfacesPtr;

    // Set the coarse interfaceFields and coefficients
    lduInterfaceFieldPtrsList* coarseInterfaceFieldsPtr =
        new lduInterfaceFieldPtrsList(interfaceFields.size());
    lduInterfaceFieldPtrsList& coarseInterfaceFields =
        *coarseInterfaceFieldsPtr;

    FieldField<Field, scalar>* coarseBouCoeffsPtr =
        new FieldField<Field, scalar>(interfaceFields.size());
    FieldField<Field, scalar>& coarseBouCoeffs = *coarseBouCoeffsPtr;

    FieldField<Field, scalar>* coarseIntCoeffsPtr =
        new FieldField<Field, scalar>(interfaceFields.size());
    FieldField<Field, scalar>& coarseIntCoeffs = *coarseIntCoeffsPtr;

    labelListList coarseInterfaceAddr(interfaceFields.size());

    // Add the coarse level: This code is not operational

    // Set the coarse ldu addressing onto the list
    lduPrimitiveMesh* coarseAddrPtr =
        new lduPrimitiveMesh
        (
            nCoarseEqns_,
            coarseOwner,
            coarseNeighbour,
            true
        );

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields[intI].interface();

            // Warning: dummy coarsening for SAMG: code missing!
            coarseInterfaces.set
            (
                intI,
                GAMGInterface::New
                (
                    *coarseAddrPtr,
                    fineInterface,
                    labelField::null(),
                    labelField::null()
                ).ptr()
            );
        }
    }

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const GAMGInterface& coarseInterface = 
                refCast<const GAMGInterface>(coarseInterfaces[intI]);

            coarseInterfaceFields.set
            (
                intI,
                GAMGInterfaceField::New
                (
                    coarseInterface,
                    interfaceFields[intI]
                ).ptr()
            );

            coarseBouCoeffs.set
            (
                intI,
                coarseInterface.agglomerateCoeffs(bouCoeffs[intI])
            );

            coarseIntCoeffs.set
            (
                intI,
                coarseInterface.agglomerateCoeffs(intCoeffs[intI])
            );

            coarseInterfaceAddr[intI] = coarseInterface.faceCells();
        }
    }

    // Add interfaces
    coarseAddrPtr->addInterfaces
    (
        *coarseInterfacesPtr,
        coarseInterfaceAddr,
        matrix_.patchSchedule()
    );

    // Matrix restriction done!

    // Set the coarse level matrix
    lduMatrix* coarseMatrixPtr = new lduMatrix(*coarseAddrPtr);
    lduMatrix& coarseMatrix = *coarseMatrixPtr;

    // Insert coarse matrix coefficients
    coarseMatrix.upper().transfer(coarseUpper);
    coarseMatrix.diag().transfer(coarseDiag);

    // Create and return amgMatrix
    return autoPtr<amgMatrix>
    (
        new amgMatrix
        (
            coarseAddrPtr,
            coarseInterfacesPtr,
            coarseMatrixPtr,
            coarseBouCoeffsPtr,
            coarseIntCoeffsPtr,
            coarseInterfaceFieldsPtr
        )
    );
}


void Foam::samgPolicy::restrictResidual
(
    const scalarField& res,
    scalarField& coarseRes
) const
{
    coarseRes = 0;
    Rptr_->dotPlus(coarseRes, res);
}


void Foam::samgPolicy::prolongateCorrection
(
    scalarField& x,
    const scalarField& coarseX
) const
{
    Pptr_->dotPlus(x, coarseX);

}


// ************************************************************************* //
