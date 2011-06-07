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
    aamgPolicy

Description
    Agglomerative AMG policy

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "aamgPolicy.H"
#include "amgMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "GAMGInterfaceField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aamgPolicy, 0);

    addToRunTimeSelectionTable(amgPolicy, aamgPolicy, matrix);

} // End namespace Foam


const Foam::scalar Foam::aamgPolicy::weightFactor_ = 0.65;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::aamgPolicy::calcChild()
{
    // Algorithm:
    // 1) Create temporary equation addressing using a double-pass algorithm.
    //    to create the offset table.
    // 2) Loop through all equations and for each equation find the best fit
    //    neighbour.  If all neighbours are grouped, add equation to best group

    // Create row-based addressing

    const label nRows = matrix_.lduAddr().size();

    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

    // For each equation get number of coefficients in a row
    labelList cols(upperAddr.size() + lowerAddr.size());
    labelList cIndex(upperAddr.size() + lowerAddr.size());
    labelList rowOffset(nRows + 1, 0);

    // Count the number of coefficients
    forAll (upperAddr, coeffI)
    {
        rowOffset[upperAddr[coeffI]]++;
    }

    forAll (lowerAddr, coeffI)
    {
        rowOffset[lowerAddr[coeffI]]++;
    }

    label nCoeffs = 0;

    forAll (rowOffset, eqnI)
    {
        nCoeffs += rowOffset[eqnI];
    }

    rowOffset[nRows] = nCoeffs;

    for (label eqnI = nRows - 1; eqnI >= 0; --eqnI)
    {
        rowOffset[eqnI] = rowOffset[eqnI  + 1] - rowOffset[eqnI];
    }

    rowOffset[0] = 0;

    // Create column and coefficient index array
    {
        // Use child array to count number of entries per row.
        // Reset the list for counting
        labelList& nPerRow = child_;
        nPerRow = 0;

        forAll (upperAddr, coeffI)
        {
            cols[rowOffset[upperAddr[coeffI]] + nPerRow[upperAddr[coeffI]]] =
                lowerAddr[coeffI];

            cIndex[rowOffset[upperAddr[coeffI]] + nPerRow[upperAddr[coeffI]]] =
                coeffI;

            nPerRow[upperAddr[coeffI]]++;
        }

        forAll (lowerAddr, coeffI)
        {
            cols[rowOffset[lowerAddr[coeffI]] + nPerRow[lowerAddr[coeffI]]] =
                upperAddr[coeffI];

            cIndex[rowOffset[lowerAddr[coeffI]] + nPerRow[lowerAddr[coeffI]]] =
                coeffI;

            nPerRow[lowerAddr[coeffI]]++;
        }

        // Reset child array
        child_ = -1;
    }


    // Calculate agglomeration

    // Get matrix coefficients
    const scalarField& diag = matrix_.diag();
    const scalarField& upper = matrix_.upper();

    labelList sizeOfGroups(nRows, 0);

    nCoarseEqns_ = 0;

    label indexU, indexG, colI, curEqn, nextEqn, groupPassI;

    scalar magRowDiag, magColDiag, weight, weightU, weightG;

    for (label eqnI = 0; eqnI < nRows; eqnI++)
    {
        if (child_[eqnI] == -1)
        {
            curEqn = eqnI;

            indexU = -1;
            indexG = -1;

            child_[curEqn] = nCoarseEqns_;

            magRowDiag = Foam::mag(diag[curEqn]);

            for (groupPassI = 1; groupPassI < groupSize_; groupPassI++)
            {
                weightU = 0;
                weightG = 0;

                indexU = -1;
                indexG = -1;

                for
                (
                    label rowCoeffI = rowOffset[curEqn];
                    rowCoeffI < rowOffset[curEqn + 1];
                    rowCoeffI++
                )
                {
                    colI = cols[rowCoeffI];

                    magColDiag = Foam::mag(diag[colI]);

                    weight = Foam::mag(upper[cIndex[rowCoeffI]])/
                        max(magRowDiag, magColDiag);

                    if (child_[colI] == -1)
                    {
                        if (indexU == -1 || weight > weightU)
                        {
                            indexU = rowCoeffI;
                            weightU = weight;
                        }
                    }
                    else if (child_[curEqn] != child_[colI])
                    {
                        if (indexG == -1 || weight > weightG)
                        {
                            indexG = rowCoeffI;
                            weightG = weight;
                        }
                    }
                }

                if
                (
                    indexU != -1
                 && (indexG == -1 || weightU >= weightFactor_*weightG)
                )
                {
                    // Found new element of group.  Add it and use as
                    // start of next search

                    nextEqn = cols[indexU];

                    child_[nextEqn] = child_[curEqn];
                    sizeOfGroups[child_[curEqn]];

                    curEqn = nextEqn;
                }
                else
                {
                    // Group full or cannot be extended
                    break;
                }
            }

            if
            (
                groupPassI > 1
             || indexG == -1
             || sizeOfGroups[child_[cols[indexG]]] > (groupSize_ + 2)
            )
            {
                sizeOfGroups[child_[eqnI]]++;
                nCoarseEqns_++;
            }
            else
            {
                child_[eqnI] = child_[cols[indexG]];
                sizeOfGroups[child_[cols[indexG]]]++;
            }
        }
    }

    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.

    if (nCoarseEqns_ > minCoarseEqns() && 3*nCoarseEqns_ <= 2*nRows)
    {
        coarsen_ = true;
    }

    reduce(coarsen_, andOp<bool>());

    if (lduMatrix::debug >= 2)
    {
        Pout << "Coarse level size: " << nCoarseEqns_;

        if (coarsen_)
        {
            Pout << ".  Accepted" << endl;
        }
        else
        {
            Pout << ".  Rejected" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aamgPolicy::aamgPolicy
(
    const lduMatrix& matrix,
    const label groupSize,
    const label minCoarseEqns
)
:
    amgPolicy(groupSize, minCoarseEqns),
    matrix_(matrix),
    child_(matrix_.lduAddr().size()),
    groupSize_(groupSize),
    nCoarseEqns_(0),
    coarsen_(false)
{
    calcChild();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aamgPolicy::~aamgPolicy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::amgMatrix> Foam::aamgPolicy::restrictMatrix
(
    const FieldField<Field, scalar>& bouCoeffs,
    const FieldField<Field, scalar>& intCoeffs,
    const lduInterfaceFieldPtrsList& interfaceFields
) const
{
    if (!coarsen_)
    {
        FatalErrorIn("autoPtr<amgMatrix> aamgPolicy::restrictMatrix() const")
            << "Requesting coarse matrix when it cannot be created"
            << abort(FatalError);
    }

    // Construct the coarse matrix and ldu addressing for the next level
    // Algorithm:
    // 1) Loop through all fine coeffs. If the child labels on two sides are
    //    different, this creates a coarse coeff. Define owner and neighbour
    //    for this coeff based on cluster IDs.
    // 2) Check if the coeff has been seen before. If yes, add the coefficient
    //    to the appropriate field (stored with the equation). If no, create
    //    a new coeff with neighbour ID and add the coefficient
    // 3) Once all the coeffs have been created, loop through all clusters and
    //    insert the coeffs in the upper order. At the same time, collect the
    //    owner and neighbour addressing.
    // 4) Agglomerate the diagonal by summing up the fine diagonal

    // Get addressing
    const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

    const label nFineCoeffs = upperAddr.size();

#   ifdef FULLDEBUG
    if (child_.size() != matrix_.lduAddr().size())
    {
        FatalErrorIn
        (
            "autoPtr<amgMatrix> aamgPolicy::restrictMatrix() const"
        )   << "Child array does not correspond to fine level. " << endl
            << " Child size: " << child_.size()
            << " number of equations: " << matrix_.lduAddr().size()
            << abort(FatalError);
    }
#   endif


    // Storage for block neighbours and coefficients

    // Guess initial maximum number of neighbours in block
    label maxNnbrs = 10;

    // Number of neighbours per block
    labelList blockNnbrs(nCoarseEqns_, 0);

    // Setup initial packed storage for neighbours and coefficients
    labelList blockNbrsData(maxNnbrs*nCoarseEqns_);

    // Create face-restriction addressing
    labelList coeffRestrictAddr(nFineCoeffs);

    // Initial neighbour array (not in upper-triangle order)
    labelList initCoarseNeighb(nFineCoeffs);

    // Counter for coarse coeffs
    label nCoarseCoeffs = 0;

    // Loop through all fine coeffs
    forAll (upperAddr, fineCoeffi)
    {
        label rmUpperAddr = child_[upperAddr[fineCoeffi]];
        label rmLowerAddr = child_[lowerAddr[fineCoeffi]];

        if (rmUpperAddr == rmLowerAddr)
        {
            // For each fine coeff inside of a coarse cluster keep the address
            // of the cluster corresponding to the coeff in the
            // coeffRestrictAddr as a negative index
            coeffRestrictAddr[fineCoeffi] = -(rmUpperAddr + 1);
        }
        else
        {
            // This coeff is a part of a coarse coeff

            label cOwn = rmUpperAddr;
            label cNei = rmLowerAddr;

            // Get coarse owner and neighbour
            if (rmUpperAddr > rmLowerAddr)
            {
                cOwn = rmLowerAddr;
                cNei = rmUpperAddr;
            }

            // Check the neighbour to see if this coeff has already been found
            bool nbrFound = false;
            label& ccnCoeffs = blockNnbrs[cOwn];

            for (int i = 0; i < ccnCoeffs; i++)
            {
                if (initCoarseNeighb[blockNbrsData[maxNnbrs*cOwn + i]] == cNei)
                {
                    nbrFound = true;
                    coeffRestrictAddr[fineCoeffi] =
                        blockNbrsData[maxNnbrs*cOwn + i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnCoeffs >= maxNnbrs)
                {
                    // Double the size of list and copy data
                    label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    // Resize and copy list
                    const labelList oldBlockNbrsData = blockNbrsData;
                    blockNbrsData.setSize(maxNnbrs*nCoarseEqns_);

                    forAll (blockNnbrs, i)
                    {
                        for (int j = 0; j < blockNnbrs[i]; j++)
                        {
                            blockNbrsData[maxNnbrs*i + j] =
                                oldBlockNbrsData[oldMaxNnbrs*i + j];
                        }
                    }
                }

                blockNbrsData[maxNnbrs*cOwn + ccnCoeffs] = nCoarseCoeffs;
                initCoarseNeighb[nCoarseCoeffs] = cNei;
                coeffRestrictAddr[fineCoeffi] = nCoarseCoeffs;
                ccnCoeffs++;

                // New coarse coeff created
                nCoarseCoeffs++;
            }
        }
    } // End for all fine coeffs


    // Renumber into upper-triangular order

    // All coarse owner-neighbour storage
    labelList coarseOwner(nCoarseCoeffs);
    labelList coarseNeighbour(nCoarseCoeffs);
    labelList coarseCoeffMap(nCoarseCoeffs);

    label coarseCoeffi = 0;

    forAll (blockNnbrs, cci)
    {
        label* cCoeffs = &blockNbrsData[maxNnbrs*cci];
        label ccnCoeffs = blockNnbrs[cci];

        for (int i = 0; i < ccnCoeffs; i++)
        {
            coarseOwner[coarseCoeffi] = cci;
            coarseNeighbour[coarseCoeffi] = initCoarseNeighb[cCoeffs[i]];
            coarseCoeffMap[cCoeffs[i]] = coarseCoeffi;
            coarseCoeffi++;
        }
    }

    forAll(coeffRestrictAddr, fineCoeffi)
    {
        if (coeffRestrictAddr[fineCoeffi] >= 0)
        {
            coeffRestrictAddr[fineCoeffi] =
                coarseCoeffMap[coeffRestrictAddr[fineCoeffi]];
        }
    }

    // Clear the temporary storage for the coarse matrix data
    blockNnbrs.setSize(0);
    blockNbrsData.setSize(0);
    initCoarseNeighb.setSize(0);
    coarseCoeffMap.setSize(0);


    // Create coarse-level coupled interfaces

    // Create coarse interfaces, addressing and coefficients

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

    // Add the coarse level

    // Set the coarse ldu addressing onto the list
    lduPrimitiveMesh* coarseAddrPtr =
        new lduPrimitiveMesh
        (
            nCoarseEqns_,
            coarseOwner,
            coarseNeighbour,
            true
        );

    // Initialise transfer of restrict addressing on the interface
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            interfaceFields[intI].interface().initInternalFieldTransfer
            (
                Pstream::blocking,
                child_
            );
        }
    }

    // Store coefficients to avoid tangled communications
    // HJ, 1/Apr/2009
    FieldField<Field, label> fineInterfaceAddr(interfaceFields.size());

    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields[intI].interface();

            fineInterfaceAddr.set
            (
                intI,
                new labelField
                (
                    fineInterface.internalFieldTransfer
                    (
                        Pstream::blocking,
                        child_
                    )
                )
            );
        }
    }

    // Create GAMG interfaces
    forAll (interfaceFields, intI)
    {
        if (interfaceFields.set(intI))
        {
            const lduInterface& fineInterface =
                interfaceFields[intI].interface();

            coarseInterfaces.set
            (
                intI,
                GAMGInterface::New
                (
                    *coarseAddrPtr,
                    fineInterface,
                    fineInterface.interfaceInternalField(child_),
                    fineInterfaceAddr[intI]
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

    // Coarse matrix diagonal initialised by restricting the
    // finer mesh diagonal
    scalarField& coarseDiag = coarseMatrix.diag();
    restrictResidual(matrix_.diag(), coarseDiag);

    // Check if matrix is assymetric and if so agglomerate both upper and lower
    // coefficients ...
    if (matrix_.hasLower())
    {
        // Get off-diagonal matrix coefficients
        const scalarField& fineUpper = matrix_.upper();
        const scalarField& fineLower = matrix_.lower();

        // Coarse matrix upper coefficients
        scalarField& coarseUpper = coarseMatrix.upper();
        scalarField& coarseLower = coarseMatrix.lower();

        forAll(coeffRestrictAddr, fineCoeffI)
        {
            label cCoeff = coeffRestrictAddr[fineCoeffI];

            if (cCoeff >= 0)
            {
                coarseUpper[cCoeff] += fineUpper[fineCoeffI];
                coarseLower[cCoeff] += fineLower[fineCoeffI];
            }
            else
            {
                // Add the fine face coefficients into the diagonal.
                coarseDiag[-1 - cCoeff] +=
                    fineUpper[fineCoeffI] + fineLower[fineCoeffI];
            }
        }
    }
    else // ... Otherwise it is symmetric so agglomerate just the upper
    {
        // Get off-diagonal matrix coefficients
        const scalarField& fineUpper = matrix_.upper();

        // Coarse matrix upper coefficients
        scalarField& coarseUpper = coarseMatrix.upper();

        forAll(coeffRestrictAddr, fineCoeffI)
        {
            label cCoeff = coeffRestrictAddr[fineCoeffI];

            if (cCoeff >= 0)
            {
                coarseUpper[cCoeff] += fineUpper[fineCoeffI];
            }
            else
            {
                // Add the fine face coefficient into the diagonal.
                coarseDiag[-1 - cCoeff] += 2*fineUpper[fineCoeffI];
            }
        }
    }

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


void Foam::aamgPolicy::restrictResidual
(
    const scalarField& res,
    scalarField& coarseRes
) const
{
    coarseRes = 0;

    forAll (res, i)
    {
        coarseRes[child_[i]] += res[i];
    }
}


void Foam::aamgPolicy::prolongateCorrection
(
    scalarField& x,
    const scalarField& coarseX
) const
{
    forAll (x, i)
    {
        x[i] += coarseX[child_[i]];
    }
}


// ************************************************************************* //
