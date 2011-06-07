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
    PriorityArray

Description
    Priority array with changing priorities for inserted elements.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

References:
    R. Sedgewick, "Algorithms", Addison-Wesley, Reading, MA, 2nd edition, 1988.

SourceFiles
    PriorityArray.C

\*---------------------------------------------------------------------------*/

#ifndef PriorityArray_H
#define PriorityArray_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PriorityArray Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Compare = std::less<Type> >
class PriorityArray
{
    // Private data

        //- Curent heap size
        label heapSize_;

        //- Is the array heaped?
        bool heaped_;

        //- Heap indices
        labelList heap_;

        //- Heaped indices
        labelList from_;

        //- Heap weights
        List<Type> weights_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PriorityArray(const PriorityArray<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const PriorityArray<Type>&);


        //- Build heap
        void buildHeap();

        void heapify(label i);

        void fixHeapUp(label i);


public:

    // Constructors

        //- Construct given capacity
        PriorityArray(const label capacity);


    // Destructor - default


    // Member Functions

        //- Return heap size
        label heapSize() const
        {
            return heapSize_;
        }

        //- Is the array empty?
        bool empty() const
        {
            return heapSize_ <= 0;
        }

        //- Return indices
        const labelList& heap() const
        {
            return heap_;
        }

        //- Return weights
        const List<Type>& weights() const
        {
            return weights_;
        }


        //- Return top index
        label top();

        //- Pop element
        void pop();

        //- Set element with weight
        void set(const label i, const Type& value);

        //- Increment weight
        void increment(const label i, const Type& delta);

        //- Return weight
        const Type& weight(const label i) const
        {
            return weights_[i];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PriorityArray.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
