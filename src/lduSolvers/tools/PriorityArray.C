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

\*----------------------------------------------------------------------------*/

#include "PriorityArray.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class Compare>
void Foam::PriorityArray<Type, Compare>::buildHeap()
{
    for (label i = 0; i < heapSize_; i++)
    {
        heap_[i] = i;
        from_[i] = i;
    }

    for (label i = heapSize_/2 - 1; i >= 0; i--)
    {
        this->heapify(i);
    }

    heaped_ = true;
}


template<class Type, class Compare>
void Foam::PriorityArray<Type, Compare>::heapify(label i)
{
    Compare cmp;

    for (;;)
    {
        // Find largest node and left and right children
        label n = i;
        label il = 2*i + 1;
        label ir = il + 1;

        if (il < heapSize_ && cmp(weights_[heap_[il]], weights_[heap_[n]]))
        {
            n = il;
        }

        if (ir < heapSize_ && cmp(weights_[heap_[ir]], weights_[heap_[n]]))
        {
            n = ir;
        }

        if (n == i) break;

        // Swap i with largest n
        Foam::Swap(heap_[i], heap_[n]);

        // Update inverse positions
        from_[heap_[i]] = i;
        from_[heap_[n]] = n;

        // Make child n into a heap
        i = n;
    }
}


template<class Type, class Compare>
void Foam::PriorityArray<Type, Compare>::fixHeapUp(label i)
{
    Compare cmp;
    int n = (i - 1)/2;

    while (i > 0 && cmp(weights_[heap_[i]], weights_[heap_[n]]))
    {
        // Swap node i and n
        Foam::Swap(heap_[i], heap_[n]);

        // Update inverse positions
        from_[heap_[i]] = i;
        from_[heap_[n]] = n;

        i = n;
        n = (i - 1)/2;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given capacity
template<class Type, class Compare>
Foam::PriorityArray<Type, Compare>::PriorityArray(const label capacity)
:
    heapSize_(capacity),
    heaped_(false),
    heap_(capacity),
    from_(capacity),
    weights_(capacity)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class Compare>
Foam::label Foam::PriorityArray<Type, Compare>::top()
{
    if (!heaped_)
    {
        this->buildHeap();
    }

    return heap_[0];
}


template<class Type, class Compare>
void Foam::PriorityArray<Type, Compare>::pop()
{
    if (!heaped_)
    {
        this->buildHeap();
    }

    if (--heapSize_ > 0)
    {
        Foam::Swap(heap_[0], heap_[heapSize_]);
        from_[heap_[0]] = 0;
        this->heapify(0);
    }
}


template<class Type, class Compare>
void Foam::PriorityArray<Type, Compare>::set(const label i, const Type& value)
{
    heaped_ = false;
    weights_[i] = value;
}


template<class Type, class Compare>
void Foam::PriorityArray<Type, Compare>::increment
(
    const label i,
    const Type& delta
)
{
    Compare cmp;

    if (!heaped_)
    {
        this->buildHeap();
    }

    weights_[i] += delta;

    if (cmp(delta, pTraits<Type>::zero))
    {
        this->fixHeapUp(from_[i]);
    }
    else
    {
        this->heapify(from_[i]);
    }
}


// ************************************************************************* //
