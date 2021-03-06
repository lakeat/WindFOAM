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

#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Limiter, template<class> class LimitFunc>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::DeferredCorrectionLimitedScheme<Type, Limiter, LimitFunc>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Correction = full interpolation minus upwinded part
    return
        surfaceInterpolationScheme<Type>::interpolate
        (
            vf,
            limitedSurfaceInterpolationScheme<Type>::weights(vf)
        )
      - upwindScheme_.interpolate(vf);
}


template<class Type, class Limiter, template<class> class LimitFunc>
tmp<surfaceScalarField>
DeferredCorrectionLimitedScheme<Type, Limiter, LimitFunc>::limiter
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tLimiter
    (
        new surfaceScalarField
        (
            IOobject
            (
                type() + "Limiter(" + phi.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );
    surfaceScalarField& lim = tLimiter();

    tmp<GeometricField<typename Limiter::phiType, fvPatchField, volMesh> >
        tlPhi = LimitFunc<Type>()(phi);

    const GeometricField<typename Limiter::phiType, fvPatchField, volMesh>&
        lPhi = tlPhi();

    GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>
        gradc(fvc::grad(lPhi));

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& pLim = lim.internalField();

    forAll(pLim, face)
    {
        label own = owner[face];
        label nei = neighbour[face];

        pLim[face] = Limiter::limiter
        (
            CDweights[face],
            this->faceFlux_[face],
            lPhi[own],
            lPhi[nei],
            gradc[own],
            gradc[nei],
            C[nei] - C[own]
        );
    }

    surfaceScalarField::GeometricBoundaryField& bLim = lim.boundaryField();

    forAll (bLim, patchi)
    {
        scalarField& pLim = bLim[patchi];

        if (bLim[patchi].coupled())
        {
            const scalarField& pCDweights = CDweights.boundaryField()[patchi];
            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];
            Field<typename Limiter::phiType> plPhiP =
                lPhi.boundaryField()[patchi].patchInternalField();
            Field<typename Limiter::phiType> plPhiN =
                lPhi.boundaryField()[patchi].patchNeighbourField();
            Field<typename Limiter::gradPhiType> pGradcP =
                gradc.boundaryField()[patchi].patchInternalField();
            Field<typename Limiter::gradPhiType> pGradcN =
                gradc.boundaryField()[patchi].patchNeighbourField();

            // Build the d-vectors
            // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
            vectorField pd = bLim[patchi].patch().delta();

            forAll(pLim, face)
            {
                pLim[face] = Limiter::limiter
                (
                    pCDweights[face],
                    pFaceFlux[face],
                    plPhiP[face],
                    plPhiN[face],
                    pGradcP[face],
                    pGradcN[face],
                    pd[face]
                );
            }
        }
        else
        {
            pLim = 1.0;
        }
    }

    return tLimiter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
