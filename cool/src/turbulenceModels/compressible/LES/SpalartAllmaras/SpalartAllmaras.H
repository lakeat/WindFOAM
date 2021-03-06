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
    Foam::compressible::LESModels::SpalartAllmaras

Description
    SpalartAllmaras for compressible flows

SourceFiles
    SpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleSpalartAllmaras_H
#define compressibleSpalartAllmaras_H

#include "LESModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                       Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmaras
:
    public LESModel
{
    // Private data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar Prt_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cv2_;
            dimensionedScalar CDES_;
            dimensionedScalar ck_;
            dimensionedScalar kappa_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;


    // Fields

        volScalarField nuTilda_;
        volScalarField dTilda_;
        volScalarField muSgs_;
        volScalarField alphaSgs_;


    // Private member functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields();

        tmp<volScalarField> fv1() const;
        tmp<volScalarField> fv2() const;
        tmp<volScalarField> fv3() const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;

        // Disallow default bitwise copy construct and assignment
        SpalartAllmaras(const SpalartAllmaras&);
        SpalartAllmaras& operator=(const SpalartAllmaras&);


public:

    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- Constructor from components
        SpalartAllmaras
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermoPhysicalModel
        );


    //- Destructor
    virtual ~SpalartAllmaras()
    {}


    // Member Functions

        tmp<volScalarField> nuTilda() const
        {
            return nuTilda_;
        }

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return sqr(muSgs()/rho()/ck_/dTilda_);
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return SGS viscosity
        virtual tmp<volScalarField> muSgs() const
        {
            return muSgs_;
        }

        //- Return SGS thermal diffusivity
        virtual tmp<volScalarField> alphaSgs() const
        {
            return alphaSgs_;
        }

        //- Return thermal conductivity
        virtual tmp<volScalarField> alphaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("alphaEff", alphaSgs_ + alpha())
            );
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoBeff() const;

        //- Returns div(rho*dev(B)).
        // This is the additional term due to the filtering of the NSE.
        virtual tmp<fvVectorMatrix> divDevRhoBeff(volVectorField& U) const;

        //- Correct nuTilda and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
