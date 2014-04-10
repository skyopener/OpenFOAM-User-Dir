/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::interfaceVPProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfaceVPProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceVPProperties_H
#define interfaceVPProperties_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceVPProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceVPProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Keep a reference to the VPProperties dictionary
        const dictionary& VPPropertyDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension
        dimensionedScalar sigma_;

        //- Contact angle for VP
        scalar thetaVP_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& alpha1_;
        const volScalarField& voidFraction_;
        const volVectorField& U_;
        surfaceScalarField nHatf_;
        volScalarField K_;

        dimensionedScalar cVoid_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfaceVPProperties(const interfaceVPProperties&);
        void operator=(const interfaceVPProperties&);

        //- Correction for the solid boundary condition: related to VP method
        void correctSolidContactAngle
        (
            surfaceVectorField& nHat,
            const surfaceVectorField& gradAlphaf
        ) const;

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::GeometricBoundaryField& nHat,
            surfaceVectorField::GeometricBoundaryField& gradAlphaf
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceVPProperties
        (
            const volScalarField& alpha1,
            const volScalarField& voidFraction,
            const volVectorField& U,
            const IOdictionary&,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        const volScalarField& K() const
        {
            return K_;
        }

        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

        tmp<volScalarField> sigmaK() const
        {
            return sigma_*K_;
        }

        void correct()
        {
            calculateK();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
