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

\*---------------------------------------------------------------------------*/

#include "Lambert.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Lambert<CloudType>::Lambert
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PendularModel<CloudType>(dict, cloud, typeName),
    surfaceTension_(readScalar(this->coeffDict().lookup("surfaceTension"))),
    contactAngle_(readScalar(this->coeffDict().lookup("contactAngle")))
{
    contactAngle_ = contactAngle_/180.*mathematical::pi;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Lambert<CloudType>::~Lambert()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::Lambert<CloudType>::evaluatePendular
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    scalar Vtot = 1.25e-10; // temporarily

    vector r_AB = (pA.position() - pB.position());

    scalar dAEff = pA.d();

    scalar dBEff = pB.d();

    scalar r_AB_mag = mag(r_AB);

    scalar normalOverlapMag = 0.5*(dAEff + dBEff) - r_AB_mag;

    scalar S = -normalOverlapMag;

    scalar Srup = (1+0.5*contactAngle_)*pow(Vtot, 1./3.);

    if (S < Srup)
    {
        //Pendular bridge formed

        vector rHat_AB = r_AB/(r_AB_mag + VSMALL);

        // Effective radius
        scalar R = 0.5*dAEff*dBEff/(dAEff + dBEff);

        // Normal force
        scalar capMag =
            4*mathematical::pi
            *R*surfaceTension_*cos(contactAngle_);

        if(S > 0)
        {
            capMag /= 1 + 1 /
                (sqrt(1+Vtot/(mathematical::pi*R*S*S))-1);
        }

        vector fN_AB = -capMag * rHat_AB;

        pA.f() += fN_AB;
        pB.f() += -fN_AB;
    }
}


// ************************************************************************* //
