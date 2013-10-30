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
    Foam::Pendular

Description

SourceFiles
    Pendular.C

\*---------------------------------------------------------------------------*/

#ifndef Pendular_H
#define Pendular_H

#include "WetModel.H"
#include "InteractionLists.H"
//#include "WallSiteData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class PendularModel;

//template<class CloudType>
//class WallModel;


/*---------------------------------------------------------------------------*\
                     Class Pendular Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Pendular
:
    public WetModel<CloudType>
{
    // Private data

        //- PendularModel to calculate the interaction between two parcels
        autoPtr<PendularModel<CloudType> > pendularModel_;

//        //- WallModel to calculate the interaction between the parcel and walls
//        autoPtr<WallModel<CloudType> > wallModel_;

        //- Interactions lists determining which cells are in
        //  interaction range of each other
        InteractionLists<typename CloudType::parcelType> il_;


    // Private member functions

        //- Pre collision tasks
        void preInteraction();

        //- Interactions between parcels
        void parcelInteraction();

        //- Interactions between real (on-processor) particles
        void realRealInteraction();

        //- Interactions between real and referred (off processor) particles
        void realReferredInteraction();

//        //- Interactions with walls
//        void wallInteraction();

        //- Post collision tasks
        void postInteraction();

        //- Calculate the pair force between parcels
        void evaluatePendular
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;

//        //- Calculate the wall forces on a parcel
//        void evaluateWall
//        (
//            typename CloudType::parcelType& p,
//            const List<point>& flatSitePoints,
//            const List<WallSiteData<vector> >& flatSiteData,
//            const List<point>& sharpSitePoints,
//            const List<WallSiteData<vector> >& sharpSiteData
//        ) const;


public:

    //- Runtime type information
    TypeName("pendular");


    // Constructors

        //- Construct from components
        Pendular(const dictionary& dict, CloudType& owner);

        //- Construct copy
        Pendular(Pendular<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<WetModel<CloudType> > clone()
        {
            return autoPtr<WetModel<CloudType> >
            (
                new Pendular<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Pendular();


    // Member Functions

        //- Indicates whether model determines wall collisions or not,
        //  used to determine what value to use for wallImpactDistance
        virtual bool controlsWallInteraction() const;

        // Collision function
        virtual void bond();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Pendular.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
