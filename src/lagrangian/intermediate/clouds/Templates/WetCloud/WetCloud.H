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
    Foam::WetCloud

Description
    Adds coolisions to kinematic clouds

SourceFiles
    WetCloudI.H
    WetCloud.C

\*---------------------------------------------------------------------------*/

#ifndef WetCloud_H
#define WetCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CollisionModel;

/*---------------------------------------------------------------------------*\
                       Class WetCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WetCloud
:
    public CloudType
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef WetCloud<CloudType> wetCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<WetCloud<CloudType> > cloudCopyPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WetCloud(const WetCloud&);

        //- Disallow default bitwise assignment
        void operator=(const WetCloud&);


protected:

    // Protected data

        // References to the cloud sub-models

            //- Collision model
            autoPtr<CollisionModel<WetCloud<CloudType> > >
                collisionModel_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();

/*
        // Cloud evolution functions

            //- Move-collide particles
            template<class TrackData>
            void moveCollide(TrackData& td, const scalar deltaT);

            //- Reset state of cloud
            void cloudReset(WetCloud<CloudType>& c);
*/

public:

    // Constructors

        //- Construct given carrier gas fields
        WetCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );
/*
        //- Copy constructor with new name
        WetCloud
        (
            WetCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        WetCloud
        (
            const fvMesh& mesh,
            const word& name,
            const WetCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new WetCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new WetCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~WetCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const WetCloud& cloudCopy() const;

            //- If the collision model controls the wall interaction,
            //  then the wall impact distance should be zero.
            //  Otherwise, it should be allowed to be the value from
            //  the Parcel.
            virtual bool hasWallImpactDistance() const;


            // Sub-models

                //- Return const access to the collision model
                inline const CollisionModel<WetCloud<CloudType> >&
                    collision() const;

                //- Return reference to the collision model
                inline CollisionModel<WetCloud<CloudType> >&
                    collision();


        // Evolution

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackData>
            void motion(TrackData& td);
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "WetCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WetCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
