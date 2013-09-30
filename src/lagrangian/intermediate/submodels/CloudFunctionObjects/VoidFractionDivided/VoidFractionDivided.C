/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "VoidFractionDivided.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFractionDivided<CloudType>::write()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFractionDivided<CloudType>::VoidFractionDivided
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudFunctionObject<CloudType>(dict, owner, typeName),
    thetaPtr_(NULL),
    refinementLevel_
    (
        readLabel
        (
            this->coeffDict().lookup("refinementLevel")
        )
    ),
    maxDistance_
    (
        readScalar
        (
            this->coeffDict().lookup("maxDistance")
        )
    )
{
    // Compute the number of the satellite particles

    DynamicList<vector> tmpSatPositions;

    // First satellite point: centre of particle
    tmpSatPositions.append(vector::zero);

    // Distribute other satellite points based on the refinement level
    scalar pi = Foam::constant::mathematical::pi;

    scalar r = scalar(0);
    scalar delR = scalar(1) / refinementLevel_;
    scalar theta = 0.25 * pi;
    scalar delTheta = 0.5 * pi;

    for(label rI = 0; rI < refinementLevel_; rI++)
    {
        r += delR;
        for(label thetaI = 0; thetaI < 4; thetaI++)
        {
            theta += delTheta * thetaI;

            // Original satellite point located in y-z plane
            vector spOrg = vector
            (
                scalar(0),
                r*cos(theta),
                r*sin(theta)
            );

            tmpSatPositions.append(spOrg);

            vector unitN;
            vector sp;

            // Rotate spOrg 90 degrees around y-axis
            unitN = vector(0, 1, 0);

            sp = unitN*(unitN&spOrg) + (spOrg-unitN*(unitN&spOrg))*cos(pi/2) - (spOrg^unitN)*sin(pi/2);

            tmpSatPositions.append(sp);

            // Rotate spPrg 90 degrees around z-axis
            unitN = vector(0, 0, 1);

            sp = unitN*(unitN&spOrg) + (spOrg-unitN*(unitN&spOrg))*cos(pi/2) - (spOrg^unitN)*sin(pi/2);

            tmpSatPositions.append(sp);
        }
    }

    satParNo_ = tmpSatPositions.size();
    satPositions_ = tmpSatPositions.xfer();

    buildInteractionList();
}


template<class CloudType>
Foam::VoidFractionDivided<CloudType>::VoidFractionDivided
(
    const VoidFractionDivided<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    thetaPtr_(NULL),
    refinementLevel_(vf.refinementLevel()),
    satParNo_(vf.satParNo()),
    satPositions_(vf.satPositions()),
    maxDistance_(vf.maxDistance())
{
	buildInteractionList();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFractionDivided<CloudType>::~VoidFractionDivided()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFractionDivided<CloudType>::preEvolve()
{
    if (thetaPtr_.valid())
    {
        thetaPtr_->internalField() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        thetaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Theta",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
}


template<class CloudType>
void Foam::VoidFractionDivided<CloudType>::postEvolve()
{
    volScalarField& theta = thetaPtr_();

    const fvMesh& mesh = this->owner().mesh();

    theta.internalField() /= mesh.time().deltaTValue()*mesh.V();

    CloudFunctionObject<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::VoidFractionDivided<CloudType>::postMove
(
    const parcelType& p,
    const label cellI,
    const scalar dt,
    const point& position0,
    bool&
)
{
/*
    volScalarField& theta = thetaPtr_();

    theta[cellI] += dt*p.nParticle()*p.volume();
*/

//    const fvMesh& mesh = this->owner().mesh();
    volScalarField& theta = thetaPtr_();

    List<vector> satParPositions
    (
        0.5 * p.d() * satPositions_
    );

    forAll(satParPositions, i)
    {
        satParPositions[i] += p.position();
    }

    scalar satParVol = p.volume()/satParNo_;

    label satCellI;
    forAll(satParPositions, i)
    {
        // satCellI = mesh.findCell(satParPositions[i]);

        satCellI = -1;

        forAll(allInteractionList_[cellI], listI)
        {
            label checkingCellI = allInteractionList_[cellI][listI];
            if(cellBbs_[checkingCellI].contains(satParPositions[i]))
            {
                // For some reason, pointInCell massively slows down simulation
                // isInside is essentially same as pointInCell, but it is faster
                //if(mesh.pointInCell(satParPositions[i], checkingCellI))
                if(isInside(satParPositions[i], checkingCellI))
                {
                    satCellI = checkingCellI;
                    break;
                }
            }
        }

        if(satCellI >= 0)
        {
            theta[satCellI] += dt*p.nParticle()*satParVol;
        }
    }


    labelHashSet lhash;
/*
    recursiveCellSearchForSatelliteParticle
    (
        p,
        cellI,
        dt,
        lhash
    );
*/
}

template<class CloudType>
void Foam::VoidFractionDivided<CloudType>::buildInteractionList()
{
    const fvMesh& mesh = this->owner().mesh();

    procBb_ = treeBoundBox(mesh.points());

    const vector interactionVec = maxDistance_*vector::one;

    extendedProcBb_ = treeBoundBox
    (
        procBb_.min() - interactionVec,
        procBb_.max() + interactionVec
    );

    allExtendedProcBbs_.setSize(Pstream::nProcs());

    allExtendedProcBbs_[Pstream::myProcNo()] = extendedProcBb_;

    Pstream::gatherList(allExtendedProcBbs_);

    Pstream::scatterList(allExtendedProcBbs_);

    cellBbs_.setSize(mesh.nCells());

    forAll(cellBbs_, cellI)
    {
        cellBbs_[cellI] = treeBoundBox
        (
            mesh.cells()[cellI].points
            (
                mesh.faces(),
                mesh.points()
            )
        );
    }

    extendedCellBbs_.setSize(mesh.nCells());

    forAll(extendedCellBbs_, cellI)
    {
        extendedCellBbs_[cellI] = treeBoundBox
        (
            cellBbs_[cellI].min() - interactionVec,
            cellBbs_[cellI].max() + interactionVec
        );
    }

    allInteractionList_.setSize(mesh.nCells());

    DynamicList<label> tmpInteractionList;
    forAll(allInteractionList_, cellI)
    {
        forAll(cellBbs_, otherCellI)
        {
            if
            (
                cellBbs_[cellI].overlaps(extendedCellBbs_[otherCellI])
            )
            {
                tmpInteractionList.append(otherCellI);
            }
        }
        allInteractionList_[cellI] = tmpInteractionList.xfer();
    }
}

template<class CloudType>
bool Foam::VoidFractionDivided<CloudType>::isInside(const point& p, label celli) const
{
    const fvMesh& mesh = this->owner().mesh();
    const labelList& f = mesh.cells()[celli];
    const labelList& owner = mesh.faceOwner();
    const vectorField& cf = mesh.faceCentres();
    const vectorField& Sf = mesh.faceAreas();

    bool inCell = true;

    forAll(f, facei)
    {
        label nFace = f[facei];
        vector proj = p - cf[nFace];
        vector normal = Sf[nFace];
        if (owner[nFace] != celli)
        {
            normal = -normal;
        }
        inCell = inCell && ((normal & proj) <= 0);
    }

    return inCell;
}

// ************************************************************************* //
