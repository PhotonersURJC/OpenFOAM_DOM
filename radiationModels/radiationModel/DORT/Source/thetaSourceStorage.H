/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Namespace
    Foam::radiation

Description
    Struct to storage source emmision power in all directions in a boundary.
	It stores a double value in each band and theta value of the quadrature.
	Every face in the boundary emits exactly the same intensities
Class
    Foam::radiation::thetaSourceStorage

Description
    

SourceFile
    thetaSourceStorage.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelThetaSourceStorage_H
#define radiationModelThetaSourceStorage_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOdictionary.H"

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class thetaSourceStorage Declaration
\*---------------------------------------------------------------------------*/

class thetaSourceStorage
{
protected:

    // Protected data

	//- patch Q in each Band and Theta Index
	List<List<scalar> > Qlist_;
	
	//quadrature theta divisions
	label nTheta_;
	
	//quadrature phi divisiones
	label nPhi_;
	

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thetaSourceStorage(const thetaSourceStorage&);

        //- Disallow default bitwise assignment
        void operator=(const thetaSourceStorage&);


public:

    // Constructors

        //- Construct from components
        thetaSourceStorage(const label nBands, const label nTheta, const label nPhi);


    //- Destructor
    ~thetaSourceStorage();


    // Member Functions

        // Edit

            //- Store Fluxes
            void setLambertQs(const List<scalar> Qlist, const scalar expFactor, const scalar correctionFactor);

        // Access

            //- Access to specific value
            scalar getValue(const label face, const label direction, const label band) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
