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
    Struct to storage source emmsion power in all directions in a boundary.
	It stores a single double value in each face, and a boolean value
	in each face and direction, identifying directions where there is emmision

Class
    Foam::radiation::booleanSourceStorage

Description
    

SourceFile
    booleanSourceStorage.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelBooleanSourceStorage_H
#define radiationModelBooleanSourceStorage_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOdictionary.H"

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class booleanSourceStorage Declaration
\*---------------------------------------------------------------------------*/

class booleanSourceStorage
{
protected:

    // Protected data

	//- patch Q in each Band
	List<scalar> Qlist_;
	
	//- Source emmision in each face and direction
	List<List<bool> > criteria_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        booleanSourceStorage(const booleanSourceStorage&);

        //- Disallow default bitwise assignment
        void operator=(const booleanSourceStorage&);


public:

    // Constructors

        //- Construct from components
        booleanSourceStorage(const label nBands, const label nFaces, const label nDirections);


    //- Destructor
    ~booleanSourceStorage();


    // Member Functions

        // Edit

            //- Store Fluxes
            void setQs(const List<scalar> Qlist);


			//- Store boolean values
			void setBools(const label face, const label direction, const bool criteria);


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
