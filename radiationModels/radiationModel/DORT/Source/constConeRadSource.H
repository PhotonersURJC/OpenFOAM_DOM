/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::radiation::constConeRadSource

Description

    Cone-shaped emmision source. Quadrature rotates to fit its Z-axis in the
	cone central direction, and sets the nearest cut between ordinates in the
	cone opening perimeter. Theta increments in each side of the cut are
	corrected to be equally spaced in their respective regions.

SourceFiles
    constConeRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceconstConeRadSource_H
#define radSourceconstConeRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class constConeRadSource Declaration
\*---------------------------------------------------------------------------*/

class constConeRadSource
:
    public radSource
{

	PtrList<booleanSourceStorage> sourceValues_;
    List<vector> patchNormals_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        constConeRadSource(const constConeRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const constConeRadSource&);


public:

    //- Runtime type information
    TypeName("constConeRadSource");


    // Constructors

        //- Construct from components
        constConeRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        constConeRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~constConeRadSource();


    // Member functions

        // Edit

	    //- Return main axis to draw the new quadrature
	    const List<vector> getAxis(label nTheta, label nPhi);

            //- Store values of emission
            virtual void generate(const List<vector> directions);
			
            //- Read radiation properties dictionary
            bool read();

        // Access
		
			//- Acces to specific stored value
            virtual scalar calculate(const word patchID, const label face, const label direction, const label band) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
