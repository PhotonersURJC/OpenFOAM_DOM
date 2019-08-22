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
    Foam::radiation::parallelRadSource

Description

    Parallel emmision source. Light is emmited in a single discrete ordinate.
	Quadrature is rotated to place a discrete ordinate near the quadrature
	equator in the light direction.
	
	Discrete ordinates near the equator are the biggest, and have a squared
	shape, been considered as "better" ordinates.

SourceFiles
    parallelRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceparallelRadSource_H
#define radSourceparallelRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class parallelRadSource Declaration
\*---------------------------------------------------------------------------*/

class parallelRadSource
:
    public radSource
{
	// Private data
	PtrList<booleanSourceStorage> sourceValues_;
    List<vector> patchNormals_;
    label quadIndex_;
    scalar omega_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        parallelRadSource(const parallelRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const parallelRadSource&);


public:

    //- Runtime type information
    TypeName("parallelRadSource");


    // Constructors

        //- Construct from components
        parallelRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        parallelRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~parallelRadSource();


    // Member functions
	
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
