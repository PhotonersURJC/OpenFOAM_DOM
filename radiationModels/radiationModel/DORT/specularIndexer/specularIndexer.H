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
    Foam::radiation::specularIndexer

Description
    Storage of specular reflection directions for each discrete Ordinate in
	a given face.

SourceFiles
    specularIndexer.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSpecularIndexer_H
#define radiationSpecularIndexer_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class specularIndexer Declaration
\*---------------------------------------------------------------------------*/

class specularIndexer
{
		//- Determines where to look for spotters, if any (else -1)
		List<short> outers_;
		
		//- Reflected incoming DOs looking at current outgoing DO
		List<List<short > > spotters_;
		
		//- Factors to add intensities of spotters
		List<List<scalar> > viewFactors_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularIndexer(const specularIndexer&);

        //- Disallow default bitwise assignment
        void operator=(const specularIndexer&);


public:

    // Constructors
		
        //- Construct from components
        specularIndexer
        (
			const List<List< label > > facesNeighbours,
			const List<vector> mainDirs,
			const List<scalar> faceMagnitudes,
			const vector faceN,
			const List<vector> quadVertex,
			const List<face> quadFaces,
			const bool secondOrder
        );

    //- Destructor
    ~specularIndexer();


    // Member functions

	    //- Return the reflected Field for a given face
            scalarField getSpec(const scalarField originalI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
