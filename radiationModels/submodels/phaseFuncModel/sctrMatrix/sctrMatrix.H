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
    Foam::radiation::sctrMatrix

Description
    Struct to store scattering factors between a pair of directions.

SourceFiles
    sctrMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef sctrMatrix_H
#define sctrMatrix_H

#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class sctrMatrix Declaration
\*---------------------------------------------------------------------------*/

class sctrMatrix
{

        //- List of specular index in patch faces
        List<scalar> sctrFactors_;
		
		label directions_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sctrMatrix(const sctrMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const sctrMatrix&);


public:

    // Constructors

        //- Construct from components
        sctrMatrix
        (
            const label directions
        );
		


    //- Destructor
    ~sctrMatrix();


    // Member functions

		//- Return non-const access to sctrFactor from dir i to dir j
			scalar& Item(const label i, const label j);
			
		//- Return const access to sctrFactor from dir i to dir j
		scalar constItem(const label i, const label j) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
