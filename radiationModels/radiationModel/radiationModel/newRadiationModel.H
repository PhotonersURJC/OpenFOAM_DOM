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
    Namespace for radiation modelling

Class
    Foam::radiation::radiationModel

Description
    Top level model for radiation modelling

SourceFiles
    newRadiationModel.C
    newRadiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef newRadiationModel_H
#define newRadiationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "fvMatricesFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes

class basicThermo;
class fvMesh;

namespace radiation
{
// Forward declaration of classes
class phaseFuncModel;
class radExtintionModel;

/*---------------------------------------------------------------------------*\
                       Class noTradiationModel Declaration
\*---------------------------------------------------------------------------*/

class newRadiationModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Radiation model on/off flag
        Switch radiation_;
		
		//- Reference to temperature field
		const volScalarField& T_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

		//- Scattering phase function
		autoPtr<phaseFuncModel> phaseFunc_;

        //- Extintion model
        autoPtr<radExtintionModel> extintion_;


private:

    // Private Member Functions
	
        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        newRadiationModel(const newRadiationModel&);

        //- Disallow default bitwise assignment
        void operator=(const newRadiationModel&);


public:

    //- Runtime type information
    TypeName("newRadiationModel");


    // Declare runtime constructor selection table

		declareRunTimeSelectionTable
        (
            autoPtr,
            newRadiationModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            newRadiationModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T
            ),
            (dict, T)
        );


    // Constructors

        //- Null constructor
        newRadiationModel(const volScalarField& T);

        //- Construct from components
        newRadiationModel(const word& type, const volScalarField& T);
		
        //- Construct from components
        newRadiationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T
        );


    // Selectors

        //- Return a reference to the selected radiation model
		static autoPtr<newRadiationModel> New(const volScalarField& T);
        //- Return a reference to the selected radiation model
        static autoPtr<newRadiationModel> New
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~newRadiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual bool correct();

            //- Solve radiation equation(s)
            virtual bool calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Radiation model on/off flag
            const Switch radiation() const
            {
                return radiation_;
            }

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField> Ru() const = 0;

            //- Energy source term
            virtual tmp<fvScalarMatrix> Sh
            (
                const basicThermo& thermo,
                const volScalarField& he
            ) const;

            //- Temperature source term
            virtual tmp<fvScalarMatrix> ST
            (
                const dimensionedScalar& rhoCp,
                volScalarField& T
            ) const;


            //- Access to phase function model
            const phaseFuncModel& phaseFunc() const;

            //- Access to extintion model
            const radExtintionModel& extintion() const;
			
			//- Set values of scattering matrix
			void sctrMatrixGen(const List<vector> directions, const label model);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToRadiationRunTimeSelectionTables(model)                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        newRadiationModel,                                                    \
        model,                                                                \
        dictionary                                                            \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        newRadiationModel,                                                    \
        model,                                                                \
        T                                                                  \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
