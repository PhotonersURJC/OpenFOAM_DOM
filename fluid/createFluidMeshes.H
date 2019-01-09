    
	const wordList fluidNames(rp["fluid"]);

    PtrList<fvMesh> fluidRegions(fluidNames.size());
	PtrList<volScalarField> fluidTemperatures(fluidNames.size());

    forAll(fluidNames, i)
    {
        Info<< "Create fluid mesh for region " << fluidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    fluidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
    // Initialise fluid field pointer lists
    PtrList<radiation::newRadiationModel> fluidRadiation(fluidRegions.size());
	
	
    // Populate fluid field pointer lists
    forAll(fluidTemperatures, i)
    {
        Info<< "*** Reading fluid radiation properties for region "
            << fluidRegions[i].name() << nl << endl;
		// Create a null temperature field to inactivate thermal calculations
		fluidTemperatures.set(
		i,
		new volScalarField
		(
		   IOobject
		   (
				"T",
				runTime.timeName(),
				fluidRegions[i],
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
		   fluidRegions[i]
		)
		);
        fluidRadiation.set
        (
            i,
            radiation::newRadiationModel::New(fluidTemperatures[i])
        );
    }

	