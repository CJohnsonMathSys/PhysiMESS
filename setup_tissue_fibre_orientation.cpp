void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;

	double meanAngle = 90;
	double variationAngle = 5;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; // position of the cell or centre of the fibre
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );

			double angle = 0.0;

			if(pC->type_name == "fibre")
			{
				double w, y1;
				double x1, x2;
				do {
					x1 = 2.0 * UniformRandom() - 1.0;
					x2 = 2.0 * UniformRandom() - 1.0;
					w = x1 * x1 + x2 * x2;
				} while ( w >= 1.0 );
					
				w = sqrt( (-2.0 * log( w ) ) / w );

				// orientation angle of the fibre (double for 2D)
				angle = meanAngle + variationAngle * x1 * w;

			}

			// if the cell type was a fibre then the angle will be non-zero
			pC->custom_data.add_variable("orientation", angle );


		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}