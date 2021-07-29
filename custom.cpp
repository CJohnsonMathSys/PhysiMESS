#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	Cell_Definition* pCD = find_cell_definition("worm"); 
	pCD->functions.update_phenotype = NULL;
	pCD->functions.custom_cell_rule = custom_worm_function;
	pCD->functions.contact_function = worm_contact_function; 
	pCD->phenotype.mechanics.attachment_elastic_constant = 0.03; 
	
	
	//cell_defaults.functions.update_phenotype = NULL; 
	//cell_defaults.functions.custom_cell_rule = custom_worm_function; 
	//ell_defaults.functions.contact_function = worm_contact_function; 
	
	//cell_defaults.phenotype.mechanics.attachment_elastic_constant = 0.03; 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

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
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );

			pC->custom_data["head"] = UniformRandom(); 
			pC->custom_data["head_start"] = pC->custom_data["head"];
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	// create some fibres 
	
	Fibre* pF;
	
	for( int n = 0 ; n < parameters.ints("number_of_fibres") ; n++ )
		{
			std::vector<double> fibre_start_position = {0,0,0}; 
			fibre_start_position[0] = Xmin + UniformRandom()*Xrange; 
			fibre_start_position[1] = Ymin + UniformRandom()*Yrange; 
			fibre_start_position[2] = Zmin + UniformRandom()*Zrange; 

			pF->assign_fibre_start_position( fibre_start_position );
                        pF->assign_fibre_radius( parameters.ints("fibre_radius");
                        pF->assign_fibre_length( parameters.ints("fibre_length");
		}
	}
	std::cout << " fibres have been initialised at random positions" << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_fibres_from_pugixml(); 
        std::cout << " fibres have been initialised from CSV file" << std::endl; 
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

// These functions are used by PhysiCell to choose a migration direction whenever 
// a cell makes a turn. 

void stretch_migration_direction( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.motility.chemotaxis_direction = -1; 
	
	phenotype.motility.migration_speed = 0.2; 
	phenotype.motility.migration_bias = 0.5; 
	phenotype.motility.persistence_time = 10; 

	return chemotaxis_function( pCell,phenotype,dt); 
}


// PhysiCell has a built-in contact function for elastic spring-like attachmetns 

void worm_contact_function( Cell* pMe, Phenotype& phenoMe, 
	Cell* pOther, Phenotype& phenoOther, double dt )
{ 
	standard_elastic_contact_function(pMe,phenoMe,pOther,phenoOther,dt);
	if( pMe->state.number_of_attached_cells() > 1 )
	{
		double head_me = pMe->custom_data["head"];
		double head_other = pOther->custom_data["head"]; 
		// make the transfer
		if( head_me > head_other )
		{
			double amount_to_transfer = dt * pMe->custom_data["transfer_rate"] 
			* (head_me - head_other ); 
			pMe->custom_data["head"] -= amount_to_transfer; 
			#pragma omp critical
			{ pOther->custom_data["head"] += amount_to_transfer; }
		}
	}

}

void custom_worm_function( Cell* pCell, Phenotype& phenotype , double dt )
{ 
	// bookkeeping 
	static int nSignal = microenvironment.find_density_index("signal");
	// look for cells to form attachments, if 0 attachments 
	int number_of_attachments = pCell->state.number_of_attached_cells();
	std::vector<Cell*> nearby = pCell->nearby_interacting_cells(); 
	if( number_of_attachments == 0 )
	{ 
		int n = 0; 
		while( number_of_attachments < (int) pCell->custom_data["max_attachments"] && n < nearby.size() )
		{
			if( nearby[n]->state.number_of_attached_cells() < nearby[n]->custom_data["max_attachments"] )
			{ 
				attach_cells( nearby[n] , pCell ); 
				number_of_attachments++;
			}
			n++; 
		}
	}
	// if no attachments, use chemotaxis
	if( number_of_attachments == 0 ){
		pCell->functions.update_migration_bias = chemotaxis_function; 
	}

	// if 1 attachment, do some logic  
	if( number_of_attachments == 1 )
	{ 
		// constant expression in end cells
		pCell->custom_data["head"] = pCell->custom_data["head_start"];
		// am I the head? 
		bool head = false; 
		if( pCell->custom_data["head"] > pCell->state.attached_cells[0]->custom_data["head"] )
		{ 
			head = true; 
		} 
		if( head )
		{ 
			pCell->functions.update_migration_bias = head_migration_direction; 
		}else{
				pCell->functions.update_migration_bias = tail_migration_direction; 
		}
		phenotype.secretion.secretion_rates[nSignal] = 100; 
	} 

	// if 2 or more attachments, use middle 
	if( number_of_attachments > 1 )
	{ 
		pCell->functions.update_migration_bias = middle_migration_direction;
		phenotype.secretion.secretion_rates[nSignal] = 1; 
	} 

	return; 
}

void head_migration_direction( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.motility.chemotaxis_direction = -1; 
	phenotype.motility.migration_speed = 0.75; 
	phenotype.motility.migration_bias = 0.5; 
	phenotype.motility.persistence_time = 60;

	return chemotaxis_function( pCell,phenotype,dt); 
}

void tail_migration_direction( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.motility.chemotaxis_direction = -1; 
	phenotype.motility.migration_speed = 0; 
	phenotype.motility.migration_bias = 0.5; 
	phenotype.motility.persistence_time = 100; 
	return chemotaxis_function( pCell,phenotype,dt); 
}

void middle_migration_direction( Cell* pCell, Phenotype& phenotype , double dt )
{ // get velocity from "Upstream" 
	Cell* pUpstream = pCell->state.attached_cells[0];

	if( pCell->state.attached_cells[1]->custom_data["head"] > pCell->state.attached_cells[0]->custom_data["head"] )
	{
		pUpstream = pCell->state.attached_cells[1]; 
	}
}

// only head cell is red
std::vector<std::string> worm_coloring_function( Cell* pCell )
{
	if( pCell->state.number_of_attached_cells() == 0 )
	{ 
		return { "grey", "black", "grey", "grey"}; 
	}

	if( pCell->state.number_of_attached_cells() == 1 && 
	pCell->custom_data["head"] > pCell->state.attached_cells[0]->custom_data["head"] )
	{ 
		return { "red", "black", "red", "red"}; 
	}

	if( pCell->state.number_of_attached_cells() == 2 )
	{ 
	return { "blue", "black", "blue", "blue"}; 
	}

	return { "yellow", "black", "yellow", "yellow" }; 
}
