/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <math.h>  
#include <chrono>
#include <random>


void create_cell_types( void ){
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
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void ){
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void ){
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true ){
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	

	// N.B following fibre assignment only works if fibres are considered as another cell agent (uses all_cells)
	bool isFibreFromFile = false; // enable the manual input

    double cell_velocity_max = 0.16666;
    double vel_adhesion = 0.0003; //0.3; THESE VALUES NEED SETTING
    double vel_contact = 0.01; //0.001; THESE VALUES NEED SETTING

	unsigned fail_count=0;

	for( int i=0; i < (*all_cells).size(); i++ )
	{
		(*all_cells)[i]->parameters.mCellVelocityMaximum= cell_velocity_max;

		if( (*all_cells)[i]->type_name == "fibre" )
		{
			isFibreFromFile = true;

			// concerned with a fibre, position provided in CSV, test whether out of bounds
	
			//set fibre length as normally distributed around 75
			//double fibreLength = NormalRandom(75.,5.);
            double fibreLength = 40; //value used for testing

			// set parameters
			(*all_cells)[i]->parameters.mLength = fibreLength/2.0;
			(*all_cells)[i]->parameters.mRadius = 2.0;
			//std::cout << " fibre length is " << fibreLength << std::endl;

			(*all_cells)[i]->parameters.mVelocityAdhesion = vel_adhesion;
			(*all_cells)[i]->parameters.mVelocityContact = vel_contact;

			//assign fibre orientation as a random vector from points on unit sphere/circle
			(*all_cells)[i]->assign_orientation();
            if( default_microenvironment_options.simulate_2D == true ) {
                (*all_cells)[i]->state.orientation = UniformOnUnitCircle();
                (*all_cells)[i]->state.orientation[2] = 0.0;
                //std::cout << " fibre orientation in 2D is " << (*all_cells)[i]->state.orientation[0] << " " << (*all_cells)[i]->state.orientation[1] << std::endl;
            }
            else{
                (*all_cells)[i]->state.orientation = UniformOnUnitSphere();
                //std::cout << " fibre orientation in 3D is " << (*all_cells)[i]->state.orientation[0] << " " << (*all_cells)[i]->state.orientation[1] << " " << (*all_cells)[i]->state.orientation[2] << std::endl;
            }

			// start and end points of a fibre are calculated from fibre center
			double xs = (*all_cells)[i]->position[0] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
			double xe = (*all_cells)[i]->position[0] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
			double ys = (*all_cells)[i]->position[1] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
			double ye = (*all_cells)[i]->position[1] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
            double zs = 0.0;
            double ze = 0.0;
            if( default_microenvironment_options.simulate_2D == true ) {
                //std::cout << " fibre endpoints in 2D are " << xs << " " << ys << " and " << xe << " " << ye << std::endl;
            }
            else if( default_microenvironment_options.simulate_2D == false ) {
                zs = (*all_cells)[i]->position[2] -
                            (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                ze = (*all_cells)[i]->position[2] +
                            (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                //std::cout << " fibre endpoints in 3D are " << xs << " " << ys << " " << zs << " and " << xe << " " << ye << " " << ze << std::endl;
            }

			// check whether a fibre end point leaves the domain and if so initialise fibre again
			// assume user placed the centre of fibre within the domain so reinitialise orientation,
            // break after 10 failures
			while ((xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax) && fail_count<10) {
				//std::cout << "fibre position is " << xs-Xmin << " " << xe-Xmin << " " << std::endl;
				std::cout << "!!! The fibre has a portion outside of the domain - trying again !!!" << std::endl;
				(*all_cells)[i]->state.orientation[0] = UniformOnUnitCircle()[0];
				xs = (*all_cells)[i]->position[0] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
				xe = (*all_cells)[i]->position[0] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
				//std::cout << "new fibre position is " << xs-Xmin << " " << xe-Xmin << " " << std::endl;
				fail_count++;
			}
		
			while ((ys < Ymin || ye > Ymax || ye < Xmin || ys > Xmax) && fail_count<10) {
				//std::cout << "fibre position is " << ys-Ymin << " " << ye-Ymin << " " << std::endl;
				std::cout << "!!! The fibre has a portion outside of the domain - trying again !!!" << std::endl;
				(*all_cells)[i]->state.orientation[1] = UniformOnUnitCircle()[1];
				ys = (*all_cells)[i]->position[1] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
				ye = (*all_cells)[i]->position[1] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
				//std::cout << "new fibre position is " << ys-Ymin << " " << ye-Ymin << " " << std::endl;
				fail_count++;
			}

            // the following needs re-writing properly to handle the 3D case
            if( default_microenvironment_options.simulate_2D == false ) {
			    while (zs < Zmin || ze > Zmax || ze < Xmin || zs > Xmax) {
                    std::cout << "!!! The fibre has a portion outside of the domain - trying again !!!" << std::endl;
                }
                (*all_cells)[i]->state.orientation[2] = UniformOnUnitSphere()[2];
			    zs = (*all_cells)[i]->position[2] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[2];
			    ze = (*all_cells)[i]->position[2] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[2];
			}

			if(fail_count>=10){
				// failed to place the fibre within the domain bounds, delete fibre
				delete_cell((*all_cells)[i]); 
			}
			fail_count=0;

		}

		//###########################################//
		//   this bit a hack for PacMan and maze	 //
		//###########################################//
		else if((*all_cells)[i]->type_name == "fibre_vertical" )
		{
			isFibreFromFile = true;

			// concerned with a fibre, position provided in CSV, test whether out of bounds
			bool isFibreOutOfBounds = false;

			//set fibre length as normally distributed around 75
			//double fibreLength = NormalRandom(75.,5.);
			double fibreLength = 40; //value used for testing

			// set parameters
			(*all_cells)[i]->parameters.mLength = fibreLength/2.0;
            (*all_cells)[i]->parameters.mRadius = 2.0;
			//std::cout << " fibre length is " << fibreLength << std::endl;

			(*all_cells)[i]->parameters.mVelocityAdhesion = vel_adhesion;
			(*all_cells)[i]->parameters.mVelocityContact = vel_contact;

			//assign fibre orientation - vertical i.e. aligned with y in xy plane
			(*all_cells)[i]->assign_orientation();
			(*all_cells)[i]->state.orientation[0] = 0.0;
			(*all_cells)[i]->state.orientation[1] = 1.0;
			(*all_cells)[i]->state.orientation[2] = 0.0;

			// start and end points of a fibre are calculated from fibre center
			double xs = (*all_cells)[i]->position[0] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
			double xe = (*all_cells)[i]->position[0] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
			double ys = (*all_cells)[i]->position[1] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
			double ye = (*all_cells)[i]->position[1] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
            double zs = 0.0;
            double ze = 0.0;
            if( default_microenvironment_options.simulate_2D == false ) {
                zs = (*all_cells)[i]->position[2] -
                            (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                ze = (*all_cells)[i]->position[2] +
                            (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
            }

			// check whether a fibre end point leaves the domain and if so initialise fibre again
			// assume user placed the centre of fibre within the domain so force delete of this fibre
			while ( (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax) && fail_count<10) {
				std::cout << "!!! VERTICAL fibre has a portion outside of the domain - deleting fibre !!!" << std::endl;
				fail_count=11;
			}
	
			while ((ys < Ymin || ye > Ymax || ye < Xmin || ys > Xmax) && fail_count<10) {
				std::cout << "!!! VERTICAL fibre has a portion outside of the domain - deleting fibre !!!" << std::endl;
				fail_count=11;
			}

            if( default_microenvironment_options.simulate_2D == false ) {
                while (zs < Zmin || ze > Zmax || ze < Xmin || zs > Xmax) {
                    std::cout << "!!! VERTICAL fibre has a portion outside of the domain - deleting fibre !!!" << std::endl;
                    fail_count=11;
                }
            }

            if(fail_count>=10){
                // failed to place the fibre within the domain bounds, delete fibre
                delete_cell((*all_cells)[i]);
            }
            fail_count=0;

            // relabel so that the rest of the code works (HACK)
            (*all_cells)[i]->type_name = "fibre";

        }

		else if((*all_cells)[i]->type_name == "fibre_horizontal" )
		{
			isFibreFromFile = true;

			// concerned withe a fibre, position provided in CSV, test whether out of bounds
			bool isFibreOutOfBounds = false;

			//set fibre length as normally distributed around 75
			//double fibreLength = NormalRandom(75.,5.);
			double fibreLength = 40; //value used for testing

			// set parameters
			(*all_cells)[i]->parameters.mLength = fibreLength/2.0;
            (*all_cells)[i]->parameters.mRadius = 2.0;
			//std::cout << " fibre length is " << fibreLength << std::endl;

			(*all_cells)[i]->parameters.mVelocityAdhesion = vel_adhesion;
			(*all_cells)[i]->parameters.mVelocityContact = vel_contact;

            //assign fibre orientation - horizontal i.e. aligned with x in xy plane
			(*all_cells)[i]->assign_orientation();
			(*all_cells)[i]->state.orientation[0] = 1.0;
			(*all_cells)[i]->state.orientation[1] = 0.0;
			(*all_cells)[i]->state.orientation[2] = 0.0;

			// start and end points of a fibre are calculated from fibre center
			double xs = (*all_cells)[i]->position[0] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
			double xe = (*all_cells)[i]->position[0] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[0];
			double ys = (*all_cells)[i]->position[1] - (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
			double ye = (*all_cells)[i]->position[1] + (*all_cells)[i]->parameters.mLength*(*all_cells)[i]->state.orientation[1];
            double zs = 0.0;
            double ze = 0.0;
            if( default_microenvironment_options.simulate_2D == false ) {
                zs = (*all_cells)[i]->position[2] -
                     (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
                ze = (*all_cells)[i]->position[2] +
                     (*all_cells)[i]->parameters.mLength * (*all_cells)[i]->state.orientation[2];
            }

			// check whether a fibre end point leaves the domain and if so initialise fibre again
			// assume user placed the centre of the fibre within the domain so force delete of this fibre
			while ((xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax) && fail_count<10) {
				std::cout << "!!! HORIZONTAL fibre has a portion outside of the domain - trying again !!!" << std::endl;
				fail_count=11;
			}
			fail_count=0;
			while ((ys < Ymin || ye > Ymax || ye < Xmin || ys > Xmax) && fail_count<10) {
				std::cout << "!!! HORIZONTAL fibre has a portion outside of the domain - trying again !!!" << std::endl;
				fail_count=11;
			}

            if( default_microenvironment_options.simulate_2D == false ) {
                while (zs < Zmin || ze > Zmax || ze < Xmin || zs > Xmax) {
                    std::cout << "!!! HORIZONTAL fibre has a portion outside of the domain - trying again !!!"
                              << std::endl;
                    fail_count = 11;
                }
            }

            if(fail_count>=10){
                // failed to place the fibre within the domain bounds, delete fibre
                delete_cell((*all_cells)[i]);
            }
            fail_count=0;

            // relabel so that the rest of the code works (HACK)
            (*all_cells)[i]->type_name = "fibre";

        }

		else
		{
			// type is a normal cell
		}
	}

	bool isAddFibres = true; // disable the manual input in favour of csv
	if(isAddFibres && !isFibreFromFile)
	{
		// fibres have not been added from the file but do want fibres
		// create some of each type of cell 
		Cell* pC;
		
		std::vector<double> position = {0, 0, 0};

		for( int k=0; k < cell_definitions_by_index.size() ; k++ ){
			Cell_Definition* pCD = cell_definitions_by_index[k]; 
			std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
			for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ ) {

                position[0] = Xmin + UniformRandom() * Xrange;
                position[1] = Ymin + UniformRandom() * Yrange;
                position[2] = Zmin + UniformRandom() * Zrange;

                pC = create_cell(*pCD);

                pC->parameters.mCellVelocityMaximum = cell_velocity_max;

                pC->assign_position(position);
            }

            if(pCD->name == "fibre"){
                for ( int nf = 0 ; nf < parameters.ints("number_of_fibres") ; nf++ ) {

                    position[0] = Xmin + UniformRandom() * Xrange;
                    position[1] = Ymin + UniformRandom() * Yrange;
                    position[2] = Zmin + UniformRandom() * Zrange;

                    pC = create_cell(*pCD);

                    pC->parameters.mCellVelocityMaximum = cell_velocity_max;

                    //set fibre length as normally distributed around 75
                    //double fibreLength = NormalRandom(75.,5.);
                    double fibreLength = 40; //value used for testing

                    // set parameters
                    pC->parameters.mLength = fibreLength/2.0;
                    pC->parameters.mRadius = 2.0;

                    pC->parameters.mVelocityAdhesion = vel_adhesion;
                    pC->parameters.mVelocityContact = vel_contact;

                    //assign fibre orientation as a random vector from points on unit sphere/circle.
                    pC->assign_orientation();
                    if( default_microenvironment_options.simulate_2D == true ) {
                        pC->state.orientation = UniformOnUnitCircle();
                        pC->state.orientation[2] = 0.0;
                    }
                    else{
                        pC->state.orientation = UniformOnUnitSphere();
                    }

                    // start and end points of a fibre are calculated from fibre center
                    double xs = position[0] - pC->parameters.mLength*pC->state.orientation[0];
                    double xe = position[0] + pC->parameters.mLength*pC->state.orientation[0];
                    double ys = position[1] - pC->parameters.mLength*pC->state.orientation[1];
                    double ye = position[1] + pC->parameters.mLength*pC->state.orientation[1];
                    double zs = 0.0;
                    double ze = 0.0;
                    if( default_microenvironment_options.simulate_2D == false ) {
                        zs = position[2] - pC->parameters.mLength * pC->state.orientation[2];
                        ze = position[2] + pC->parameters.mLength * pC->state.orientation[2];
                    }

                    // check whether a fibre end point leaves the domain and if so initialise fibre again
                    while (xs < Xmin || xe > Xmax || xe < Xmin || xs > Xmax) {
                        std::cout << "!!! The fibre has a portion outside of the domain - trying again !!!" << std::endl;
                        position[0] = Xmin + UniformRandom() * Xrange;
                        xs = position[0] - pC->parameters.mLength*pC->state.orientation[0];
                        xe = position[0] + pC->parameters.mLength*pC->state.orientation[0];
                    }

                    while (ys < Ymin || ye > Ymax || ye < Xmin || ys > Xmax) {
                        std::cout << "!!! The fibre has a portion outside of the domain - trying again !!!" << std::endl;
                        position[1] = Ymin + UniformRandom() * Yrange;
                        ys = position[1] - pC->parameters.mLength*pC->state.orientation[1];
                        ye = position[1] + pC->parameters.mLength*pC->state.orientation[1];
                    }

                    if( default_microenvironment_options.simulate_2D == false ) {
                        while (zs < Zmin || ze > Zmax || ze < Xmin || zs > Xmax) {
                            std::cout << "!!! The fibre has a portion outside of the domain - trying again !!!"
                                      << std::endl;
                            position[2] = Zmin + UniformRandom() * Zrange;
                            zs = position[2] - pC->parameters.mLength * pC->state.orientation[2];
                            ze = position[2] + pC->parameters.mLength * pC->state.orientation[2];
                        }
                    }

                    pC->assign_position(position);

                }

			}
		}

	}

	std::cout << std::endl; 
	
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
