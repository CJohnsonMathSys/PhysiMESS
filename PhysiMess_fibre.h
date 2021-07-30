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
# THIS NEW AGENT CLASS WAS CREATED BY                                         #
#                                                                             #
# Zoe Bell; Temitope Benson; Connah Johnson; Cicely Macnamara; James O'Neill; #
# Robyn Shuttleworth; Niki Tavakoli                                           #
#                                                                             #
# as part of WS2021 (Cluster 5)                                               #
###############################################################################
*/

#ifndef __PhysiCell_fibre_h__
#define __PhysiCell_fibre_h__

#include "./PhysiCell_custom.h" 

#include "../BioFVM/BioFVM.h"
#include "./PhysiCell_phenotype_fibre.h"
#include "./PhysiCell_fibre_container.h" // do we actually need a different container?
#include "./PhysiCell_constants.h"

#include "../modules/PhysiCell_settings.h" 

#include "./PhysiCell_standard_models.h"

using namespace BioFVM; 

namespace PhysiMess{
class Fibre_Container;

class Fibre_Parameters
{
 private:
 public:
	// store any fibre parameters we may need (not sure if length and radius go here yet)
	//double fibre_radius; // fibre radius
	//double fibre_length; // fibre length 
	
	Fibre_Parameters(); 
}; 

class Fibre_Definition
{
 private:
 public: 
	int type; 
	std::string name; 
 
	Microenvironment* pMicroenvironment; 
	
	Fibre_Parameters parameters; 
	Custom_Fibre_Data custom_data; 
	Fibre_Functions functions; 
	Fibre_Phenotype phenotype; // not sure if we want this as Phenotype_Fibre instead as per Robyns files

	Fibre_Definition();  // done 
	Fibre_Definition( Fibre_Definition& cd ); // copy constructor 
	Fibre_Definition& operator=( const Fibre_Definition& cd ); // copy assignment 
};

extern Fibre_Definition fibre_defaults; 

class Fibre_State
{
 private:
 public:
	std::vector<Fibre*> attached_fibre; 

	std::vector<Fibre*> neighbors; // not currently tracked! 
	std::vector<double> orientation;
	
	double simple_pressure; 
	
	int number_of_attached_fibres( void ); 
	
	Fibre_State(); 
};

class Fibre : public Basic_Agent 
{
 private: 
	Fibre_Container * container;
	int current_mechanics_voxel_index;
	int updated_current_mechanics_voxel_index; // keeps the updated voxel index for later adjusting of current voxel index
		
 public:
	std::string type_name; 
 
	Custom_Fibre_Data custom_data;
	Fibre_Parameters parameters;
	Fibre_Functions functions; 

	Fibre_State state; 
	Fibre_Phenotype phenotype; // not sure if we want this as Phenotype_Fibre instead
	
	void update_motility_vector( double dt_ );
	void advance_bundled_phenotype_functions( double dt_ ); 
	
	void add_potentials(Cell*);       // Add repulsive and adhesive forces.
	void set_previous_velocity(double xV, double yV, double zV);
	int get_current_mechanics_voxel_index();
	void turn_off_reactions(double); 		  // Turn off all the reactions of the cell
	
	bool is_out_of_domain;
	bool is_movable;
	
	void flag_for_division( void ); // done 
	void flag_for_removal( void ); // done 
	
	void start_death( int death_model_index ); 
	void lyse_fibre( void ); 

	Fibre* divide( void );
	void die( void ); 
	void step(double dt);
	Fibre();
	
	~Fibre(); 
	
	bool assign_position(std::vector<double> new_position);
	bool assign_position(double, double, double);
	void set_total_volume(double);
	
	double& get_total_volume(void); // NEW
	
	void set_target_volume(double); 
	void set_target_radius(double); 
	void set_radius(double); 
	
	// mechanics 
	void update_position( double dt ); //
	std::vector<double> displacement; // this should be moved to state, or made private  

	
	void assign_orientation();  // if set_orientaion is defined, uses it to assign the orientation
								// otherwise, it assigns a random orientation to the cell.
	
	void copy_function_pointers(Fibre*);
	
	void update_voxel_in_container(void);
	void copy_data(Fibre *);
	
	void ingest_fibre( Fibre* pFibre_to_eat ); // for use in predation, e.g., immune cells 

	void attach_fibre( Fibre* pAddMe ); // done 
	void detach_fibre( Fibre* pRemoveMe ); // done 
	void remove_all_attached_fibres( void ); // done 

	// I want to eventually deprecate this, by ensuring that 
	// critical BioFVM and PhysiCell data elements are synced when they are needed 
	
	void set_phenotype( Fibre_Phenotype& phenotype ); // no longer needed?
	void update_radius();
	Fibre_Container * get_container();
	
	std::vector<Fibre*>& fibres_in_my_container( void ); 
	std::vector<Fibre*> nearby_fibres( void ); 
	std::vector<Fibre*> nearby_interacting_fibres( void );  
	
	void convert_to_fibre_definition( Fibre_Definition& cd ); 
};

Cell* create_fibre( void );  
Cell* create_fibre( Fibre_Definition& cd );  


void delete_fibre( int ); 
void delete_fibre( Cell* ); 
void save_all_fibres_to_matlab( std::string filename ); 

//function to check if a neighbor voxel contains any cell that can interact with me
bool is_neighbor_voxel(Fibre* pFibre, std::vector<double> myVoxelCenter, std::vector<double> otherVoxelCenter, int otherVoxelIndex);  


extern std::unordered_map<std::string,Fibre_Definition*> fibre_definitions_by_name; 
extern std::unordered_map<int,Fibre_Definition*> fibre_definitions_by_type; 
extern std::vector<Fibre_Definition*> fibre_definitions_by_index; // works 

void display_fibre_definitions( std::ostream& os ); // done 
void build_fibre_definitions_maps( void ); // done 

Fibre_Definition* find_fibre_definition( std::string search_string ); // done 
Fibre_Definition* find_fibre_definition( int search_type );  

Fibre_Definition& get_fibre_definition( std::string search_string ); // done 
Fibre_Definition& get_fibre_definition( int search_type );  

Fibre_Definition* initialize_fibre_definition_from_pugixml( pugi::xml_node cd_node ); 
void initialize_fibre_definitions_from_pugixml( pugi::xml_node root ); 
void initialize_fibre_definitions_from_pugixml( void );

extern std::vector<double> (*fibre_division_orientation)(void);

void attach_fibres( Fibre* pFibre_1, Cell* pFibre_2 );
void detach_fibres( Fibre* pFibre_1 , Cell* pFibre_2 );

std::vector<Fibre*> find_nearby_fibres( Fibre* pFibre ); // new in 1.8.0
std::vector<Fibre*> find_nearby_interacting_fibres( Fibre* pFibre ); // new in 1.8.0

};

#endif
