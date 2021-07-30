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

#ifndef __PhysiMess_fibre_h__
#define __PhysiMess_fibre_h__

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
	// store any fibre parameters we may need 
	
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
        std::vector<Fibre*> attached_fibre; //- currently not going to use as fibres will not attach
	std::vector<Fibre*> attached_cell;

	std::vector<Fibre*> neighbors; // not currently tracked! 
	std::vector<double> orientation;
	
	double simple_pressure; 
	
	int number_of_attached_fibres( void ); //- currently not going to use fibre won't attach yet
	int number_of_attached_cells( void ); 
	
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
	
	void add_potentials(Fibre*);       // Add repulsive and adhesive forces.
	void set_previous_velocity(double xV, double yV, double zV);
	int get_current_mechanics_voxel_index();
	void turn_off_reactions(double); 		  // Turn off all the reactions of the fibre
	
	bool is_out_of_domain;
	bool is_movable;
	
	void flag_for_division( void ); //- currently not going to use fibre will not divide
	void flag_for_removal( void ); // - currently not going to use fibre will remain in domain
	
	void start_death( int death_model_index ); //- currently not going to use
	void lyse_fibre( void ); //- currently not going to use

	Fibre* divide( void ); //- currently not going to use
	void die( void ); //- currently not going to use
        void step(double dt); //- currently not going to use
	Fibre();
	
	~Fibre(); 
	
	bool assign_position(std::vector<double> new_position);
	bool assign_position(double, double, double);
	void set_total_volume(double);
	
	double& get_total_volume(void);
	
	void set_target_volume(double); //- currently not going to use fibre will be static shouldn't be needed
	void set_target_radius(double); //- currently not going to use fibre will be static shouldn't be needed
	
	void set_radius(double);
	void set_length(double);
	
	// mechanics 
	void update_position( double dt ); //- currently not going to use fibre will be static
	std::vector<double> displacement; // this should be moved to state, or made private  

	
	void assign_orientation();  // if set_orientaion is defined, uses it to assign the orientation
								// otherwise, it assigns a random orientation to the cell.
	
	void copy_function_pointers(Fibre*);
	
	void update_voxel_in_container(void);
	void copy_data(Fibre *);
	
	void ingest_fibre( Fibre* pFibre_to_eat ); //- currently not going to use unlikely to use 

	void attach_fibre( Fibre* pAddMe ); //- currently not going to use but later fibres may crosslink later
	void detach_fibre( Fibre* pRemoveMe ); //- currently not going to use but later fibres may crosslink later
	void remove_all_attached_fibres( void ); //- currently not going to use but later fibres may crosslink later

	void attach_cell( Fibre* pAddMe ); 
	void detach_cell( Fibre* pRemoveMe );  
	void remove_all_attached_cells( void ); 

	// I want to eventually deprecate this, by ensuring that 
	// critical BioFVM and PhysiCell data elements are synced when they are needed 
	
	void set_phenotype( Fibre_Phenotype& phenotype ); // do we need the Fibre flag
	void update_radius(); //- currently not going to use 
	Fibre_Container * get_container();
	
	std::vector<Fibre*>& fibres_in_my_container( void ); //- currently not going to use
	std::vector<Fibre*> nearby_fibres( void ); //- currently not going to use
	std::vector<Fibre*> nearby_interacting_fibres( void );//- currently not going to use

	std::vector<Fibre*>& cells_in_my_container( void ); 
	std::vector<Fibre*> nearby_cells( void ); 
	std::vector<Fibre*> nearby_interacting_cells( void ); 
	
	void convert_to_fibre_definition( Fibre_Definition& cd ); 
};

Fibre* create_fibre( void );  
Fibre* create_fibre( Fibre_Definition& cd );  

void delete_fibre( int ); //- currently not going to use
void delete_fibre( Fibre* ); //- currently not going to use
void save_all_fibres_to_matlab( std::string filename ); 

//function to check if a neighbor voxel contains any fibre that can interact with me 
 bool is_neighbor_voxel(Fibre* pFibre, std::vector<double> myVoxelCenter, std::vector<double> otherVoxelCenter, int otherVoxelIndex); //- currently not going to use

 //function to check if a neighbor voxel contains any cell that can interact with me
bool is_neighbor_voxel(Cell* pCell, std::vector<double> myVoxelCenter, std::vector<double> otherVoxelCenter, int otherVoxelIndex);

extern std::unordered_map<std::string,Fibre_Definition*> fibre_definitions_by_name; 
extern std::unordered_map<int,Fibre_Definition*> fibre_definitions_by_type; 
extern std::vector<Fibre_Definition*> fibre_definitions_by_index;  

void display_fibre_definitions( std::ostream& os ); 
void build_fibre_definitions_maps( void ); 

Fibre_Definition* find_fibre_definition( std::string search_string ); 
Fibre_Definition* find_fibre_definition( int search_type );  

Fibre_Definition& get_fibre_definition( std::string search_string ); 
Fibre_Definition& get_fibre_definition( int search_type );  

Fibre_Definition* initialize_fibre_definition_from_pugixml( pugi::xml_node cd_node ); 
void initialize_fibre_definitions_from_pugixml( pugi::xml_node root ); 
void initialize_fibre_definitions_from_pugixml( void );

extern std::vector<double> (*fibre_division_orientation)(void); //- currently not going to use fibre potentially may 

void attach_fibres( Fibre* pFibre_1, Fibre* pFibre_2 ); //- currently not going to use
void detach_fibres( Fibre* pFibre_1 , Fibre* pFibre_2 ); //- currently not going to use

void attach_fibres( Fibre* pFibre, Cell* pCell );
void detach_fibres( Fibre* pFibre , Cell* pCell );

std::vector<Fibre*> find_nearby_fibres( Fibre* pFibre ); //- currently not going to use
std::vector<Fibre*> find_nearby_interacting_fibres( Fibre* pFibre ); //- currently not going to use

std::vector<Fibre*> find_nearby_cells( Cell* pCell ); 
std::vector<Fibre*> find_nearby_interacting_cells( Cell* pCell ); 

};

#endif
