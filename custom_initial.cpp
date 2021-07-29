// set the initial condition 
	


void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	// initialize BioFVM 
	
	initialize_microenvironment(); 	

    double x_coordinate_fibre = 100.0;
    //std::vector< double > position = {0.0,0.0,0.0};
    for( unsigned int n=0; n < microenvironment.number_of_voxels() ; n++ )
    {
        if (
            //dist(microenvironment.mesh.voxels[n].center,position) >= parameters.doubles["initial_tumor_rad"].value
            microenvironment.mesh.voxels[n].center[0]-x_coordinate_fibre <=1.0
            )
        {
            //microenvironment.add_dirichlet_node(n,default_microenvironment_options.Dirichlet_condition_vector); 
            microenvironment.density_vector(n) = 1.0;
        }
        else{
            microenvironment.density_vector(n) = 0.0;
        }
         
    }
	return; 
}