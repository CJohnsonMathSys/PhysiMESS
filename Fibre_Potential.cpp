Fibre_Potential(std::vector<double>* pVelocity, Cell* thisCell, Cell* thatCell)
{

    double scalarFactor =1.0; 
    double vel_adhesion = 0.03; \\ will need putting in parameters
    double vel_contact = 0.001; \\ will need putting in parameters
    double cell_velocity_max = 0.16666; \\ will need putting in parameters

    std::vector<double>vectorFactor((*y).size(),1.0);

    if(thisCell->type_name == "bacteria" && thatCell->type_name == "bacteria")
    {
        //do nothing more !!!!!!
    }
    else if(thisCell->type_name == "bacteria" && thatCell->type_name == "fibre")
    {
        double velocity_dot_direction = 0.;  
        for (unsigned int j=0; j<params.dimension; j++) 
        {
            cell_velocity += cell.velocity[j]*cell.velocity[j];
        }
        cell_velocity = sqrt(cell_velocity);
    
        // note that for now let's fudge fibre direction/orientation based on fibres being vertical
        double fibre_orientation[3];
        fibre_oreintation[0] = 0.0;
        fibre_orientation[1] = 1.0;
        fibre_orientation[2] = 0.0;
        
        for (unsigned int j=0; j<3; j++) 
        {
            velocity_dot_direction += fibre_orientation[j]*cell.velocity[j];   
        }
        
        for (unsigned int j=0; j<params.dimension; j++) 
        {   
            double xip = fabs(velocity_dot_direction)/(cell_velocity+1e-8);
            double xiq = (1-xi*xi);

            // fibre_adhesion
            fibre_adhesion[j] += vel_adhesion * xip * (1 - cell_velocity/cell_velocity_max) * fibre_orientation[j];						   
            // fibre_repulsion 
            fibre_repulsion[j] += -vel_contact * xiq * cell.velocity[j];      
        }
  }

  // force updates velocity
  for (unsigned int j=0; j<params.dimension; j++) {
    cell.vel[j] += (fibre_adhesion[j] + fibre_repulsion[j]);
  }
    
    }
    else if(thisCell->type_name == "fibre" && thatCell->type_name == "bacteria")
    {
        //do nothing at this time  !!!!!!
    }
    else if(thisCell->type_name == "fibre" && thatCell->type_name == "fibre")
    {
        //do nothing at this time  !!!!!!
    }
    else
    {
        // wierd cell types? Do nothing  !!!!!!


    }


    for( unsigned int i=0; i < (*y).size() ; i++ )
    {
        // dereference the velocity pointer to access ith direction
        (*pVelocity)[i] += scalarFactor * vectorFactor[i] * x[i] ; 
    }



    return;
}
