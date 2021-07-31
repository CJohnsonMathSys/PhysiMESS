Fibre_Potential(Cell* other_agent)
{

    double vel_adhesion = 0.03; \\ will need putting in parameters
    double vel_contact = 0.001; \\ will need putting in parameters
    double cell_velocity_max = 0.16666; \\ will need putting in parameters

    std::vector<double>cell_velocity(velocity.size(),0.0);

    if(other_agent->type_name == "fibre")
    {
        
        double velocity_dot_direction = 0.;  
        for (unsigned int j=0; j<velocity.size(); j++) 
        {
            cell_velocity += velocity[j]*velocity[j];
        }
        cell_velocity = sqrt(cell_velocity);
    
        // note that for now let's fudge fibre direction/orientation based on fibres being vertical
        double fibre_length = 1.0; // may be the unit length, notsure? options are (1,5,17) not sure
        bool isSetFibreOrientation = false;

        std::vector<double> fibre_orientation(velocity.size(),0.0);

        if(isSetFibreOrientation)
        {
            for(unsigend i-0; i<velocity.size();i++)
            {
                fibre_orientation[i] = this->state.orientation[i];
            }
        }
        else{
            // set as vertical 
            fibre_orientation[1] = 1.0;
        }

        
        
        for (unsigned int j=0; j<3; j++) 
        {
            velocity_dot_direction += fibre_orientation[j]*cell_velocity[j];   
        }
        
        for (unsigned int j=0; j<velocity.size(); j++) 
        {   
            double xip = fabs(velocity_dot_direction)/(cell_velocity+1e-8);
            double xiq = (1-xi*xi);

            // fibre_adhesion
            fibre_adhesion[j] += vel_adhesion * xip * (1 - cell_velocity/cell_velocity_max) * fibre_orientation[j];						   
            // fibre_repulsion 
            fibre_repulsion[j] += -vel_contact * xiq * velocity[j];      
        }

        // force updates velocity
        for (unsigned int j=0; j<params.dimension; j++) 
        {
            velocity[j] += (fibre_adhesion[j] + fibre_repulsion[j]);
        }
    }



    return;
}
