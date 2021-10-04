Fibre_Potential(Cell* other_agent)
{

    // update the velocity of the cell if the other agent is a fibre
    std::vector<double>cell_velocity(velocity.size(),0.0);

    if(other_agent->type_name == "fibre")
    {
        
        double velocity_dot_direction = 0.;  
        for (unsigned int j=0; j<velocity.size(); j++) 
        {
            cell_velocity += velocity[j]*velocity[j];
        }
        cell_velocity = sqrt(cell_velocity);
    
        // note that for now let's fudge fibre direction/orientation based on fibres being vertical (is this still the case? - Connah (4th October 2022))
        
        for (unsigned int j=0; j<3; j++) 
        {
            velocity_dot_direction += other_agent->state.orientation[i]*cell_velocity[j];   
        }
        
        for (unsigned int j=0; j<velocity.size(); j++) 
        {   
            double xip = fabs(velocity_dot_direction)/(cell_velocity+1e-8);
            double xiq = (1-xi*xi);

            // fibre_adhesion
            fibre_adhesion[j] += other_agent->parameters.mVelocityAdhesion * xip * (1 - cell_velocity/this->parameters.mCellVelocityMaximum) * other_agent->state.orientation[i];						   
            // fibre_repulsion 
            fibre_repulsion[j] += -other_agent->parameters.mVelocityContact * xiq * velocity[j];      
        }

        // force updates velocity
        for (unsigned int j=0; j<params.dimension; j++) 
        {
            velocity[j] += (fibre_adhesion[j] + fibre_repulsion[j]);
        }
    }

    return;
}