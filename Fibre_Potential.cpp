Fibre_Potential(std::vector<double>* pVelocity, Cell* thisCell, Cell* thatCell)
{

    double scalarFactor =1.0;

    std::vector<double>vectorFactor((*y).size(),1.0);

    if(thisCell->type_name == "bacteria" && thatCell->type_name == "bacteria")
    {

    }
    else if(thisCell->type_name == "bacteria" && thatCell->type_name == "fibre")
    {

    }
    else if(thisCell->type_name == "fibre" && thatCell->type_name == "bacteria")
    {

    }
    else if(thisCell->type_name == "fibre" && thatCell->type_name == "fibre")
    {

    }
    else
    {
        // wierd cell types? Do nothing


    }


    for( unsigned int i=0; i < (*y).size() ; i++ )
    {
        // dereference the velocity pointer to access ith direction
        (*pVelocity)[i] += scalarFactor * vectorFactor[i] * x[i] ; 
    }



    return;
}