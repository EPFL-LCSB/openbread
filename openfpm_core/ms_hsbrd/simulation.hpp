#include <cstdlib>

// OPENFPM
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "VCluster/VCluster.hpp"
#include "math.h"

// HSBRD
#include "particles.hpp"
#include "reaction.hpp"
#include "collision.hpp"
#include "additional_operators.hpp"


 #ifndef SIMULATION_H
 #define SIMULATION_H

// Abstract
typedef std::vector<Species*> SpeciesList;
typedef std::vector<size_t> InitialConditions;

typedef std::map<std::string, std::map<int, double> > td_result;
typedef std::vector<std::map<std::string, std::map<int, double> > > td_results;
//typedef std::vector<double> td_time;
//typedef std::tuple<td_time,td_results> td_simulation_result;
// simulations spce
class SimulationSpace {
    //
public:
    size_t boundary_conditions[3];
    Box<3,float> box;
    float box_size;
    SimulationSpace();
    SimulationSpace(float size, size_t _type);
};

// Constructor for quadratic box
SimulationSpace::SimulationSpace(){} // EMPTY CONSTRUCTOR
SimulationSpace::SimulationSpace(float size, size_t _type)
{
    box_size =  size;
    boundary_conditions[0] = _type;
    boundary_conditions[1] = _type;
    boundary_conditions[2] = _type;
    box = Box<3,float>({0.0,0.0,0.0},{size,size,size});
}




// Simulator object

class Simulator {
    //
public:
    td_particle_list* particle_list;
    double time;

    SpeciesList* species_list;
    InitialConditions *initial_conditions;
    Crowding* crowding;
    Reactions* reactions;
    Collisions* collisions;

    SimulationSpace* simulation_space;
    Vcluster *virtual_cluster;

    bool is_hardsphere;
    bool is_constant;

    CellList<3, double, Mem_fast, shift<3, double>,
                openfpm::vector<long unsigned int> > nearest_neighbours;
    Ghost<3,float> ghost;

    // Random numbers
    std::mt19937 random_number_generator;
    std::normal_distribution<double> normal_distributed_random_number;
    std::uniform_real_distribution<double> uniform_distributed_random_number;
    std::lognormal_distribution<double> log_normal_distributed_random_number;

    //
    double cutoff_radius;
    // Sampling
    int n_samples_1st_order;

    // log speces_numbers
    std::map<int,double> species_number;
    std::map<int,double> species_number_prev;

    Simulator();
    Simulator(  SpeciesList *in_species_list,
                InitialConditions *in_initial_conditions,
                Reactions *in_reactions,
                Crowding *in_crowding,
                SimulationSpace *in_simulation_space,
                int random_seed,
                bool _is_hardspher,
                int _n_sample_1st_order,
                bool _is_constant);


    virtual void step(double delta_t, bool is_reactive );
    virtual td_result log();
    virtual td_results simulate( double delta_t,
                                           double t_max,
                                           int i_log );
    virtual void equilibrate( double delta_t,
                              double t_max );

};

// Empty class
Simulator::Simulator(){}

Simulator::Simulator( SpeciesList *in_species_list,
                      InitialConditions *in_initial_conditions,
                      Reactions *in_reactions,
                      Crowding *in_crowding,
                      SimulationSpace *in_simulation_space,
                      int random_seed,
                      bool _is_hardsphere = true,
                      int  _n_sample_1st_order = 10,
                      bool _is_constant = true)
{

    species_list = in_species_list;
    initial_conditions = in_initial_conditions;
    reactions = in_reactions;
    crowding = in_crowding;
    simulation_space = in_simulation_space;

    is_hardsphere = _is_hardsphere;
    is_constant = _is_constant;
    n_samples_1st_order = _n_sample_1st_order;

    cutoff_radius = crowding->max_r*2.0;

    // create the virutal clsuter
    virtual_cluster = &create_vcluster();

    ghost = Ghost<3,float>(cutoff_radius);

    // init particle_list
    particle_list = new td_particle_list(0, simulation_space->box,
                                        simulation_space->boundary_conditions,
                                        ghost);

    nearest_neighbours = particle_list->getCellListSym(cutoff_radius);

    int parallel_seed = (virtual_cluster->getProcessUnitID()+random_seed)*8365+1238239;

    // Init random numbers
    random_number_generator = std::mt19937(parallel_seed);
        normal_distributed_random_number
            = std::normal_distribution<double>(0.0,1.0);
    uniform_distributed_random_number
            = std::uniform_real_distribution<double>(0.0,1.0);
    log_normal_distributed_random_number
            = std::lognormal_distribution<double>(crowding->mu, crowding->sig);



    crowding->init( *particle_list,
                    *virtual_cluster,
                    log_normal_distributed_random_number,
                    random_number_generator,
                    1e6);


    // Add particles
    particle_list->map();
    add_particles_rnd(  *particle_list,
                        *virtual_cluster,
                        *species_list,
                        *initial_conditions,
                        0.5);

    // Init collsion operaor
    collisions = new Collisions(reactions);

    particle_list->map();
    particle_list->ghost_get<>();

    time = 0.0;


    // Init species numbers
    auto it5 = particle_list->getDomainIterator();
    while (it5.isNext())
    {
        auto p = it5.get();
        int this_id = particle_list->getProp<id>(p.getKey());
        species_number[this_id] += 1.0;
        species_number_prev[this_id] += 1.0;
        ++it5;
    }


}



void Simulator::step(double delta_t, bool _is_reactive = true)
{
    // Map particles and get ghosts
    particle_list->map();
    particle_list->template ghost_get<>();

	auto it3 = particle_list->getDomainIterator();

    // Propagate the particles
    //std::cout << "Propagate particles "  <<   v_cl.getProcessUnitID()  << std::endl;
    while (it3.isNext())
    {
        auto particle = it3.get();
        // Reset the collision flg for both particles
        particle_list->getProp<collision_flg>(particle.getKey()) = false;
        // First order reaction on the particles
        if (_is_reactive)
        {
            reactions->react_first_order(*particle_list,
                                         nearest_neighbours,
                                         particle.getKey(),
                                         delta_t,
                                         random_number_generator,
                                         uniform_distributed_random_number,
                                         n_samples_1st_order);
        }

        // Propagate particles
        propagate_particles(*particle_list,
                            particle.getKey(),
                            delta_t,
                            random_number_generator,
                            normal_distributed_random_number);
        ++it3;
    }

    if (_is_reactive and is_constant)
    {
        reactions->react_zeroth_order(*particle_list,
                                      *virtual_cluster,
                                      nearest_neighbours,
                                      cutoff_radius,
                                      delta_t);
    }

    // Update the verlet list
    particle_list->map();
    particle_list->template ghost_get<>();
    // Get the Cell list structure
    particle_list->updateCellListSym(nearest_neighbours);

	collisions->doCollisions(*particle_list,
                             nearest_neighbours,
                             *virtual_cluster,
                             delta_t,
                             random_number_generator,
                             uniform_distributed_random_number,
                             is_hardsphere,
                             _is_reactive);
    // Overwrite ghosts and update
	particle_list->ghost_put<replace_if_collsion_>();
    particle_list->map();
    particle_list->ghost_get<>();
    particle_list->updateCellListSym(nearest_neighbours);
    reactions->update_n_particles(*particle_list,
                                  *virtual_cluster);
    if (_is_reactive)
    {
        reactions->react_zeroth_order(*particle_list,
                                      *virtual_cluster,
                                      nearest_neighbours,
                                      cutoff_radius,
                                      delta_t);
    }

    openfpm::vector<size_t> particles_to_delete;

    for(SpeciesList::iterator it_species = species_list->begin() ;
                             it_species != species_list->end() ;
                             it_species++)
    {
        species_number_prev[(*it_species)->id] = species_number[(*it_species)->id];
        species_number[(*it_species)->id] = 0;

    }
    auto it5 = particle_list->getDomainIterator();
    while (it5.isNext())
    {
        auto p = it5.get();

        //Fetch all particles to be deleted
        if(particle_list->getProp<id>(p.getKey()) < 0 )
        {
            //std::cout <<  p.getKey() << std::endl;
            particles_to_delete.add(p.getKey());
            ++it5;
            continue;
        }
        else
        {
            int this_id = particle_list->getProp<id>(p.getKey());
            species_number[this_id] += 1.0;
        }

        ++it5;
    }

    // Delte particles
    particles_to_delete.sort();
    particle_list->remove(particles_to_delete);

}

td_result Simulator::log()
{
    td_result result;
    result["species"]   = species_number;
    result["acceptance"] = reactions->log_acceptance(*virtual_cluster,species_number_prev);
    result["collisions"] = reactions->log_collisions(*virtual_cluster);

    return result;
}

td_results Simulator::simulate( double delta_t,
                                double t_max,
                                int i_log)
{
    td_results results;

    for (size_t i = 1; i < t_max/delta_t ; i++)
    {
        time = i*delta_t;
        step(delta_t,true);

        // Reorder particle list for chache friendlyness
        if (i % int(5e2) == 0)
        {
            particle_list->reorder(5);
        }
        if (i % i_log == 0)
        {
            td_result result = log();
            std::map<int,double> log_time;
            log_time[0]=time;
            result["time"] = log_time;
            results.push_back(result);
        }

    }
    return results;
}


void Simulator::equilibrate( double delta_t,
                             double t_max )
/*
    Run until time T without any reaction to obtain an
    equilibrated simulation environment
*/
{
    for (size_t i = 1; i < t_max/delta_t ; i++)
    {
        time = i*delta_t;
        // Run
        step(delta_t,false);
        // Reorder particle list for chache friendlyness
        if (i % int(5e2) == 0)
        {
            particle_list->reorder(5);
        }
    }
}

 #endif //SIMULATION_H
