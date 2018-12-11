// System libs
#include <cstdlib>
#include <string>
#include <map>
// HSBRD
#include "./ms_hsbrd/wrappers.hpp"
#include "openfpm_core.hpp"

HSBRDSimulation::HSBRDSimulation(   td_species_wrapper species,
                                    td_reaction_wrapper reactions,
                                    td_initial_condition_wrapper initial_conditions,
                                    td_crowding_wrapper crowding ,
                                    td_simulation_space_warpper simulation_space,
                                    int random_seed = 1,
                                    bool is_hardsphere = true,
                                    int n_sample_1st_order = 100,
                                    bool is_constant = true
                                )
{

    _species_list           = wrap_species_list(species);
    _initial_conditions     = wrap_initial_conditions(initial_conditions,species);
    _crowding               = wrap_crowding(crowding);
    _simulation_space       = wrap_simulation_space(simulation_space);
    //TODO REMOVE TIME FROM REACTION CONSTRUCTOR!!!
    _reactions = new Reactions(*_species_list,*_initial_conditions);
    wrap_reactions(_reactions,reactions);

    _simulator = new Simulator( _species_list,
                                _initial_conditions,
                                _reactions,
                                _crowding,
                                _simulation_space,
                                random_seed,
                                is_hardsphere,
                                n_sample_1st_order,
                                is_constant);
}



HSBRDSimulation::~HSBRDSimulation(){}


// Simulation run a step with delta t
void HSBRDSimulation::step(double delta_t)
{
    _simulator->step(delta_t);
}

// return current results
td_result HSBRDSimulation::log()
{
    td_result result = _simulator->log();
    return result;
}

// Run simulation with delta t until t_max and log every t_log
td_results HSBRDSimulation::simulate(double delta_t, double t_max, int i_log)
{
    td_results results = _simulator->simulate(delta_t,t_max,i_log);
    // Deconstruct itself
    return results;
}


void HSBRDSimulation::equilibrate(double delta_t, double t_max)
{
    /*
        equilibrate simulation
    */
    _simulator->equilibrate(delta_t,t_max);
}



void py_openfpm_finalize(void)
{
    openfpm_finalize();
}


void py_openfpm_init(void)
{
    // Init virtual cluster
    int argc = 1;
    char* program_name = "python_call";
    char** argv= {&program_name,};

    openfpm_init(&argc, &argv);
}

void py_openfpm_init(char* program_name)
{
    // Init virtual cluster
    int argc = 1;
    char** argv= {&program_name,};

    openfpm_init(&argc, &argv);
}
