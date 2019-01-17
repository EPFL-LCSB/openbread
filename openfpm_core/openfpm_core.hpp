#include <cstdlib>
#include <string>
#include <map>
#include <vector>

#ifndef PYHSBRD_H
#define PYHSBRD_H
// Abstract type definitions
typedef std::map< std::string , std::map<std::string, double > > td_species_wrapper;
typedef std::map< std::string , std::map<std::string, std::vector<double> > > td_reaction_wrapper;
typedef std::map< std::string ,int > td_initial_condition_wrapper;
typedef std::map<std::string, double > td_crowding_wrapper;
typedef std::map<std::string, double > td_simulation_space_warpper;


class Species;
typedef std::vector<Species*> SpeciesList;
class Reactions;
typedef std::vector<size_t> InitialConditions;
class Crowding;
class Simulator;
class SimulationSpace;

typedef std::map<std::string, std::map<int, double> > td_result;
typedef std::vector<std::map<std::string, std::map<int, double> > > td_results;
typedef std::vector<double> td_time;
//typedef std::tuple<td_time,td_results> td_simulation_result;


class HSBRDSimulation {
public:
    HSBRDSimulation(    td_species_wrapper species,
                        td_reaction_wrapper reactions,
                        td_initial_condition_wrapper initial_conditions,
                        td_crowding_wrapper crowding,
                        td_simulation_space_warpper simulation_space,
                        int random_seeed,
                        bool is_hardsphere,
                        int n_sample_1st_order,
                        bool is_constant
                    );

    ~HSBRDSimulation();

    void step(double delta_t);
    td_result log();
    td_results simulate(double delta_t,
                        double t_max,
                        int i_log,
                        bool is_reactive);
    void equilibrate(  double delta_t,
                       double t_max);


private:
    SpeciesList*         _species_list;
    InitialConditions*   _initial_conditions;
    Reactions*           _reactions;
    Crowding*            _crowding;
    SimulationSpace*     _simulation_space;
    Simulator* _simulator;

};


void py_openfpm_init(void);
void py_openfpm_init(char* _id);

void py_openfpm_finalize(void);


#endif //
