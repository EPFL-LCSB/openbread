
#include <cstdlib>
#include <map>
#include <vector>
#include <string>

//#include "./reaction.hpp"
//#include "./particles.hpp"
//#include "./species.hpp"
#include "./simulation.hpp"

#ifndef WRAPPERS_H
#define WRAPPERS_H


// Define SpciesList derived derived from Species
typedef std::vector<Species*> SpeciesList;
typedef std::vector<size_t> InitialConditions;

// Wrappers
typedef std::map< std::string , std::map<std::string, double > > td_species_wrapper;
typedef std::map< std::string , std::map<std::string, std::vector<double> > > td_reaction_wrapper;
typedef std::map< std::string ,int > td_initial_condition_wrapper;
typedef std::map<std::string, double > td_crowding_wrapper;
typedef std::map<std::string, double > td_simulation_space_warpper;



SpeciesList* wrap_species_list(td_species_wrapper species){
	auto species_list = new SpeciesList() ;

	td_species_wrapper::iterator this_species;

	for ( this_species = species.begin(); this_species!= species.end(); this_species++)
	{
		std::string name = this_species->first;
		std::map<std::string, double > species_data = this_species->second;

		int    id = (int) species_data["id"];
		double diff_constant = species_data["diff"];
		double coll_radius = species_data["rad"];
		double mass = species_data["mass"];
		// Construct species object
		Species* new_species = new Species(id, coll_radius, diff_constant, mass);

		// add to species list
		species_list->push_back(new_species);

	}
	return species_list;
}


// wrap the reactions from python dict to
//Second order
// {ReactionName : {Educts: [0,1], Products: [0.0], rate: [0.1]} }
//First order
// {ReactionName : {Educts: [0], Products: [0.0], rate: [0.1]} }
//Zeroth order const concetration
// {ReactionName : {Educts: [], Products: [0.0], rate: [0,1000]} }
//Zeroth order const flux
// {ReactionName : {Educts: [], Products: [0.0], rate: [0.1]} }

typedef std::map<std::string, std::vector<double>> td_reaction_data;

void wrap_reactions(Reactions* _reactions, td_reaction_wrapper reactions_dict)
{ td_reaction_wrapper::iterator this_reaction;
	for (this_reaction = reactions_dict.begin();
			 this_reaction != reactions_dict.end();
		 	 this_reaction++)
	{
		std::string name = this_reaction->first;
		td_reaction_data reaction_data = this_reaction->second;
		auto educts = reaction_data["educts"];
		auto products = reaction_data["products"];
		auto rates = reaction_data["rates"];

		if(educts.size() == 2){
			std::vector<int> int_educts(educts.begin(),educts.end());
			std::vector<int> int_products(products.begin(),products.end());
			double rate_constant = rates[0];
			_reactions->add_second_order_reaction(int_educts,int_products,rate_constant);
		}
		else if(educts.size() == 1){
			int int_educt = (int) educts[0];
			std::vector<int> int_products(products.begin(),products.end());
			double rate_constant = rates[0];
			_reactions->add_first_order_reaction(int_educt,int_products,rate_constant);
		}
		else if(educts.empty()){
			int int_product = (int) products[0];
			if (rates.size() == 2){
				double number_of_particles = rates[1];
				_reactions->add_zeroth_order_reaction(int_product,0,number_of_particles);
			}
			else if(rates.size() == 1){
				double rate_constant = rates[0];
				_reactions->add_zeroth_order_reaction(int_product,rate_constant,-1);
			}
		}


	}

}

// wrap the crowding from python dict
// {"mu": parameter, "sigma": reactions_dict }
Crowding* wrap_crowding( td_crowding_wrapper crowding_dict)
{
	int crowding_id = (int) crowding_dict["id"];
	double mu =  crowding_dict["mu"];
	double sigma =  crowding_dict["sigma"];
	double vol_fract = crowding_dict["volume_fraction"];
	double viskosity = crowding_dict["viscosity"];
	double kbT = crowding_dict["kBT"];
	double max_particle_size = crowding_dict["max_radius"];
	// {"mu": parameter, "sigma": reactions_dict }
	auto crowding = new Crowding (	crowding_id,
																	mu,
																	sigma,
																	vol_fract,
																	viskosity,
																	kbT,
																	max_particle_size);
	return crowding;
}


InitialConditions* wrap_initial_conditions(td_initial_condition_wrapper init_con_dict,
																					td_species_wrapper species)
{
	auto initial_conditions= new InitialConditions();
	td_species_wrapper::iterator this_species;

	for ( this_species = species.begin(); this_species!= species.end(); this_species++)
	{
		std::string name = this_species->first;
		//initial_conditions = this_species->second;

		initial_conditions->push_back(init_con_dict[name]);
	}

	return initial_conditions;
}

SimulationSpace* wrap_simulation_space(td_simulation_space_warpper space_dict)
{
	float box_size = space_dict["box_size"];
	auto simulation_space = new SimulationSpace(box_size,PERIODIC);
	return simulation_space;
}

#endif //WRAPPERS_H
