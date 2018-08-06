// HSBRD brd for an enzymatic reaction.
#include "../../ms_hsbrd/simulation.hpp"

// Unit System:
// 		Concentration	Diffusion 	uni-reactions	bi-reactions 		adsorption rates
// tV   10 μM 			10 μm^2s–1	1 s–1			10^5 M-1L-1			1 μm s–1

// μm-s	6000 μm–3		10 μm^2s–1	1 s–1 			1.7x10–4 μm3s–1		1 μm s–1
//
//


int main(int argc, char* argv[])
{

	openfpm_init(&argc, &argv);
	// Define the simulation
	double t_max = 1e-5;
	double dt = 1e-9; // 1.6e-9 trunc normal dist
  // Maximal possible particle size
	double max_particle_size = 10e-3;

	/*
		Simulations Space
	*/
	double box_size = 0.1;
	SimulationSpace simulation_space(box_size,PERIODIC);

	/*
		Species List
	*/
	double 	r_E  	= 3.0e-3;
	double 	r_ES 	= 3.0e-3;
	double 	r_S 	= 1.0e-3;
	double 	r_P 	= 1.0e-3;
	double 	D_E 	= 100 ;
	double 	D_ES 	= 100 ;
	double 	D_S 	= 300 ;
	double 	D_P 	= 300 ;
	double 	m_E 	= 60;
	double 	m_ES 	= 60.2;
	double 	m_S 	= 0.2;
	double 	m_P 	= 0.2;
	Species E(0,r_E,D_E,m_E);
	Species ES(1,r_ES,D_ES,m_ES);
	Species S(2,r_S,D_S,m_S);
	Species P(3,r_S,D_S,m_S);
	Species CRWD(13,0,0,0);
	// Species list
	SpeciesList species_list = {&S, &E, &ES, &P, &CRWD};

	/*
		Crowding
	*/
	double vol_fract = 0.6;
	double mu  = log(31.9);
	double sig = 0.825;
	double kbT = 1.38e-23*298.15;
	double viskosity = 1e-3;
	Crowding crowding(CRWD.id,mu,sig,vol_fract,viskosity,kbT,max_particle_size);


	/*
		Initial conditions
	*/
	size_t N_S 		= 	1350;
	size_t N_E 		= 	120;
	size_t N_ES 	= 	2;
	size_t N_P 		= 	934;
	InitialConditions initial_conditions = {N_S,N_E,N_ES,N_P,0};

	/*
		Reactions
	*/
	/* TODO: Remove the delta t from the constructor to allow dynamic time step
		 addaptions
	*/
	Reactions reactions(species_list,initial_conditions);
	//reactions.add_zeroth_order_reaction(S.id,0,N_S);
	//reactions.add_zeroth_order_reaction(E.id,0,N_E);
	//reactions.add_zeroth_order_reaction(ES.id,0,N_ES);
	//reactions.add_zeroth_order_reaction(P.id,0,N_P);
	// Reaction rates
	double k_1f = 1.01685795e-04;
	double k_1b = 4.49346405e+03;
	double k_2f = 4.49346405e+03;
	double k_2b	= 1.22023123e-04;
	// First order reactions
	reactions.add_first_order_reaction(ES.id,{E.id,S.id},k_1b);
	reactions.add_first_order_reaction(ES.id,{E.id,P.id},k_2f);
	// Second order reactions
	reactions.add_second_order_reaction({E.id,S.id}, {ES.id,} ,k_1f);
	reactions.add_second_order_reaction({E.id,P.id}, {ES.id,} ,k_2b);



	int random_seed = 1;

	Simulator simulator(&species_list,
											&initial_conditions,
										 	&reactions,
											&crowding,
											&simulation_space,
											random_seed);
  std::cout << "equilibrate" << std::endl;
  simulator.equilibrate(1e-9, 1e-7);
  std::cout << "simulate" << std::endl;
	td_results result = simulator.simulate(1e-9, 1e-7, 10);
  std::cout << "finish" << std::endl;
  openfpm_finalize();

	return 0;

}
