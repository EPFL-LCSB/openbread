// HSBRD brd for an enzymatic reaction.

// Original file
#include "Vector/vector_dist.hpp"

#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "VCluster/VCluster.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"
#include "math.h"

#include <cstdlib>

// Particle types and functions
#include "../ms_hsbrd/particles.hpp"
//Reaction operators
#include "../ms_hsbrd/reaction.hpp"
//Collision operatrors
#include "../ms_hsbrd/collision.hpp"
// Additional operators
#include "../ms_hsbrd/additional_operators.hpp"
// Logger
#include "../ms_hsbrd/logger.hpp"

// Unit System:
// 		Concentration	Diffusion 	uni-reactions	bi-reactions 		adsorption rates
// tV   10 μM 			10 μm^2s–1	1 s–1			10^5 M-1L-1			1 μm s–1

// μm-s	6000 μm–3		10 μm^2s–1	1 s–1 			1.7x10–4 μm3s–1		1 μm s–1
//
//


int main(int argc, char* argv[])
{

	/// Time and step size

	std::string folder = "0.6_crowding";
	double t_max = 1e-5;

	double dt = 1e-9; // 1.6e-9 trunc normal dist

	// Cut of radius
	double max_particle_size = 10e-3;
	double r_cut = 2.0*max_particle_size;


	// Log the time and the particle number
	std::vector<double> t_log;
	std::vector<std::vector<size_t> > n_particles_log;

	// Init virtual cluster
	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

	// domain
	float box_size = 0.1;

	Box<3,float> box({0.0,0.0,0.0},{box_size,box_size,box_size});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);


	// Total number of particles
	// Particle object Position, Properties[species_id,radius,Diffusion constant]
	td_particle_list particle_list(0,box,bc,ghost);

	size_t n_cpu = v_cl.getProcessingUnits();
	int id_cpu = v_cl.getProcessUnitID();


	//
	// Comon code to be wraped in a function
	//


	// Pseudo random number generator
	// Find the logic behind that!
	int seed = v_cl.getProcessUnitID()*8365+ 1238239;
	std::mt19937 mt_rgn(seed);

	// normal dist mean = 0 std = 1
	std::normal_distribution<double> rd_normal(0.0,1.0);
	// Create a usniform distribution
	std::uniform_real_distribution<double> uniform_rnd(0.0,1.0);



	/*
		Species  conditions
	*/


	// Parameters
	double 	r_E  	= 3.0e-3;
	double 	r_ES 	= 3.0e-3;
	double 	r_S 	= 1.0e-3;
	double 	r_P 	= 1.0e-3;

	double 	D_E 	= 100 ;
	double 	D_ES 	= 100 ;
	double 	D_S 	= 300 ;
	double 	D_P 	= 300 ;

	// Approx in kDa
	double 	m_E 	= 60;
	double 	m_ES 	= 60.2;
	double 	m_S 	= 0.2;
	double 	m_P 	= 0.2;

	// Enzyme
	Species E(0,r_E,D_E,m_E);


	// Complex
	Species ES(1,r_ES,D_ES,m_ES);


	// Metabolites
	Species S(2,r_S,D_S,m_S);
	Species P(3,r_S,D_S,m_S);


	Species CRWD(13,0,0,0);

	// Species list
	std::vector<Species*> species_list = {&S, &E, &ES, &P, &CRWD};



	// Fluid
	double kbT = 1.38e-23*298.15;
	double viskosity = 1e-3;

	// Log noramal distribution
	// Crowding (Disutribution from ecoli)
	double mu  = log(31.9);
	double sig = 0.825;

	std::lognormal_distribution<double> rd_lognormal(mu,sig);


	// Occupied volume crwd
	double vol_fract = 0.6;
	// Object
	Crowding crowding(CRWD.id,mu,sig,vol_fract,viskosity,kbT,max_particle_size);

	// Initialize crowding where 1e6 is that maximal number of particles allowed to add
	// by one processor

	size_t N_CRWD = crowding.init(particle_list,v_cl,rd_lognormal,mt_rgn,1e6);


	// Initial conditions from stoch model.
	size_t N_S 		= 	1350;
	size_t N_E 		= 	120;
	size_t N_ES 	= 	2;
	size_t N_P 		= 	934;


	// Reactions
	Reactions _reactions(species_list,{N_S,N_E,N_ES,N_P,0}, dt);

	//std::cout << "Crowding " << id_cpu << std::endl;

	// Create the Species from reactions and crowding
	//add_particles(particle_list,v_cl,_reactions);
	// Update particle numbers get the real ones

	particle_list.map();

	// Add particles
	add_particles_rnd(particle_list,v_cl,species_list,{N_S,N_E,N_ES,N_P,0},0.5);

	//Init Collsion operator
	Collisions _collisions(&_reactions);

	//std::cout << "Init reactions updated " << id_cpu << std::endl;
	// Zeroth order reactions (Keep the state constant)

	_reactions.add_zeroth_order_reaction(S.id,0,N_S);
	_reactions.add_zeroth_order_reaction(E.id,0,N_E);
	_reactions.add_zeroth_order_reaction(ES.id,0,N_ES);
	_reactions.add_zeroth_order_reaction(P.id,0,N_P);

	// Reaction rates
	double k_1f = 1.01685795e-04;
	double k_1b = 4.49346405e+03;
	double k_2f = 4.49346405e+03;
	double k_2b	= 1.22023123e-04;

	// First order reactions
	_reactions.add_first_order_reaction(ES.id,{E.id,S.id},k_1b);
	_reactions.add_first_order_reaction(ES.id,{E.id,P.id},k_2f);


	// Second order reactions
	_reactions.add_second_order_reaction({E.id,S.id}, {ES.id,} ,k_1f);
	_reactions.add_second_order_reaction({E.id,P.id}, {ES.id,} ,k_2b);

	//std::cout << "Reactions added " << id_cpu << std::endl;

	// Measurements

	std::map<int,double> species_MSD;
	std::map<int,double> species_number;

	for(int this_id = 0 ; this_id < 4 ; this_id++)
	{
		species_MSD[this_id] = 0;
		species_number[this_id] = 0;
	}


	//Timer
	timer tsim;

	// Init verlet list
	// (Ghosts need to be present ! )
	particle_list.map();
	particle_list.ghost_get<>();

	// Write the initial configuration
	size_t i = 0;
	particle_list.deleteGhost();
	particle_list.write(folder+"/particle_list",i);
	particle_list.ghost_get<>();


	// Symetric cell list
	auto NN = particle_list.getCellListSym(r_cut);


	// BRD time stepping
	tsim.start();

	// T = 0.5
	for (size_t i = 1; i < t_max/dt ; i++)
	{

		// Current time step
		double time = i*dt;

		//std::cout << "Get iterator"  <<   v_cl.getProcessUnitID()  << std::endl;
		// Get the iterator

		//
        // Update the verlet list
        particle_list.map();
        particle_list.template ghost_get<>();

		auto it3 = particle_list.getDomainIterator();

		// Propagate the particles
		//std::cout << "Propagate particles "  <<   v_cl.getProcessUnitID()  << std::endl;
		while (it3.isNext())
		{

			auto p = it3.get();

			// Reset the collision flg for both particles
			particle_list.getProp<collision_flg>(p.getKey()) = false;


			// First order reaction on the particles
			_reactions.react_first_order(particle_list,NN,p.getKey(),dt,mt_rgn,uniform_rnd,100);

			// Propagation the other particles
			double sqrt2Ddt = sqrt(2.0*particle_list.getProp<diffusion>(p)*dt);

			// here we calculate x(t + 1)s
			for(size_t j = 0; j < 3 ;j++)
			{
				double dx = sqrt2Ddt*rd_normal(mt_rgn);
				// Set pos0 as the position before the displacement
				particle_list.getProp<pos0>(p.getKey())[j] = particle_list.getPos(p.getKey())[j];
				// Displace the particle
				particle_list.getPos(p.getKey())[j] += dx;
				// Assign the velocity
				particle_list.getProp<velocity>(p.getKey())[j] = dx/dt;

			}


			++it3;
		}

        // Update the verlet list
        particle_list.map();
        particle_list.template ghost_get<>();
        // Get the Cell list structure
        particle_list.updateCellListSym(NN);


        //std::cout << time << std::endl;
		_collisions.doCollisions(particle_list,NN,v_cl,dt,mt_rgn,uniform_rnd,true);

		// Overwrite the ghosts if the collision flg is set!
		particle_list.ghost_put<replace_if_collsion_>();


		// Calculate Meansquared displacements

		auto it4 = particle_list.getDomainIterator();

		while (it4.isNext())
		{
			auto p = it4.get();

			auto this_id = particle_list.getProp<id>(p.getKey());

			double squared_disp = 0;


			for(size_t j = 0; j < 3 ; j++)
			{
				auto dx = particle_list.getProp<pos0>(p.getKey())[j] - particle_list.getPos(p.getKey())[j];

				// periodic boundary conditions
				if (dx > box_size/2.0 ){dx  -= box_size; }
				if (dx < -box_size/2.0){dx += box_size; }

				squared_disp += dx*dx;
			}

			species_MSD[this_id] += squared_disp;

			++it4;
		}

		//
		//Zeroth order reaction
		//
		// Update ghost and neighbours (because of possible collisions)
		particle_list.map();
		particle_list.ghost_get<>();
		particle_list.updateCellListSym(NN);

		// Update n_particles
		_reactions.update_n_particles(particle_list,v_cl);


		_reactions.react_zeroth_order(particle_list,v_cl,NN,r_cut,dt);
		//particle_list.map();


		// Delete the marked particles and count particles
		//auto it4 = particle_list.getDomainAndGhostIterator();
		// Clear particle numberlist
		//

		openfpm::vector<size_t> particles_to_delete;

		for(std::vector<Species*>::iterator it_species = species_list.begin() ; it_species != species_list.end() ; it_species++)
		{
			species_number[(*it_species)->id] = 0;
		}

		auto it5 = particle_list.getDomainIterator();
		while (it5.isNext())
		{
			auto p = it5.get();

			//
			if(particle_list.getProp<id>(p.getKey()) < 0 )
			{
				//std::cout <<  p.getKey() << std::endl;
				particles_to_delete.add(p.getKey());
				++it5;
				continue;
			}
			else
			{
				int this_id = particle_list.getProp<id>(p.getKey());
				species_number[this_id] += 1;
			}

			++it5;
		}

		// Delte particles
		particles_to_delete.sort();
		particle_list.remove(particles_to_delete);


		// Reodering for cache friendlyness
		if (i % int(5e2) == 0)
		{
			particle_list.reorder(5);
		}


		// Start loging after physical equlibration
		if (time > 1e-6	)
		{
			//Output log collisions
			_reactions.log_collisions(folder+"/collisions", v_cl , time );
			//Output log acceptance
			_reactions.log_acceptance(folder+"/acceptance", v_cl , time );

			// Output MSD
			log_qty_species(folder+"/mean_squared_disp", v_cl , time, species_MSD );
			// Output Species number
			log_qty_species(folder+"/species_number", v_cl , time, species_number );

		}




	}

	// Print out simulation time
	tsim.stop();
	std::cout << "Time: " << tsim.getwct() << std::endl;

	openfpm_finalize();
	return 0;
}
