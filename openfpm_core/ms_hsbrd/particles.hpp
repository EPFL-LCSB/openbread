// Simple brownian motion simulation propagate the particles
//

#include "Vector/vector_dist.hpp"

//#include "vector_dist_test.hpp"

#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "VCluster/VCluster.hpp"

#include "species.hpp"

#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "timer.hpp"
#include "math.h"

#include <cstdlib>


#ifndef PARTICLE_H
#define PARTICLE_H

constexpr int id = 0;
// Partile properties
constexpr int radius = 1;
constexpr int diffusion = 2;
constexpr int mass = 3;

// Particle integration variables
// last particle positions
constexpr int pos0 = 4;
constexpr int velocity = 5;
constexpr int force = 6;

// Identifier for the corract reaplce
constexpr int collision_flg = 7;

//Mesurements
constexpr int displacement = 8;

const size_t DOMAIN_DIM = 3;

// Redefinition of the type particle list
typedef vector_dist<DOMAIN_DIM,double,
					aggregate<int,double,double,double,
										double[DOMAIN_DIM],
										double[DOMAIN_DIM],
										double[DOMAIN_DIM],
										bool,
										double[DOMAIN_DIM]> > td_particle_list;




// Particle propagator
void propagate_particles(   td_particle_list &particle_list,
                            size_t particle_key,
                            double delta_t,
                            std::mt19937 &mt_rng,
                            std::normal_distribution<double> &rd_normal)
{
    // Propaate  the other particles
    double sqrt2Ddt
            = sqrt(2.0*particle_list.getProp<diffusion>(particle_key)*delta_t);
    // here we calculate x(t + 1)s
    for(size_t j = 0; j < 3 ;j++)
    {
        double dx = sqrt2Ddt*rd_normal(mt_rng);
        // Set pos0 as the position before the displacement
        particle_list.getProp<pos0>(particle_key)[j]
                                   = particle_list.getPos(particle_key)[j];
        // Displace the particle
        particle_list.getPos(particle_key)[j] += dx;
        // Assign the velocity
        particle_list.getProp<velocity>(particle_key)[j] = dx/delta_t;
    }

}



// Crowding
struct Crowding
{
	int n_crowders;
	int crwd_id; // species id for the particles
	double sig; // parameters of the log-normal distribution for the mass
	double mu;
	double max_r; // Maximal radius in mum
	// Volume occupacy of the crowders
	double phi;
	int seed;
	//Fluid
	double kbT;
	double visk;

    //Constrauctros
    Crowding(){}
	Crowding(	int _id,
				double _mu,
				double _sig,
				double _phi,
				double _visk,
				double _kbT,
				double _max_r)
	{
		crwd_id = _id;
		sig =_sig;
		mu = _mu;
		phi = _phi;
		max_r = _max_r;

		visk = _visk;
		kbT  = _kbT;

		n_crowders = 0;

	}

	~Crowding(){}

	double add_crowder( td_particle_list & _particle_list,
		 									std::mt19937 & _mt_rgn ,
											std::lognormal_distribution<double> & _rd_lognormal )
	{

		// Init log normal distribution
		//std::lognormal_distribution<double>  rd_lognormal(mu,sig);
        double m_rnd;

		// Draw the molcular weight (in kDa)
        if (sig > 0)
        {
            m_rnd = _rd_lognormal(_mt_rgn);
        }
        else
        {
            m_rnd = exp(mu);
        }


        // Calculate the hydrodynamic radius in nm from "Biologist formula"
        // Kalwarczyk, T., M. Tabaka and R. Holyst (2012) Bioinformatics 28(22)
        // Mass is in Da (Dist in kDa)!!!
		// Formula is in nanometer thus x1e-3 to get mum
		double r_rnd = pow((m_rnd*1000.0),0.392)*0.0515*1e-3;

		// If particle is larger than max radius
		if (r_rnd > max_r)
		{
			double sc = max_r/r_rnd;
			r_rnd = max_r;
			// Scale mass
			m_rnd = m_rnd*sc*sc*sc;

		}

		// Viskosity and kbT in Si Units ..
		// r_ranf in mum
		// Diffusion output mum^2/s
		// Proper scaling TBD
		//
		double diff_rnd = kbT/(6.0*M_PI*r_rnd*1e-6*visk)*pow(10.0,12);

		SpaceBox<DOMAIN_DIM,float> box = _particle_list.getDecomposition().getDomain();

		// Volume raction
		double phi = r_rnd*r_rnd*r_rnd*M_PI*4.0/3.0/box.getVolume();

		// Random positon
		auto rnd_pos = box.rnd();

		// Add the particle
		_particle_list.add();

		// Dimensional poroperties
		for(size_t i = 0 ; i < DOMAIN_DIM ; i++)
		{
			_particle_list.getLastPos()[i] = rnd_pos[i];
		}
		// Dimensional poroperties
		for(size_t i = 0 ; i < DOMAIN_DIM ; i++)
		{
			_particle_list.template getLastProp<pos0>()[i] = 0.0;
			_particle_list.template getLastProp<velocity>()[i] = 0.0;
			_particle_list.template getLastProp<displacement>()[i] = 0.0;
		}

		_particle_list.template getLastProp<id>() = crwd_id;
		_particle_list.template getLastProp<radius>() = r_rnd;
		_particle_list.template getLastProp<diffusion>() = diff_rnd;
		_particle_list.template getLastProp<mass>() = m_rnd;

		_particle_list.template getLastProp<collision_flg>() = false;

		// Output is the contribution to volume fraction of this particle
		return phi;

	}

	// initialize the corwing speciefied
	size_t init( td_particle_list & _particle_list,
				 Vcluster & _v_cl,
				 std::lognormal_distribution<double>  & _rd_lognormal,
				 std::mt19937 & _mt_rgn,
				 int n_max)

	{

		int n_cpu = _v_cl.getProcessingUnits();
		// Total volume fraction added
		double phi_total = 0.0;
		double phi_total_loc = 0.0;

		size_t n_total = 0;

		openfpm::vector<double> v;
   		_v_cl.allGather(phi_total_loc,v);
    	_v_cl.execute();


		// Add particles until reached occupied volume fraction is reached
		for(int counter = 0 ;  (phi_total < (float)phi/n_cpu  && counter < n_max) ; counter++ )
		{
			// Clear

			 v.clear();
			// Add a crowder
			//std::cout << "Add particle  " << id_cpu << std::endl;
			phi_total_loc += add_crowder(_particle_list,_mt_rgn,_rd_lognormal);
			//std::cout << "Added particle " << id_cpu <<std::endl;
			phi_total = 0;


			n_total++;

   			_v_cl.allGather(phi_total_loc,v);

    		_v_cl.execute();

    		for(size_t i = 0 ; i < v.size() ; i++)
    		{
    			phi_total += v.get(i);
    		}

			//std::cout << "end  " << phi_total << "  " << id_cpu << std::endl;

		}
		//std::cout << "return " << phi_total << std::endl;

   		_v_cl.sum(n_total);
    	_v_cl.execute();

		return n_total;

	}

};





// Particle creators
void add_n_particles_no_check(td_particle_list & _particle_list,
					 										Species _species,
					 										size_t N_p)
{

	for(size_t i = 0 ; i < N_p ; i++ )
	{
		// Create a random point inside the SpaceBox of the Domain
		SpaceBox<DOMAIN_DIM,float>  box = _particle_list.getDecomposition().getDomain();
		auto rnd_pos = box.rnd();

		_particle_list.add();

		// Dimensional poroperties
		for(size_t i = 0 ; i < DOMAIN_DIM ; i++)
		{
			_particle_list.getLastPos()[i] = rnd_pos[i];
		}
		// Dimensional poroperties
		for(size_t i = 0 ; i < DOMAIN_DIM ; i++)
		{
			_particle_list.template getLastProp<pos0>()[i] = 0.0;
			_particle_list.template getLastProp<velocity>()[i] = 0.0;
			_particle_list.template getLastProp<displacement>()[i] = 0.0;
		}

		_particle_list.template getLastProp<id>() = _species.id;
		_particle_list.template getLastProp<radius>() = _species.radius;
		_particle_list.template getLastProp<diffusion>() = _species.diff;
		_particle_list.template getLastProp<mass>() = _species.mass;

		_particle_list.template getLastProp<collision_flg>() = false;

		// Roll a random position in space !!! Box size has to go here
		// Go through all subdomains and try to add a particle
		//std::cout << "Create Particle position: "<<  _particle_list.getLastPos()[0] << std::endl;
	}


};

// Add particles rnd (no )
void add_particles_rnd(td_particle_list & _particle_list,
					  Vcluster & _v_cl,
					  std::vector<Species*> _species_list,
					  std::vector<size_t> _species_numbers,
					  double _dt)
{
	//Add particles with out checking for collisions
	std::vector<Species*>::iterator it;
	std::vector<size_t>::iterator it_n;

	it_n = _species_numbers.begin();
	it = _species_list.begin();

	int n_cpu = _v_cl.getProcessingUnits();

	int id_cpu = _v_cl.getProcessUnitID();

	double min_r = FLT_MAX;
	double max_r = 0.0;

	size_t N_total = 0;

	while(it != _species_list.end())
	{
		// Eache processor adds  number of particles / n_cpu
		add_n_particles_no_check(_particle_list,(**it), *it_n / n_cpu);

		// First processor adds the rest i.e. number of particles % n_cpu
		if (id_cpu == 0)
		{
			add_n_particles_no_check(_particle_list,(**it), *it_n % n_cpu);
		}

		// get the maximal radius
		if (max_r < (*it)->radius)
		{
			max_r = (*it)->radius;
		}
		if ((*it)->radius < min_r)
		{
			min_r = (*it)->radius;
		}

		N_total += *it_n;

		it++ ;
		it_n ++;
	}

}


// Function used in reactions
void add_particle_pos(td_particle_list & _particle_list,
					 Species _species,
					 Point<3,double> _position)
{
	_particle_list.add();

	_particle_list.template getLastProp<id>() = _species.id;
	_particle_list.template getLastProp<radius>() = _species.radius;
	_particle_list.template getLastProp<diffusion>() = _species.diff;
	_particle_list.template getLastProp<mass>() = _species.mass;

	_particle_list.template getLastProp<collision_flg>() = false;

	for (size_t i = 0; i < DOMAIN_DIM; i++ )
	{
		_particle_list.getLastPos()[i]  = _position[i];
		_particle_list.template getLastProp<pos0>()[i] = 0.0;
		_particle_list.template getLastProp<velocity>()[i] = 0.0;
		_particle_list.template getLastProp<displacement>()[i] = 0.0;
	}

}



#endif //PARTICLE_H
