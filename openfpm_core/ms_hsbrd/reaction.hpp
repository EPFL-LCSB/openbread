#include "Vector/vector_dist.hpp"
#include "Vector/vector_dist_key.hpp"

#include "VCluster/VCluster.hpp"

#include "math.h"
#include <fstream>
#include <cstdlib>
#include <list>
#include <map>


#include <unistd.h>
#include <string>

#include <random>

#include <stdexcept>
#include <algorithm>    // std::random_shuffle
#include <iterator>

// Dependencies
#include "Space/SpaceBox.hpp"

#include "particles.hpp"
#include "geometry.hpp"



#define N_RANDOM_POINT 1024


#ifndef REACTION_H
#define REACTION_H

//Template


// Helper functions

// Collision test:
// test for collision of a hypothetical particle in the neightbour hood of a particle
// This function is used in first
template<typename CellList>
inline bool test_collision_reacions(td_particle_list  & _particle_list,
				   		CellList & _NN,
				   		size_t _p_key,
				   		Point<3,double> _position,
				   		double _radius)
{

  //std::cout << "Test 1 " << std::endl;
        auto Np = _NN.template getNNIterator<SAFE>(_NN.getCell(_position));



	// Postion of the hypothetical particle to be checked

	Point<3,double> xp = _position;

	// Radius of the particle to be cheked
	double r_p = _radius;



	while (Np.isNext())
	{
		// ... q
		auto q_key = Np.get();

		// skip the particle p
		if (q_key == _p_key) {++Np; continue;};

		// Get the position of q
		Point<3,double> xq = _particle_list.getPos(q_key);

		// Get the distance between p and q
		Point<3,double> r = xp - xq;

		// take the norm of this vector
		double rn = norm2(r);

		// get the radius of q
		double r_q = _particle_list.getProp<radius>(q_key);

		double t_test = (r_q + r_p)*(r_q + r_p);

		// Test for collsions or if the particle is deleted
		if (rn > t_test or _particle_list.getProp<id>(q_key) < 0 )
		  {++Np; continue;} // No Collision
		else 

		// A collision occured
		  {return true;}
	}

	// No collsion occured
	return false;
}


// Collision test:
// test collision with a singel particle!
bool test_collision_domain(td_particle_list  & _particle_list,
				   			size_t _q_key,
				   			Point<3,double> _postion,
				   			double _radius)
{
	// Postion of the hypothetical particle to be checked
	Point<3,double> xp = _postion;

	// Radius of the particle to be cheked
	double r_p = _radius;

	// skip the particle p

	// Get the position of q
	Point<3,double> xq = _particle_list.getPos(_q_key);

	// Get the distance between p and q
	Point<3,double> r = xp - xq;

	// take the norm of this vector
	double rn = norm2(r);

	// get the radius of q
	double r_q = _particle_list.getProp<radius>(_q_key);

	if (_particle_list.getProp<id>(_q_key) < 0){return false;}
	// No collision particle will be deleted

	// Test for collsions
	if (rn <= (r_q + r_p)*(r_q + r_p)) {return true;} //// A collision occured

	// No collsion occured
	return false;
}





// Zero order reaction
// To be implemented for input fluxes (reduced probablity based on the available volume to the particle that is to
// blaced in the box

struct ZerothOrderReaction
{
	int product ;
	double rate;
	int N_p;

	// Constructor
	ZerothOrderReaction(int _product, double _rate ,int _N_p){product = _product; rate = _rate; N_p = _N_p; }

	//ZerothOrderReaction(int _product, double _rate ){product = _product; rate = _rate; N_p = -1;}
	//ZerothOrderReaction(int _product, int _N_p){product = _product; rate = 0.0; N_p = _N_p;}

	~ZerothOrderReaction(){}
	template<typename CellList>
	inline bool operator()( td_particle_list & _particle_list,
													Vcluster & _v_cl,
													CellList & _NN,
													std::map<int,size_t> & _n_particles,
													std::map<int,int> & _delta_n_particles,
													std::map< int, Species* > _species_list,
													double _dt,
													std::mt19937 & _mt_rng)
	{

		if(N_p != -1) // Constant concetration >  number of particles
		{
				// Claculate the difference between current state and set velue

			int delta_n_t = N_p - _n_particles[product];

			//std::cout << "CPU:" << _v_cl.getProcessUnitID() << " " <<  delta_n_t <<  std::endl;

			// Add delta_n_t particles
			if (delta_n_t > 0)
			{
				// Go through all subdomains and try to add a particle
				auto decomp			= _particle_list.getDecomposition();

				// Number of subhyper cubes the space is divided into
				auto n_hc			= decomp.getNSubDomain();

				// List of particles to add on the current processor
				std::vector<Point<DOMAIN_DIM,float>> particles_to_add;

				// Global and local counter for the particle to add
				int add_particle = 0;
				int add_particle_loc =0;

				for(int counter = 0; (add_particle < delta_n_t  and counter < 100 * delta_n_t); counter++ )
				{

					//try yo add a particle
					bool check_next = false;

					//std::cout << counter <<  std::endl;

					//Select a subdomain randomly
					size_t hc = rand() % n_hc;

					// Get the subdomain of the processor

				    // Get the margin of the sub cubes
					SpaceBox<DOMAIN_DIM,float> subdomain_box = decomp.getSubDomain(hc);

					// Roll a point that is in the spacebox
				 	Point<DOMAIN_DIM,float> pos_1 = subdomain_box.rnd();

				    //std::cout << pos_1[0] <<  pos_1[1] << pos_1[2] << std::endl;

					// Test for all particle in the hyper cube
					auto it_dom = _particle_list.getDomainIterator();

					// Check all particles in the domain!

					while (it_dom.isNext() and !check_next and add_particle < delta_n_t )
					{

							auto p = it_dom.get();

							//std::cout  << test_collision_reacions(_particle_list, _NN, p.getKey(), pos_1, _species_list[product]->radius) << std::endl;
							// If no collision occured
							if (!test_collision_domain(_particle_list, p.getKey(), pos_1, _species_list[product]->radius))
							{
								// Check next particle
								check_next = true;
								//continue;

								// Count up the particles to add
								add_particle_loc += 1;
								add_particle = add_particle_loc;

								// Communicate the particle added
								_v_cl.sum(add_particle);
								_v_cl.execute();

								// add to list
								particles_to_add.push_back(pos_1);


							}
							else
							{
								add_particle = add_particle_loc;
								// Communicate
								_v_cl.sum(add_particle);
								_v_cl.execute();

							}
							++it_dom;
					}

					// Syncronize the proceesoes
					add_particle = add_particle_loc;

					// Communicate the particle added
					_v_cl.sum(add_particle);
					_v_cl.execute();


				}


				// Add a particle
				add_particle = add_particle_loc;
				_v_cl.sum(add_particle);
				_v_cl.execute();

				while(add_particle > delta_n_t)
				{
					int random_nr_int =  0;
					int n_cpu = _v_cl.getProcessingUnits();


					// Decide which particles to add
					// Processor n = 0 draws  a random number
					if (_v_cl.getProcessUnitID() == 0)
					{
						random_nr_int =  rand()%(n_cpu) ;
						//std::cout <<"NCPU: " <<n_cpu <<  std::endl;
						//std::cout << "Random number: " <<random_nr_int <<  std::endl;
						//std::cout << "Add X particles: " << add_particle << std::endl;
						//std::cout << "Particles to add: " << delta_n_t << std::endl;
					}
					else
					{
						random_nr_int =  0;
					}

					_v_cl.sum(random_nr_int);
					_v_cl.execute();

					if (_v_cl.getProcessUnitID() == random_nr_int and add_particle_loc > 0)
					{
						add_particle_loc--;
						add_particle = add_particle_loc;
					}
					else
					{
						add_particle = add_particle_loc;
					}

					_v_cl.sum(add_particle);
					_v_cl.execute();

				}

				// Add particles
				//

				int test_num = 0;
				for (int i=0 ; i < add_particle_loc; i++)
				{
					add_particle_pos( _particle_list , *_species_list[product] , particles_to_add[i] );
					// Update the particle list (this is on each proccesor)
					test_num++;

				}
				//std::cout << "Aded xx"<< test_num << std::endl;

				return true;

				// Create the partilce at a position

			}

			// Remove delta_n_t particles
			else if  (delta_n_t < 0)
			{
				int rm_particle = 0;
				int rm_particle_loc =0;

				_v_cl.sum(rm_particle);
				_v_cl.execute();

				//std::cout << delta_n_t <<  std::endl;

				// Particles to remove
				std::vector<size_t> particles_to_rm;
				
				//Create a product list
				auto it_dom = _particle_list.getDomainIterator();
				std::vector<size_t> product_list;
				while(it_dom.isNext() )
				  {
				    auto p = it_dom.get();
				    if (_particle_list.getProp<id>(p.getKey()) == product)
				      {
					product_list.push_back(p.getKey());
				      }
				    ++it_dom;
				  }

				//Shuffle vector
				std::shuffle(std::begin(product_list), std::end(product_list), _mt_rng);
				auto itv = product_list.begin();

				//Draw random particles from the product list
				while( -rm_particle > delta_n_t and itv != product_list.end())
				{
				        // Get particle key
					size_t _rand_particle_key = *itv;
				
					//std::cout << "Particle key " << _rand_particle_key << std::endl;

					// Add particle to the remove list
					rm_particle_loc ++;
					rm_particle = rm_particle_loc;
					_v_cl.sum(rm_particle);
					_v_cl.execute();

					particles_to_rm.push_back(_rand_particle_key);

					itv++;
				       
				}

				//std::cout << "DELTE N " << rm_particle_loc << std::endl;
				//std::cout << "DELTE N " << rm_particle << std::endl;
				//std::cout << "DELTA " << delta_n_t << std::endl;
				//std::cout << "IS NEXT" <<  it_dom.isNext()  << std::endl;

				while(-rm_particle < delta_n_t)
				{
					int random_nr_int =  0;
					int n_cpu = _v_cl.getProcessingUnits();
					// Decide which particles to add
					// Processor n = 0 draws  a random number
					if (_v_cl.getProcessUnitID() == 0)
					{
						random_nr_int =  rand()%(n_cpu);
					}
					else
					{
						random_nr_int =  0;
					}

					_v_cl.sum(random_nr_int);
					_v_cl.execute();

					if (_v_cl.getProcessUnitID() == random_nr_int && rm_particle_loc > 0)
					{
						rm_particle_loc--;
						rm_particle = rm_particle_loc;
					}
					else
					{
						rm_particle = rm_particle_loc;
					}

					_v_cl.sum(rm_particle);
					_v_cl.execute();

				}

				//std::cout << "DELTE N " << rm_particle_loc << std::endl;
				int test_num = 0;
				// Mark particles for deletion
				for (int i=0 ; i < rm_particle_loc ; i++)
				{
					// Mark particles for deletion
					_particle_list.getProp<id>(particles_to_rm[i]) = -1;
					// Update the particle list (this is on each proccesor)
					test_num--	;

					//std::cout << "Delete" << std::endl;


				}

				//std::cout << "Removed xx"<< test_num << std::endl;
				return true;

			}
			// nothing todo
			else
			{
				return false;
			}

			return false;
		}
		// Constant influx with volume rejection
		// i.e.: if a reaction takes place and it is colliding with another particle
		// reject the move
		else
		{
			// Go through all subdomains and try to add a particle
			auto decomp			= _particle_list.getDecomposition();

			// Number of subhyper cubes the space is divided into
			auto n_hc			= decomp.getNSubDomain();

			auto volume_domain = decomp.getDomain().getVolume();

			//
			// Core of the Idea of paralelization of zeroth order reaction
			//
			// Execute for all subcubes with a volume scaled rate!
			//

			bool success = false;

			// Outer loop
			for(size_t hc = 0 ; hc < n_hc ;  hc++)
			{
				// Get the margin of the sub cubes
				SpaceBox<DOMAIN_DIM,float> subdomain_box = decomp.getSubDomain(hc);

				// WTF TYPE IS THIS !!! > hc_margins
				bool check_next = false;
				// Volume of the subdomain / Total volume
				// Relative scling of the flux in to the total domain!
				double volume_scaling =  (double)subdomain_box.getVolume()/volume_domain;

				// Create a uniform distributed pseudo random number

				// Pseudo random number generator
				int seed = rand();
				std::mt19937 mt_rgn(seed);
		  		std::uniform_real_distribution<double> uniform_rnd(0.0,1.0);

				double r_i =  uniform_rnd(mt_rgn);

				// Try to create a single particle and check if it is colliding
				if(r_i <= (1 - exp( -_dt * rate * volume_scaling) ) )
				{
					// Roll a point that is in the spacebox
				    Point<DOMAIN_DIM,float> pos_1 = subdomain_box.rnd();

					// Test for all particle in the hyper cube
					auto it_dom = _particle_list.getDomainIterator();

					// Check all particles in the domain!
					while (it_dom.isNext() && !check_next)
						{
							auto p = it_dom.get();
							if (test_collision_reacions(_particle_list, _NN, p.getKey(), pos_1, _species_list[product]->radius))
							{
								// Go to the next subdomain
								check_next = true;
								//continue;
							}
							++it_dom;
						}

					// Create the partilce at pos_1
					add_particle_pos( _particle_list , *_species_list[product] , pos_1 );
					// Update the particle list
					_delta_n_particles[product]++	;
					// At least one reaction happend == sucess
					success = true;
				}

			}
			return success;

		}

	 return false;
	}

};


// First order reaction class
// struct with the information of the
// Method react is executed for every particle key
struct FirstOrderReaction
{
	// Defined from input
	int educt;
 	std::vector<int> products;
	double reaction_rate;
	// From constructor
	size_t num_prod;

	// Reaction ID
	size_t rxn_id;

	// Construtctor
	FirstOrderReaction(int _educt,  std::vector<int>_products, double _rate)
	{
		rxn_id = 0;
		educt = _educt;
		products = _products;
		reaction_rate = _rate;
		//std::cout << reaction_rate <<std::endl;
		num_prod = _products.size();
		if (num_prod > 2) {throw std::invalid_argument( "Product list exeds two products, not allowed!" );}

	}
	// Deconstructor
	~FirstOrderReaction(){}

	template<typename CellList>
	inline bool operator()(td_particle_list & _particle_list,
											   CellList & _NN,
											   std::map< int, Species* > _species_list,
											   std::map<int,int> & _n_particles,
											   size_t _p_key,
											   double _dt,
											   double _r_i)
	{

		double p_i = 1.0 - exp( -_dt * reaction_rate);
		// Try to create a single particle and check if it is colliding
		if(_r_i <= p_i)
			{
				// Try to react
				switch (num_prod)
				{
					// No products (decay)
					case 0:
					{
						// Delete the _> Problems!!!!!! do not remove particle because of ghosts
						// instead use  species = -1 to mark for deleletion !
						 //_particle_list.remove(_p_key);

						_particle_list.getProp<id>(_p_key) = -1;
						// Update the list
						_n_particles[educt]--;

						// Debugging time
						//std::cout <<  "1st order Reaction" << std::endl;
						//std::cout <<  educt << std::endl;
						//std::cout <<  "-1 "     << std::endl;

						return true;
					}


					// For product list lengt == 1
					case 1:
					{
						// Change the proprties of the particle
						_particle_list.getProp<id>(_p_key) = _species_list[products[0]]->id;
						_particle_list.getProp<radius>(_p_key) = _species_list[products[0]]->radius;
						_particle_list.getProp<diffusion>(_p_key) = _species_list[products[0]]->diff;
						_particle_list.getProp<mass>(_p_key) =_species_list[products[0]]->mass;

						// for (size_t i = 0 ; i < DOMAIN_DIM ; i++){
						// 	double dx = sqrt2Ddt*rd_normal(mt_rgn);
						// 	// Set pos0 as the position before the displacement
						// 	_particle_list.getProp<pos0>(_p_key)[i] = _particle_list.getPos(_p_key)[i];
						// 	// Displace the particle
						// 	_particle_list.getPos(_p_key)[i] += dx;
						// 	// Assign the velocity
						// 	_particle_list.getProp<velocity>(_p_key)[i] = dx/_dt;
						// }

						// Update the number of particle list
						_n_particles[educt]-- ;
						_n_particles[products[0]]++ ;


						// Debugging time
						//std::cout <<  "1st order Reaction" << std::endl;
						//std::cout <<  educt << " " <<  products[0]<< std::endl;
						//std::cout <<  "-1 "     << " 1 " << std::endl;

						return true;
					}


					// For procuct list length == 2
					case 2:
					{
						// Roll orientation vector
						Point<3,double> u_vec = random_unit_vector();

						// COM of the educt
						Point<3,double> pos_0 = _particle_list.getPos(_p_key);

						double sigma = _species_list[products[0]]->radius +  _species_list[products[1]]->radius;
						double mass_1 = _species_list[products[0]]->mass;
						double mass_2 = _species_list[products[1]]->mass;

						Point<3,double> pos_1 = pos_0 +  u_vec*(mass_2/(mass_2 + mass_1) * sigma);
						double radius_1 = _species_list[products[0]]->radius;
						// Check collsion
						if (test_collision_reacions(_particle_list, _NN, _p_key, pos_1, radius_1)){return false;}

						Point<3,double> pos_2 = pos_0 - u_vec * (mass_1/(mass_1+mass_2)) * sigma;
						double radius_2 = _species_list[products[1]]->radius;
						// Check collsion
						if (test_collision_reacions(_particle_list, _NN, _p_key, pos_2, radius_2)){return false;}

						// Change the proprties of the particle reacting to product[0]
						_particle_list.getProp<id>(_p_key) = _species_list[products[0]]->id;
						_particle_list.getProp<radius>(_p_key) = _species_list[products[0]]->radius;
						_particle_list.getProp<diffusion>(_p_key) = _species_list[products[0]]->diff;
						_particle_list.getProp<mass>(_p_key) =_species_list[products[0]]->mass;

						// Optimized directly done in propagte
						for (size_t i = 0 ; i < DOMAIN_DIM ; i++){
							_particle_list.getPos(_p_key)[i] = pos_1[i];
						}

						// Add a particle for the product[1]
						add_particle_pos( _particle_list , *_species_list[products[1]] , pos_2 );

						// Update the number of particle list
						_n_particles[educt]-- ;
						_n_particles[products[0]]++	;
						_n_particles[products[1]]++	;

						return true;
					}

				}

			}
		return false;
	}


	template<typename CellList>
	inline bool sample(td_particle_list & _particle_list,
									   CellList& _NN,
									   std::map< int, Species* > _species_list,
									   std::map<int,int> & _n_particles,
									   size_t _p_key)
	{
		// Try to react
		switch (num_prod)
		{
			// No products (decay)
			case 0:
			{
				return true;
			}

			// For product list lengt == 1
			case 1:
			{
				return true;
			}

			// For procuct list length == 2
			case 2:
			{
				// Roll orientation vector
				Point<3,double> u_vec = random_unit_vector();

				//std::cout << "u" << u_vec[0] << " " <<   u_vec[1] << " " << u_vec[2] << std::endl;

				// COM of the educt
				Point<3,double> pos_0 = _particle_list.getPos(_p_key);

				double sigma = _species_list[products[0]]->radius +  _species_list[products[1]]->radius;
				double mass_1 = _species_list[products[0]]->mass;
				double mass_2 = _species_list[products[1]]->mass;

				Point<3,double> pos_1 = pos_0 +  u_vec*(mass_2/(mass_2 + mass_1) * sigma);
				double radius_1 = _species_list[products[0]]->radius;
				
				Point<3,double> pos_2 = pos_0 - u_vec * (mass_1/(mass_1+mass_2)) * sigma;
				double radius_2 = _species_list[products[1]]->radius;

				// Check collsion
				if (test_collision_reacions(_particle_list, _NN, _p_key, pos_1, radius_1) or 
				    test_collision_reacions(_particle_list, _NN, _p_key, pos_2, radius_2)){return false;}

				return true;
			}

		}


		return false;
	}
};


// Second order reactions
struct SecondOrderReaction
{

	// Defined from input
	std::vector<int> educts;
 	std::vector<int> products;
	double reaction_rate;

	// From constructor
	size_t num_prod;
	size_t num_educ;

	SecondOrderReaction(std::vector<int> _educts,
			     std::vector<int>_products,
			     double _rate,
			     std::map< int, Species* > _species_list)
	{
		educts = _educts;
		products = _products;
		reaction_rate = _rate;

		num_prod = _products.size();
		num_educ = _educts.size();

		if (num_prod > 2) {throw std::invalid_argument( "Product list exeed two products, not allowed!" );}

	}
	// Deconstructor
	~SecondOrderReaction(){}


	// Reaction method
	template<typename CellList>
	inline bool operator()(	td_particle_list & _particle_list,
													CellList & _NN,
										   		std::map< int, Species* > _species_list,
										   		std::map<int,int> & _n_particles,
										   		size_t _p_key,
										   		size_t _q_key,
										   		double _dt,
										   		double _r_i)
	{
		// Bi mol rxn scaling
		double R =  _species_list[educts[0]]->radius + _species_list[educts[1]]->radius;

		double sig2 =4*_dt*( _species_list[educts[0]]->diff + _species_list[educts[1]]->diff);
		double sig = sqrt(sig2);

		double exp4rsig = exp(-4.0*R*R/sig2) ;

		////Calcualte monster integral
		double I_part1 = ((6.0*R*R*sig - sig2*sig + 4.0*sqrt(M_PI)*R*R*R*erfc(2.0*R/sig)));
		double I_part2 = exp4rsig*((sig2*sig)-2.0*R*R*sig);
		double I_res   = 4.0*M_PI*(I_part1 + I_part2)/12.0/sqrt(M_PI);

		//std::cout << "Hardsphere: " << 4.0*M_PI*I_res << std::endl;

		double scaled_reaction_rate = reaction_rate/I_res;

		double p_i = 1.0 - exp( -1.0 * _dt * scaled_reaction_rate);



		// Does a reaction occur acc. to the microrate law?
		if( (_r_i <= p_i) or (reaction_rate == -1) )
			{

			// std::cout <<  "2nd order Reaction " <<  std::endl;
			// std::cout <<  "R " << R <<  std::endl;
			// std::cout <<  "sig " << sig <<  std::endl;
			// std::cout <<  "prop " << p_i <<  std::endl;
			// std::cout <<  "scaled" << scaled_reaction_rate <<  std::endl;
			// std::cout <<  "_r_i" << _r_i <<  std::endl;
			
			// Log the reaction collision
			// Does a reaction occur acc. to the microrate law?
			// Difflimited

			// assign keys according to the order of the given educts (important for later)
			if (_particle_list.getProp<id>(_p_key) != educts[0])
			{
				auto temp_key = _p_key;
				_p_key = _q_key;
				_q_key = temp_key;
			}

			// Try to react
			switch (num_prod)
			{
				// No products (decay)
				case 0:
				{
					// Delete the educt

				 	//_particle_list.remove(_p_key);
				 	//_particle_list.remove(_q_key);

				 	// Mark for deletion
				 	_particle_list.getProp<id>(_p_key) = -1;
				 	_particle_list.getProp<id>(_q_key) = -1;

				 	// Update the number of particle list
					_n_particles[educts[0]]-- ;
					_n_particles[educts[1]]-- ;

					// Debugging time
					//std::cout <<  "2nd order Reaction" << std::endl;
					//std::cout <<  educts[0] << " " << educts[1] << std::endl;
					//std::cout <<  "-1 " <<  "-1 "  << std::endl;

				 	return true;
				}


				// For product list lengt == 1
				case 1:
				{
					// Roll orientation vector
					Point<3,double> u_vec = random_unit_vector();

					// COM of the educts
					Point<DOMAIN_DIM,double> pos_0_p = _particle_list.getPos(_p_key);
					Point<DOMAIN_DIM,double> pos_0_q = _particle_list.getPos(_q_key);

					double m_p = _particle_list.getProp<mass>(_p_key);
					double m_q = _particle_list.getProp<mass>(_q_key);

					// Calculation of the COM of the educts
				  Point<DOMAIN_DIM,double> pos_0 = (pos_0_p*m_p + pos_0_q*m_q )/(m_p+m_q);

					// Update the position of p to the COM
					for (size_t i = 0 ; i < DOMAIN_DIM ; i++){
						_particle_list.getPos(_p_key)[i] = pos_0[i];
					}

					// Change the proprties of the particle p to match the educt
					_particle_list.getProp<id>(_p_key) = _species_list[products[0]]->id;
					_particle_list.getProp<radius>(_p_key) = _species_list[products[0]]->radius;
					_particle_list.getProp<diffusion>(_p_key) = _species_list[products[0]]->diff;
					_particle_list.getProp<mass>(_p_key) =_species_list[products[0]]->mass;

					// Mark q for deletion
					_particle_list.getProp<id>(_q_key) = -1;

					// Update the number of particle list
					_n_particles[educts[0]]-- ;
					_n_particles[educts[1]]-- ;
					_n_particles[products[0]]++;

					// Debugging time
					//std::cout <<  "2nd order Reaction" << std::endl;
					//std::cout <<  educts[0] << " " << educts[1] << " " <<  products[0]<< std::endl;
					//std::cout <<  "-1 " <<  "-1 "  <<  " 1 " << std::endl;
					//std::cout <<  "saved key " << _p_key << std::endl;
					//std::cout <<  "delted key " << _q_key << std::endl;
					//std::cout <<  "delted id " << _particle_list.getProp<id>(_q_key)  << std::endl;

					return true;
				}

				// For procuct list length == 2
				case 2:
				{
					// Roll orientation vector
					Point<3,double> u_vec = random_unit_vector();

					// COM of the educts
					Point<3,double> pos_0_p = _particle_list.getPos(_p_key);
					Point<3,double> pos_0_q = _particle_list.getPos(_q_key);

					double m_p = _particle_list.getProp<mass>(_p_key);
					double m_q = _particle_list.getProp<mass>(_q_key);

					// Calculation of the COM of the educts
					Point<3,double> pos_0 = (pos_0_p*m_p + pos_0_q*m_q )/(m_p+m_q);

					double sigma = _species_list[products[0]]->radius +
								   _species_list[products[1]]->radius;

					double mass_1 = _species_list[products[0]]->mass;
					double mass_2 = _species_list[products[1]]->mass;

					Point<3,double> pos_1 = pos_0 + u_vec * mass_2/(mass_1+mass_2) * sigma;
					double radius_1 = _species_list[products[0]]->radius;
					// Check collsion
					if (test_collision_reacions(_particle_list, _NN, _p_key, pos_1, radius_1))
						{return false;}

					Point<3,double> pos_2 = pos_0 + u_vec * mass_2/(mass_1+mass_2) * sigma;
					double radius_2 = _species_list[products[0]]->radius;
					// Check collsion
					if (test_collision_reacions(_particle_list, _NN, _p_key, pos_1, radius_1))
						{return false;}

					// Only change the proterties of both particles to match their traget products

					// Change the proprties of the particle reacting to product[0]
					_particle_list.getProp<id>(_p_key) = _species_list[products[0]]->id;
					_particle_list.getProp<radius>(_p_key) = _species_list[products[0]]->radius;
					_particle_list.getProp<diffusion>(_p_key) = _species_list[products[0]]->diff;
					_particle_list.getProp<mass>(_p_key) =_species_list[products[0]]->mass;

					for (size_t i = 0 ; i < DOMAIN_DIM ; i++){
						_particle_list.getPos(_p_key)[i] = pos_1[i];
					}

					// Change the proprties of the particle reacting to product[1]
					_particle_list.getProp<id>(_q_key) = _species_list[products[1]]->id;
					_particle_list.getProp<radius>(_q_key) = _species_list[products[1]]->radius;
					_particle_list.getProp<diffusion>(_q_key) = _species_list[products[1]]->diff;
					_particle_list.getProp<mass>(_q_key) =_species_list[products[1]]->mass;

					for (size_t i = 0 ; i < DOMAIN_DIM ; i++){
						_particle_list.getPos(_q_key)[i] = pos_2[i];
					}


					// Update the delta number of particle list
					_n_particles[educts[0]]--   ;
					_n_particles[educts[1]]--	;
					_n_particles[products[0]]++ ;
					_n_particles[products[1]]++	;

					// Debugging time
					//std::cout <<  "2nd order Reaction" << std::endl;
					//std::cout <<  educts[0] << " " << educts[1] << " " <<  products[0]<< " " <<  products[1] << std::endl;
					//std::cout <<  "-1 " <<  "-1 "  <<  " 1 " << " 1 " << std::endl;

					return true;
				}

			}
		}

		return false;


	}
};



// The main reaction object containing all rules


struct Reactions
{
	// Properties
	std::vector<Species*> 		species_list;
	std::map< int, Species* > 	species_map;

	// Zero order reaction list
	std::vector<ZerothOrderReaction> zero_react_list;

	// First order raction list as has table to
	// for fast lookup
	std::map< int, std::vector<FirstOrderReaction> > first_react_map;
	std::map<double, double > log_accecptance;


	// Second order raction list as has table to
	// for fast lookup
	std::map< int, std::vector<SecondOrderReaction> > second_react_map;
	std::map<int, size_t > log_collisions_cpu;


	// Number of particles
	std::map< int, size_t>  n_particles;

	// Delta Vectors
	std::map< int, int> delta_n_particles;

	double dt;

	// Constructor
	Reactions(){} //empty Constructor
	Reactions(std::vector<Species*> _species_list, std::vector<size_t> _species_numbers)
		{
			// Create all possible first order maps
			std::vector<Species*>::iterator it;
			std::vector<size_t>::iterator it_n;

			it_n = _species_numbers.begin();
			it = _species_list.begin();

			while(it != _species_list.end())
			{
				std::vector<FirstOrderReaction> list;
				first_react_map[(*it)->id] = list;

				//
				// Assing number of spceies to the map species particles
				//

				n_particles[(*it)->id] = *it_n;
				delta_n_particles[(*it)->id] = 0;

				//Map species ids to species pointers
				species_map[(*it)->id] = (*it);

				it++ ;
				it_n ++;
			}

			if (_species_list.size() > 1)
			{
				// Create all possible second order maps
				std::vector<Species*>::iterator it2;
				it = _species_list.begin();

				while(it != _species_list.end())
				{
					it2 = _species_list.begin();

					while(it2 != it){
						int s_q = (*it2)->id;
						int s_p = (*it)->id;

						//std::cout << "Test make second order " << s_q << std::endl;
						//std::cout << "Test make second order " << s_p << std::endl;

						//int reaction_key;
						auto reaction_key = (s_q+s_p)*(s_q+s_p+1)+s_p;

						std::vector<SecondOrderReaction> list;

						second_react_map[reaction_key] = list;

						// Init Log collisions
						log_collisions_cpu[reaction_key] = 0;

						it2++;
					}

					// Crappy solution for it = it2
					int s_q = (*it2)->id;
					int s_p = (*it)->id;

					//std::cout << "Test make second order " << s_q << std::endl;
					//std::cout << "Test make second order " << s_p << std::endl;

					//int reaction_key;
					auto reaction_key = (s_q+s_p)*(s_q+s_p+1)+s_p;

					std::vector<SecondOrderReaction> list;

					second_react_map[reaction_key] = list;

					it++;
				}
			}

			// Assing the species list
			species_list = _species_list;

		}
	// Deconstructor
	~Reactions(){}

	// Methods

	// Add zeros order reaction
	// input flux of the system (creates particles)
	void add_zeroth_order_reaction(int _educt, double _rate, double _N)
	{
		ZerothOrderReaction reaction(_educt, _rate, _N);
		zero_react_list.push_back(reaction);
	}


	void add_first_order_reaction(int _educt, std::vector<int>_products, double _rate)
	{
		// Add the reaction to the species map
		FirstOrderReaction reaction(_educt,_products,_rate);

		// Count the id of added reactions up
		size_t num_reactions = first_react_map[_educt].size();
		reaction.rxn_id += num_reactions;

		first_react_map[_educt].push_back(reaction);

	}

	// Add second order reaction
	void add_second_order_reaction( std::vector<int>_educts,
								    std::vector<int>_products,
								   	double _rate)

	{
		if(_educts.size()>2){throw std::invalid_argument( "Educt list exeds two, this is not a second order reaction!" );}

		int s_q = _educts[0];
		int s_p = _educts[1];
		int reaction_key;


		//Make the problem sysmetric
		if(s_q < s_p)
		{ // Make Key [p,q]
			reaction_key = (s_q+s_p)*(s_q+s_p+1)+s_p;
		}
		else
		{ // Make Key [q,p]
			reaction_key = (s_q+s_p)*(s_q+s_p+1)+s_q;
		}

		//std::cout << "Reaction key" << reaction_key << std::endl;
		// Assing the reaction to the reaction list!
		SecondOrderReaction reaction(_educts,_products,_rate,species_map);
		second_react_map[reaction_key].push_back(reaction);
	}




	// Executed by each processor in every loop
	template<typename CellList>
	inline void react_zeroth_order(	td_particle_list & _particle_list,
																	Vcluster&  _v_cl,
																	CellList & _NN,
																	double _r_gskin,
																	double _dt,
																	std::mt19937 & _mt_rng)
	{

		if( zero_react_list.empty() )
		{
			return;
		}
		else
		{

			// Loop over the reaction list
			std::vector<ZerothOrderReaction>::iterator it;

			bool success = false;
			//std::cout << "Test" << std::endl;
			for (it = zero_react_list.begin(); it != zero_react_list.end() ; ++it)
			{
				auto reaction = *it;

				_particle_list.map();
				_particle_list.ghost_get<>();
				//update neigbours
				_particle_list.updateCellListSym(_NN);
				// Attempt the reaction
				success = reaction(_particle_list,_v_cl,_NN,n_particles,delta_n_particles,species_map,_dt, _mt_rng);


			}
		}

	}


	// Executed for every particle
	template<typename CellList>
	inline bool react_first_order( td_particle_list & _particle_list,
															   CellList & _NN,
															   size_t _p_key,
															   double _dt,
															   std::mt19937 & _mt_rng,
															   std::uniform_real_distribution<double> & _uniform_rnd,
															   int _n_sample = 0 )
	{

		int particle_id =  _particle_list.getProp<id>( _p_key );

		// React the possibile reactions in a random order
		//std::cout << first_react_map[ particle_id ].empty() << std::endl;

		if( first_react_map[ particle_id].empty() )
		{
			return false;
		}
		else
		{
			std::vector<FirstOrderReaction>::iterator it;
			std::vector<FirstOrderReaction> reaction_list = first_react_map[ particle_id];

			std::random_shuffle ( reaction_list.begin(), reaction_list.end() );

			//std::cout << "Particle " <<  _p_key << " ... Test" << std::endl;
			//std::cout << "ID " <<  particle_id << std::endl;
			//std::cout << 	   <<  particle_id << std::endl;

			bool success = false;
			//std::cout << "Test" << std::endl;
			for (it = reaction_list.begin(); (it != reaction_list.end() and success == false); ++it)
			{
				//std::cout << "test"<< std::endl;
				auto reaction = *it;
				// Execute the reaction if not occured
				if (not success){
					double this_rand = _uniform_rnd(_mt_rng);
					success = reaction(_particle_list,_NN,species_map,delta_n_particles,_p_key,_dt,this_rand);
				}

				// Sample the space arround the reactive particle with n
				if(_n_sample > 0)
				{
					int accepted = 0;
					for( int i = 0; i < _n_sample ; i ++)
					{
						bool is_accepted = reaction.sample(_particle_list,_NN,species_map,delta_n_particles,_p_key);

						if(is_accepted)
						{
							accepted += 1;
						}

					}
					// The acceptance log need a key that contains information about the
					// species and
					// TODO this should be initiliazed with the reaction object !!!
					double product_key = 0;
					int counter = 0;
					std::vector<int>::iterator product_it;
					for (product_it =  reaction.products.begin();
							 product_it != reaction.products.end();
						   product_it++)
					{
						product_key += (double) *product_it * std::pow(10.0,(2*counter));
						counter++;
					}

					double this_rxn_id = 1e6*(particle_id+1) + product_key;
					log_accecptance[this_rxn_id] += (double)accepted/(double)_n_sample;

					//Debuging
					//std::cout << "Out of " << _n_sample << " accepted "<< accepted << " a/n= " << log_accecptance[this_rxn_id] << std::endl;
				}



			}
			return success;
		}
	}


	// Executed for every collision
	template<typename CellList>
	inline bool react_second_order( td_particle_list & _particle_list,
															    CellList & _NN,
															    size_t _p_key,
															    size_t _q_key,
															    double _dt,
															    std::mt19937 & _mt_rng,
															    std::uniform_real_distribution<double> & _uniform_rnd
															    )
	{
		// Characterize the collision:
		int s_q = _particle_list.getProp<id>(_q_key);
		int s_p = _particle_list.getProp<id>(_p_key);
		int reaction_key;
		// Do not react with deleted particles id == -1
		if (s_q >= 0 && s_p >= 0)
		{
			//Make the problem sysmetric
			if (s_q < s_p)
			{ // Make Key [p,q]
				reaction_key = (s_q+s_p)*(s_q+s_p+1)+s_p;
			}
			else
			{ // Make Key [q,p]
				reaction_key = (s_q+s_p)*(s_q+s_p+1)+s_q;
			}

			// React the possibile reactions in a random order
			if ( second_react_map[reaction_key].empty() )
			{
				return false;
			}
			else
			{
				// do the reaction
				std::vector<SecondOrderReaction>::iterator it;
				std::vector<SecondOrderReaction> reaction_list = second_react_map[reaction_key];
				std::random_shuffle ( reaction_list.begin(), reaction_list.end() );

				bool success = false;
				for (it = reaction_list.begin(); (it != reaction_list.end() & success == false) ; ++it)
				{
					auto reaction = *it;

					// Draw a uniform radom number
					double this_rand = _uniform_rnd(_mt_rng);
					success	 = reaction(_particle_list,_NN,species_map,delta_n_particles,_p_key,_q_key,_dt,this_rand);

				}

				// log the collision attempt
				log_collisions_cpu[reaction_key] += 1;


				return success;



			}
		}

		return false;

	}


void update_n_particles(td_particle_list & _particle_list, Vcluster& _v_cl)
	{
		// summ al the deltas
		std::vector<Species*>::iterator it;
		it = species_list.begin();

		bool has_changed = false;

		std::map<int,double> this_cpu_species_number;

		// Reset the local species number
		while(it != species_list.end())
		{
			this_cpu_species_number[(*it)->id] = 0;
			it++ ;
		}

		// Count the particles in the box
		auto it_dom = _particle_list.getDomainIterator();
		while (it_dom.isNext())
		{
			auto p = it_dom.get();


			if(_particle_list.getProp<id>(p.getKey()) < 0 )
			{
				++it_dom;
				continue;
			}
			else
			{
				int this_id = _particle_list.getProp<id>(p.getKey());
				this_cpu_species_number[this_id] += 1;
			}

			++it_dom;
		}

		it = species_list.begin();
		while(it != species_list.end())
		{
			//
			// Assing number of spceies to the map species particles
			//
			int species_number = this_cpu_species_number[(*it)->id];
			//std::cout  << "Difference in particles" << delta_n << std::endl;

			_v_cl.sum(species_number);
			_v_cl.execute();


			n_particles[(*it)->id] = species_number;

			//std::cout << "Change in  species " << this_id << "  " << delta_n << std::endl;
			//std::cout << "Now species " << this_id << " " << n_particles[this_id] << std::endl;



			it++ ;
		}

	}


	//*******************************************************************************************
	// Logging functions
	//*******************************************************************************************


	// Log collisions:
	std::map<int,double> log_collisions(Vcluster& _v_cl)
	{
		std::map<int,double> all_collsions;
		std::map<int,size_t>::iterator it_map;
		for(it_map = log_collisions_cpu.begin(); it_map != log_collisions_cpu.end(); it_map++ )
			{
				size_t this_collisions =it_map->second;
				_v_cl.sum(this_collisions);
				_v_cl.execute();
				all_collsions[it_map->first] = (double) this_collisions;
				// Rest the collision counter
				it_map->second = 0;
			}
		return all_collsions;
	}


	// Log collisions TODO A nicer solution currently the acceptance is normed
	std::map<int,double> log_acceptance(Vcluster& _v_cl,
																			std::map<int,double> species_number)
	{
		std::map<int,double> all_acceptance;
		std::map<double,double>::iterator it_map;

		for(it_map = log_accecptance.begin(); it_map != log_accecptance.end(); it_map++ )
		{
			double this_accecptance = it_map->second;
			_v_cl.sum(this_accecptance);
			_v_cl.execute();
			this_accecptance /= (double)_v_cl.getProcessingUnits();

			int this_educt_id = (int)it_map->first/1e6 - 1 ;
			//int this_reaction_id = (int)it_map->first - 1e6*this_educt_id;
			// Acceptance normed by the current number of species
			all_acceptance[it_map->first] = (double) this_accecptance
																		 /(double) species_number[this_educt_id];
			// Rest the acceptance counter
			it_map->second = 0;
		}

		return all_acceptance;
	}



};


#endif //REACTION_H
