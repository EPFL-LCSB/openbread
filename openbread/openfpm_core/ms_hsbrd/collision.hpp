#include "Vector/vector_dist.hpp"
#include "Vector/vector_dist_key.hpp"

#include "math.h"
//
#include <mpi.h>

#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "VCluster/VCluster.hpp"
#include "Packer_Unpacker/Packer_util.hpp"

#include <cstdlib>
#include <list>
#include <vector>

#include <typeinfo>

#include "reaction.hpp"

#ifndef COLLISION_H
#define COLLISION_H

// Hardsphere elastic factor
const int k_fact = 1;


// Collsion detection
struct coll
{
	size_t key_q;
	size_t key_p;
	double dist;
	//Constructor
	coll(size_t _key_q, size_t _key_p, double _dist){key_q = _key_q; key_p = _key_p; dist = _dist;}
	//Deconstructor
	~coll(){}

	// to sort a list of coll's
	bool operator<(coll const &other){return dist < other.dist;}

};



// Hard sphere collision dunciton compputes the new particle positions
// from a detected collision
bool hardSphereCollision(td_particle_list  & _particle_list,
				   		size_t q,
				   		size_t p,
				   		double _dt)
{

	//std::cout << q << std::endl;
	//std::cout << p << std::endl;
	//std::cout << _particle_list.getProp<mass>(p) << std::endl;

	double mq = 2.0*_particle_list.getProp<mass>(q)/(_particle_list.getProp<mass>(q)+_particle_list.getProp<mass>(p));
	double mp = 2.0*_particle_list.getProp<mass>(p)/(_particle_list.getProp<mass>(q)+_particle_list.getProp<mass>(p));

	// velocities before the collsion
	Point<3,double> xq = _particle_list.getProp<pos0>(q);
	Point<3,double> xp = _particle_list.getProp<pos0>(p);


	Point<3,double> vq = _particle_list.getProp<velocity>(q);
	Point<3,double> vp = _particle_list.getProp<velocity>(p);

	double rq = _particle_list.getProp<radius>(q);
	double rp = _particle_list.getProp<radius>(p);

	double sig2 = (rq+rp)*(rq+rp);

	Point<3,double> dx =  xp - xq;
	Point<3,double> dv =  vp - vq;

	//
	double XX = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	double VV = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
	double XV = dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2];

	if (XV >= 0 ){return false;}

	//std::cout << "XX " << XX << std::endl;
	//std::cout << "VV " << XX << std::endl;
	//std::cout << "XV " << XX << std::endl;

	double det = XV*XV - (XX - sig2)*VV;

	if (det <= 0 ){return false;}

	// Collsion time
	double collision_time = -(sqrt(det)+XV)/VV;

	if (collision_time > _dt){return false;}

	// distance vector at the time of collision
	Point<3,double> dx_tc = dx+dv*collision_time;


	//helper vars
	double dx_tc_2  = dx_tc[0]*dx_tc[0] + dx_tc[1]*dx_tc[1] + dx_tc[2]*dx_tc[2];
	double dx_tc_dv = dx_tc[0]*dv[0]    + dx_tc[1]*dv[1]    + dx_tc[2]*dv[2];

	// Velocity after the collision
	Point<3,double> vq_tc = vq + dx_tc*(mp*(dx_tc_dv/dx_tc_2))*k_fact;
	Point<3,double> vp_tc = vp - dx_tc*(mq*(dx_tc_dv/dx_tc_2))*k_fact;

	//std::cout << "w/o corr " << std::endl;
	//std::cout << "xq " << xp[0] << std::endl;
	//std::cout << "xp " << xq[0] << std::endl;

	// Postion fter the collision
	xq = xq +  collision_time*vq + (_dt - collision_time) * vq_tc ;
	xp = xp +  collision_time*vp + (_dt - collision_time) * vp_tc ;

	//std::cout << "corr " << std::endl;
	//std::cout << "xq " << xp[0] << std::endl;
	//std::cout << "xp " << xq[0] << std::endl;

	// Assing the values to the particles
	_particle_list.getPos(p)[0] = xp[0];
	_particle_list.getPos(p)[1] = xp[1];
	_particle_list.getPos(p)[2] = xp[2];

	_particle_list.getPos(q)[0] = xq[0];
	_particle_list.getPos(q)[1] = xq[1];
	_particle_list.getPos(q)[2] = xq[2];

	return true;

};



struct Collisions
{
	// Collison list
	//openfpm::vector<coll> coll_list;
	// Reactions
	Reactions *all_reactions;

	// Init
	Collisions(Reactions *_all_reactions) { all_reactions = _all_reactions; }

	//Deconstructor
	~Collisions(){}

	// Get the collison list
	template<typename CellList>
	inline void doCollisions(	td_particle_list  & _particle_list,
														CellList& _NN,
														Vcluster&  _v_cl,
														double _dt,
														std::mt19937 & _mt_rng,
														std::uniform_real_distribution<double> & _uniform_rnd,
														bool is_hs_collision = true,
													  bool is_reactive = true)
	{

		//Particle itterator for symmetric crossing scheme (HOW TO THIS WITH CELL LISTS ??? )
		//auto it2 = _particle_list.getParticleIteratorCRS(_NN.getInternalCellList());

		// Get a particle itterator
		auto it2 = _particle_list.getDomainIterator();


		// Private Collision lists
		std::list<coll> coll_list;


		// for each particle
		while (it2.isNext())
		{
			// Get p
			auto p = it2.get().getKey();

			// skip if the particle to chek is marked for deletion
			if ( _particle_list.getProp<id>(p) < 0) {++it2; continue;};


			// Get a symmetric neighborhood iterator

			//auto Np = _NN.getIterator(_NN.getCell(_particle_list.getPos(p)));
			auto Np = _NN.template getNNIteratorSym<NO_CHECK>(_NN.getCell(_particle_list.getPos(p)),p,_particle_list.getPosVector());

			// Get the position of p
			Point<3,double> xp = _particle_list.getPos(p);

			// get the radius of p
			double r_p = _particle_list.getProp<radius>(p);

			while (Np.isNext())
			{
				// ... q
				auto q = Np.get();
				// skip it self or if the neighbour is already marked for deletion
				//if (q == p ) {++Np; continue;};

				// Get the position of q
				Point<3,double> xq = _particle_list.getPos(q);

				// take the norm of this vector
				double rn = norm2(xp - xq);

				// get the radius of q
				double r_q = _particle_list.getProp<radius>(q);

				int this_id = _particle_list.getProp<id>(q);

				//std::cout << "Test for Collision" << std::endl;
				//std::cout << rn << std::endl;
				//double r_test = (r_q + r_p) * (r_q + r_p);
				//std::cout << ( r_test) << std::endl;

				// Test for collsions
				if ( (q == p) or this_id < 0 or  rn >= (r_q + r_p) * (r_q + r_p)) {++Np;continue;} // No Collision

				// If collsion add the collsion
				//std::cout << "Found a Collision" << std::endl;

				coll collsion = coll(q,p,rn);

				coll_list.push_back(collsion);

				// next neighbourgh
				++Np;

			}
			// Next Particle
			++it2;
		}

		// Sort the list
		//std::cout << "Sort" << std::endl;
		if (coll_list.size() == 0){return;}

		coll_list.sort();

		//std::cout << "Sorted" << std::endl;

		bool succes;

		// Init an id vector iterator
		//std::vector<coll>::iterator coll_it;
		std::list<coll>::iterator coll_it = coll_list.begin();



		size_t q_key = coll_it->key_q;
		size_t p_key = coll_it->key_p;
		double dist = coll_it->dist;


		// Allways do the first collision
		// Set collision flg for both particles
		_particle_list.getProp<collision_flg>(q_key) = true;
		_particle_list.getProp<collision_flg>(p_key) = true;

		//Log hypothetical bimolcular reactions
		bool is_reacted = false;
		if(is_reactive)
		{
			is_reacted = all_reactions->react_second_order(_particle_list, _NN, p_key, q_key, _dt ,_mt_rng, _uniform_rnd);
		}
		// Make HS collision
		if(is_hs_collision  && !is_reacted)
		{
			succes = hardSphereCollision(_particle_list,q_key,p_key,_dt);
		}

		//std::cout << "Coll distance (init)" << coll_it->dist << std::endl;

		// Check the rest of the collisions

		//std::cout << "Got the keys" << std::endl;
		// Itterate over all collision
		while(coll_it != coll_list.end())
		{

			// Collision will be seen with in the Domain as double > skip the collision
			// when it was done before!

			if(  coll_it->dist == dist or
			  ( _particle_list.getProp<id>(coll_it->key_q) < 0) or
			  ( _particle_list.getProp<id>(coll_it->key_p) < 0)  )

				{
					//std::cout << "Coll distance (skip) " << coll_it->dist << std::endl;
					//std::cout << "Coll distance (skip) " << q_key << std::endl;
					//std::cout << "Coll distance (skip) " << p_key << std::endl;
					++coll_it;
					continue;
				}


			//std::cout << "Coll distance (do)" << coll_it->dist << std::endl;

			//else do collisions

			// get distance and keys
			dist = coll_it->dist;
			q_key = coll_it->key_q;
			p_key = coll_it->key_p;

			//std::cout << "Execute collision" << std::endl;
			//std::cout << "p: "<< p_key << std::endl;
			//std::cout << "q: "<< q_key << std::endl;

			// Do a collision (reaction or HS)

			// Set collision flg for both particles
			_particle_list.getProp<collision_flg>(q_key) = true;
			_particle_list.getProp<collision_flg>(p_key) = true;

			//Log hypothetical bimolcular reactions
			bool is_reacted = all_reactions->react_second_order(_particle_list, _NN, p_key, q_key,_dt ,_mt_rng, _uniform_rnd);

			// Make HS collision
			if(is_hs_collision  && !is_reacted)
			{
				succes = hardSphereCollision(_particle_list,q_key,p_key,_dt);
			}

			++coll_it;
		}


	}// DoCollions


};// End struct

#endif //COLLISION_H
