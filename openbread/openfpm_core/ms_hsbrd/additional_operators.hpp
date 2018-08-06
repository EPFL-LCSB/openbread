 /*
  * aggregate_copy.hpp
  *
  *  Created on:
  *      Author: Daniel
 */

 #ifndef ADD_OP_H
 #define ADD_OP_H

template<typename Tdst, typename Tsrc>
struct replace_if_collsion_
{
	static inline void operation(Tdst & dst, const Tsrc & src)
    {

    	auto coll_src = src. template get<collision_flg>();
    	auto id_src = src. template get<id>();

    	auto coll_dst = dst. template get<collision_flg>();
    	auto id_dst = dst. template get<id>();

    	//std::cerr << "the check variable: " <<stuff << std::endl;
    	// If ghost
    	if( id_src < 0 )
    	{
    		dst = src;
    	}
    	// If collision and not deleted
    	if(coll_src && id_dst>= 0 )
    	{
    		dst = src;
    	}

        // Reverse if
        if(id_dst< 0 )
        {
            src = dst;
        }

        // If collision and not deleted
        if(coll_dst && id_src>= 0 )
        {
            src = dst;
        }

    }
 };

 #endif //ADD_OP_H
