 /*
  * aggregate_copy.hpp
  *
  *  Created on:
  *      Author: Daniel
 */
#include <cstdlib>


#ifndef ADD_OP_H
#define ADD_OP_H

template<typename Tdst, typename Tsrc>
struct delete_parent_
{
	static inline void operation(Tdst & dst, const Tsrc & src)
    {


      //if (src < 0 and dst < 0)
      //    throw std::runtime_error("Ghost and particle collide that should never happen!");

        if (src < 0 and dst >= 0)
            dst = INT_MIN/3;

    }

 };

template<typename Tdst, typename Tsrc>
struct overwrite_bool_
{
	static inline void operation(Tdst & dst, const Tsrc & src)
    {

        if (src == true and dst == false)
            dst = true;

    }

 };


#endif //ADD_OP_H
