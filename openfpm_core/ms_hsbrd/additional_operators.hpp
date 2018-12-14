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

        if (src < 0)
            dst = INT_MIN/3;

    }
 };

#endif //ADD_OP_H
