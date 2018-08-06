

%module pyhsbrd
%include "std_string.i"
%include "std_map.i"
%include "std_vector.i"


namespace std {
%template(vector_double) vector<double>;
%template(map_string_int) map<string ,int >;
%template(map_string_double) map<string, double> ;
%template(map_string_vector_double) map<string, vector<double> > ;
%template(map_string_map_string_double) map< string,  map<string, double> > ;
%template(map_string_map_string_vector_double) map<string,map<string,vector<double> > >;
%template(vector_uint) std::vector<size_t> ;
%template(map_int_double) map<int, double> ;
%template(map_string_map_int_double) map<string, map<int, double> > ;
%template(vector_map_string_map_int_double) vector< map<string, map<int, double> > >;
}

%{
    #include <string>
    #include <map>
    #include <cstdlib>
    #include "pyhsbrd.hpp"
%}

%include "pyhsbrd.hpp"
%include "pyhsbrd.cpp"
