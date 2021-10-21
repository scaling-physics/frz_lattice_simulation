# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <array>

const int L=10;
class lattice
{
private:




public:
   std::array<int,pow(L,2)> grid; //0 if empty, 1 if occupied but diffuse; orientation if occupied and bound
   std::vector<int> diffuse_pos; //reference to position of diffuse particles
   std::vector<int> bound_pos; //reference to position of bound particles
   std::vector<int> neighbors[L];


   lattice():neighbors{}{neighbors[0]=1;
   neighbors[1]=-1;
   neighbors[2]=L+1;
   neighbors[3]=L-1;
   neighbors[4]=L;
   neighbors[5]=-L;
   };
};
