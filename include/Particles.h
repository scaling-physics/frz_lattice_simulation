#ifndef PARTICLES_H
#define PARTICLES_H

# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
#include <array>
#include <algorithm>
#include "lattice.h"

class Particles

{
private:

public:
    std::array<short,Nxy> grid{0}; //0 if empty, 1 if occupied but diffuse; orientation if occupied and bound
    std::vector<int> diffuse_pos; //reference to position of diffuse particles
    std::vector<int> bound_pos; //reference to position of bound particles

    void get_diffuse_pos()
    {
        for (int j=0; j<Nxy; j++)
        {
            if (grid[j]==1)
            {
                if(!(std::find(diffuse_pos.begin(), diffuse_pos.end(), j) != diffuse_pos.end()))
                {
                    diffuse_pos.push_back(j);
                }
            }
        }
        //sort(diffuse_pos.begin(),diffuse_pos.end());
    }

    void get_bound_pos()
    {
        for (int j=0; j<Nxy; j++)
        {
            if (grid[j]>1)
            {
                if(!(std::find(bound_pos.begin(), bound_pos.end(), j) != bound_pos.end()))
                {
                    bound_pos.push_back(j);
                }
            }
        }
        sort(bound_pos.begin(),bound_pos.end());
    }
};


#endif // PARTICLES_H
