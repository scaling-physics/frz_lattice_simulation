#ifndef lattice_H
#define lattice_H



# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
#include <array>
#include <algorithm>
#include <set>
# include <chrono>
# include <random>


class Lattice
{
private:




public:
    const int Nx=4;
    const int Ny=4;
    const int L=4;
    const int Lsq = Nx*Ny;


    std::array<short,16> grid{0}; //0 if empty, 1 if occupied but diffuse; orientation if occupied and bound
    std::vector<int> diffuse_pos={1,2,3}; //reference to position of diffuse particles
    std::vector<int> bound_pos; //reference to position of bound particles


    std::array<int,2> get_index(int const pos)
    {
        //transform 1d array into 2d notation
        int x=pos%Nx;
        int y=pos/Nx;

        std::array<int, 2> coord{x,y};


        return coord;
    }

    int get_single_index(int x, int y)
    {
        int index = y*Ny+x;
        return index;

    }

    bool is_in_lattice(std::array<int,2> coord)
    {
        int x=coord[0];
        int y= coord[1];
        return ( (unsigned int)x<Nx && (unsigned int)y<Ny );//negative values are looped around to big numbers by the cast to unsigned. In this way, only one inequality is required.
    }

    std::vector<int> get_neighbors(int pos)
    {
        std::vector<int> neighbors;
        std::array<int,2> coord;
        coord = get_index(pos);

        int x=coord[0];
        int y= coord[1];



        if(y%2==0)
        {
            if(x-1>=0)
            {
                neighbors.push_back(get_single_index(x-1,y));
                neighbors.push_back(get_single_index(x-1,(y+1)%Ny));
                neighbors.push_back(get_single_index(x-1,(Ny+y-1)%Ny));
            }

        if(x+1<Nx)
        {
            neighbors.push_back(get_single_index(x+1,y));
        }

        neighbors.push_back(get_single_index(x,(y+1)%Ny));
        neighbors.push_back(get_single_index(x,(Ny+y-1)%Ny));

    }


    else
    {
        if(x-1>=0)
        {
            neighbors.push_back(get_single_index(x-1,y));
        }

        if(x+1<Nx)
        {
            neighbors.push_back(get_single_index(x+1,y));
            neighbors.push_back(get_single_index(x+1,(y+1)%Ny));
            neighbors.push_back(get_single_index(x+1,(Ny+y-1)%Ny));
        }

        neighbors.push_back(get_single_index(x,(y+1)%Ny));
        neighbors.push_back(get_single_index(x,(Ny+y-1)%Ny));



    }
    sort(neighbors.begin(),neighbors.end());
    return neighbors;
}


};

#endif // lattice_H
