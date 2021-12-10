#ifndef LATTICE_H
#define LATTICE_H



# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
#include <array>
#include <algorithm>

//declare global dimensions
const int Nx=50;
const int Ny=50;
const int Nxy = Nx*Ny;

int mod(int a, int b)
{
    int r=a%b;
    return r>=0 ? r: std::abs(b)+r;
}


struct Neighbours
{
    std::vector<int> positions;
    std::vector<int> slopes;
};


class Lattice
{
private:




public:

    inline std::array<int,2> get_coord(const int pos) const
    {
        //transform 1d array into 2d notation
        int x=pos%Nx;
        int y=pos/Nx;

        std::array<int, 2> coord{x,y};


        return coord;
    }

    inline int get_single_index(const int x, const int y) const
    {
        int index = y*Nx+x;
        return index;
    }

    inline bool is_in_lattice(const std::array<int,2> coord) const
    {
        int x=coord[0];
        int y= coord[1];
        return ( (unsigned int)x<Nx && (unsigned int)y<Ny );//negative values are looped around to big numbers by the cast to unsigned. In this way, only one inequality is required.
    }

    Neighbours get_neighbors(const int pos) const
    {
        Neighbours n;

        std::array<int,2> coord{get_coord(pos)};

        int x=coord[0];
        int y=coord[1];


        if(y%2==0) //if the position is in an even row, which protude on the left side
        {
            if(x-1>=0)
            {
                n.positions.emplace_back(get_single_index(x-1,y));
                n.slopes.emplace_back(0);
                n.positions.emplace_back(get_single_index(x-1,mod(y+1,Ny)));
                n.slopes.emplace_back(-1);
                n.positions.emplace_back(get_single_index(x-1,mod(y-1,Ny)));
                n.slopes.emplace_back(1);
            }

        if(x+1<Nx)
        {
            n.positions.emplace_back(get_single_index(x+1,y));
            n.slopes.emplace_back(0);
        }

        n.positions.emplace_back(get_single_index(x,mod(y+1,Ny)));
        n.slopes.emplace_back(1);

        n.positions.emplace_back(get_single_index(x,mod(y-1,Ny)));
        n.slopes.emplace_back(-1);
    }


    else //if the position is in an odd row, which protude on the right side
    {
        if(x-1>=0)
        {
            n.positions.emplace_back(get_single_index(x-1,y));
            n.slopes.emplace_back(0);
        }

        if(x+1<Nx)
        {
            n.positions.emplace_back(get_single_index(x+1,y));
            n.slopes.emplace_back(0);
            n.positions.emplace_back(get_single_index(x+1,mod(y+1,Ny)));
            n.slopes.emplace_back(1);
            n.positions.emplace_back(get_single_index(x+1,mod(y-1,Ny)));
            n.slopes.emplace_back(-1);
        }

        n.positions.emplace_back(get_single_index(x,mod(y+1,Ny)));
        n.slopes.emplace_back(-1);
        n.positions.emplace_back(get_single_index(x,mod(y-1,Ny)));
        n.slopes.emplace_back(1);



    }
    //sort(neighbors.begin(),neighbors.end());
    return n;
}



};

#endif // LATTICE_H
