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
const int Nx=4;
const int Ny=4;
const int Nxy = Nx*Ny;

int mod(int a, int b)
{
    int r;
    return r=a%b>=0 ? a%b: std::abs(b)+a%b;
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
    std::array<int,6> slope;



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

    Neighbours get_neighbors(int pos)
    {
        Neighbours n;

        std::array<int,2> coord{get_index(pos)};

        int x=coord[0];
        int y= coord[1];



        if(y%2==0)
        {
            if(x-1>=0)
            {
                n.positions.push_back(get_single_index(x-1,y));
                n.slopes.push_back(0);
                n.positions.push_back(get_single_index(x-1,mod(y+1,Ny)));
                n.slopes.push_back(-1);
                n.positions.push_back(get_single_index(x-1,mod(y-1,Ny)));
                n.slopes.push_back(1);
            }

        if(x+1<Nx)
        {
            n.positions.push_back(get_single_index(x+1,y));
            n.slopes.push_back(0);
        }

        n.positions.push_back(get_single_index(x,mod(y+1,Ny)));
        n.slopes.push_back(1);

        n.positions.push_back(get_single_index(x,mod(y-1,Ny)));
        n.slopes.push_back(-1);
    }


    else
    {
        if(x-1>=0)
        {
            n.positions.push_back(get_single_index(x-1,y));
            n.slopes.push_back(0);
        }

        if(x+1<Nx)
        {
            n.positions.push_back(get_single_index(x+1,y));
            n.slopes.push_back(0);
            n.positions.push_back(get_single_index(x+1,mod(y+1,Ny)));
            n.slopes.push_back(1);
            n.positions.push_back(get_single_index(x+1,mod(y-1,Ny)));
            n.slopes.push_back(-1);
        }

        n.positions.push_back(get_single_index(x,mod(y+1,Ny)));
        n.slopes.push_back(-1);
        n.positions.push_back(get_single_index(x,mod(y-1,Ny)));
        n.slopes.push_back(1);



    }
    //sort(neighbors.begin(),neighbors.end());
    return n;
}



};

#endif // LATTICE_H
