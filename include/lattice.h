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
const int L=4;
const int Nxy = Nx*Ny;

int mod(int a, int b)
{
    int r;
    return r=a%b>=0 ? a%b: std::abs(b)+a%b;
}

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

    std::vector<int> get_neighbors(int pos)
    {
        int j=0;
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
                slope[j]=0;
                neighbors.push_back(get_single_index(x-1,mod(y+1,Ny)));
                slope[j+1]=-1;
                neighbors.push_back(get_single_index(x-1,mod(y-1,Ny)));
                slope[j+2]=1;
                j=3;
            }

        if(x+1<Nx)
        {
            neighbors.push_back(get_single_index(x+1,y));
            slope[j]=0;
            j++;
        }

        neighbors.push_back(get_single_index(x,mod(y+1,Ny)));
        slope[j]=1;
        j++;
        neighbors.push_back(get_single_index(x,mod(y-1,Ny)));
        slope[j]=-1;

    }


    else
    {
        if(x-1>=0)
        {
            neighbors.push_back(get_single_index(x-1,y));
            slope[j]=0;
            j++;
        }

        if(x+1<Nx)
        {
            neighbors.push_back(get_single_index(x+1,y));
            slope[j]=0;
            j++;
            neighbors.push_back(get_single_index(x+1,mod(y+1,Ny)));
            slope[j]=1;
            j++;
            neighbors.push_back(get_single_index(x+1,mod(y-1,Ny)));
            slope[j]=-1;
            j++;
        }

        neighbors.push_back(get_single_index(x,mod(y+1,Ny)));
        slope[j]=-1;
        j++;
        neighbors.push_back(get_single_index(x,mod(y-1,Ny)));
        slope[j]=1;



    }
    //sort(neighbors.begin(),neighbors.end());
    return neighbors;
}



};

#endif // LATTICE_H
