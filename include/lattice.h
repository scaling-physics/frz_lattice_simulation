#ifndef LATTICE_H
#define LATTICE_H



# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
# include <sstream>      // std::stringstrea
# include <vector>
# include <array>
# include <algorithm>

//declare global dimensions
int Nx=70;
int Ny=30;
int Nxy = Nx*Ny;

inline int mod(int a, int b)
{
    int r=a%b;
    return r>=0 ? r: std::abs(b)+r;
}


//struct Neighbours
//{
//    std::vector<int> positions;
//    std::vector<int> slopes;
//};

struct Neighbour
{
    int position;
    int slope;
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
//        int x=pos/Ny;
//        int y=pos%Ny;

        std::array<int, 2> coord{x,y};
        return coord;
    }

    inline int get_single_index(const int x, const int y) const
    {
        int index = y*Nx+x;
//        int index = x*Ny+y;
        return index;
    }

//    inline bool is_in_lattice(const std::array<int,2> coord) const
//    {
//        int x=coord[0];
//        int y= coord[1];
//        return ( (unsigned int)x<Nx && (unsigned int)y<Ny );//negative values are looped around to big numbers by the cast to unsigned. In this way, only one inequality is required.
//    }


    std::vector<Neighbour> get_neighbors2(const int pos) const
    {
        std::vector<Neighbour> n;

        std::array<int,2> coord{get_coord(pos)};

        int x=coord[0];
        int y=coord[1];
////////// if counting rowwise
////////        if(unsigned(x-1)<Nx-2)
////////        {
////////            if(y%2==0)
////////            {
////////                n= {{get_single_index(x-1,y),4},{get_single_index(x-1,mod(y+1,Ny)),3},{get_single_index(x-1,mod(y-1,Ny)),5},{get_single_index(x+1,y),1},{get_single_index(x,mod(y+1,Ny)),2},{get_single_index(x,mod(y-1,Ny)),0}};
////////            }
////////            else
////////            {
////////                n= {{get_single_index(x-1,y),4},{get_single_index(x,mod(y+1,Ny)),3},{get_single_index(x,mod(y-1,Ny)),5},{get_single_index(x+1,y),1},{get_single_index(x+1,mod(y+1,Ny)),2},{get_single_index(x+1,mod(y-1,Ny)),0}};
////////            }
////////        }
////////        else if(x==0)
////////        {
////////            if(y%2==0)
////////            {
////////                n= {{get_single_index(x+1,y),1},{get_single_index(x,mod(y+1,Ny)),2},{get_single_index(x,mod(y-1,Ny)),0}};
////////            }
////////            else
////////            {
////////                n= {{get_single_index(x,mod(y+1,Ny)),3},{get_single_index(x,mod(y-1,Ny)),5},{get_single_index(x+1,y),1},{get_single_index(x+1,mod(y+1,Ny)),2},{get_single_index(x+1,mod(y-1,Ny)),0}};
////////            }
////////        }
////////        else
////////        {
////////            if(y%2==0)
////////            {
////////                n= {{get_single_index(x-1,y),4},{get_single_index(x-1,mod(y+1,Ny)),3},{get_single_index(x-1,mod(y-1,Ny)),5},{get_single_index(x,mod(y+1,Ny)),2},{get_single_index(x,mod(y-1,Ny)),0}};
////////            }
////////            else
////////            {
////////                n= {{get_single_index(x-1,y),4},{get_single_index(x,mod(y+1,Ny)),3},{get_single_index(x,mod(y-1,Ny)),5}};
////////            }
////////        }
        if(unsigned(x-1)<Nx-2)
        {
            if(y%2==0)
            {
                n= {{get_single_index(x-1,y),4},{get_single_index(x-1,mod(y+1,Ny)),3},{get_single_index(x-1,mod(y-1,Ny)),5},{get_single_index(x+1,y),1},{get_single_index(x,mod(y+1,Ny)),2},{get_single_index(x,mod(y-1,Ny)),0}};
            }
            else
            {
                n= {{get_single_index(x-1,y),4},{get_single_index(x,mod(y+1,Ny)),3},{get_single_index(x,mod(y-1,Ny)),5},{get_single_index(x+1,y),1},{get_single_index(x+1,mod(y+1,Ny)),2},{get_single_index(x+1,mod(y-1,Ny)),0}};
            }
        }
        else if(x==0)
        {
            if(y%2==0)
            {
                n= {{get_single_index(x+1,y),1},{get_single_index(x,mod(y+1,Ny)),2},{get_single_index(x,mod(y-1,Ny)),0}};
            }
            else
            {
                n= {{get_single_index(x,mod(y+1,Ny)),3},{get_single_index(x,mod(y-1,Ny)),5},{get_single_index(x+1,y),1},{get_single_index(x+1,mod(y+1,Ny)),2},{get_single_index(x+1,mod(y-1,Ny)),0}};
            }
        }
        else
        {
            if(y%2==0)
            {
                n= {{get_single_index(x-1,y),4},{get_single_index(x-1,mod(y+1,Ny)),3},{get_single_index(x-1,mod(y-1,Ny)),5},{get_single_index(x,mod(y+1,Ny)),2},{get_single_index(x,mod(y-1,Ny)),0}};
            }
            else
            {
                n= {{get_single_index(x-1,y),4},{get_single_index(x,mod(y+1,Ny)),3},{get_single_index(x,mod(y-1,Ny)),5}};
            }
        }

//
//        if(unsigned(x-1)<Nx-2)
//        {
//            if(y%2==0)
//            {
//                n= {{get_single_index(x-1,y),0},{get_single_index(x-1,mod(y+1,Ny)),-1},{get_single_index(x-1,mod(y-1,Ny)),1},{get_single_index(x+1,y),0},{get_single_index(x,mod(y+1,Ny)),1},{get_single_index(x,mod(y-1,Ny)),-1}};
//            }
//            else
//            {
//                n= {{get_single_index(x-1,y),0},{get_single_index(x,mod(y+1,Ny)),-1},{get_single_index(x,mod(y-1,Ny)),1},{get_single_index(x+1,y),0},{get_single_index(x+1,mod(y+1,Ny)),1},{get_single_index(x+1,mod(y-1,Ny)),-1}};
//            }
//        }
//        else if(x==0)
//        {
//            if(y%2==0)
//            {
//                n= {{get_single_index(x+1,y),0},{get_single_index(x,mod(y+1,Ny)),1},{get_single_index(x,mod(y-1,Ny)),-1}};
//            }
//            else
//            {
//                n= {{get_single_index(x,mod(y+1,Ny)),-1},{get_single_index(x,mod(y-1,Ny)),1},{get_single_index(x+1,y),0},{get_single_index(x+1,mod(y+1,Ny)),1},{get_single_index(x+1,mod(y-1,Ny)),-1}};
//            }
//        }
//        else
//        {
//            if(y%2==0)
//            {
//                n= {{get_single_index(x-1,y),0},{get_single_index(x-1,mod(y+1,Ny)),-1},{get_single_index(x-1,mod(y-1,Ny)),1},{get_single_index(x,mod(y+1,Ny)),1},{get_single_index(x,mod(y-1,Ny)),-1}};
//            }
//            else
//            {
//                n= {{get_single_index(x-1,y),0},{get_single_index(x,mod(y+1,Ny)),-1},{get_single_index(x,mod(y-1,Ny)),1}};
//            }
//        }

        return n;
    }



};

#endif // LATTICE_H
