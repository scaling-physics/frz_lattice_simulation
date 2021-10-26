#ifndef lattice_H
#define lattice_H



# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <array>

const int L=10;
int Lsq = L*L
class lattice
{
private:




public:
   std::array<short,Lsq> grid{0}; //0 if empty, 1 if occupied but diffuse; orientation if occupied and bound
   std::vector<int> diffuse_pos; //reference to position of diffuse particles
   std::vector<int> bound_pos; //reference to position of bound particles


  std::vector<int> get_neighbors(int const &pos, int const &L)
        {std::vector<int> neighbors;

        //corners
        std::array<int,4> corners={0,L-1, (L-1)*L,L*L-1};
        //edges
        std::vector<int> edge_top(L-2);
            for (int i=0;i<L-2;i++){edge_top[i]=i+1;}
        std::vector<int> edge_bottom(L-2);
            for (int i=0;i<L-2;i++){edge_bottom[i]=(L-1)*L+i+1;}
        std::vector<int> edge_left(L-2);
            for (int i=0;i<L-2;i++){edge_left[i]=(i+1)*L;}
        std::vector<int> edge_right(L-2);
            for (int i=0;i<L-2;i++){edge_right[i]=(i+2)*L-1;}

        for (int j= 0; j<2; j++)
        {cout<<edge_left[j] << '\n';}
        for (int j= 0; j<2; j++)
        {cout<<edge_right[j] << '\n';}


        if(pos==corners[0]){
            neighbors={pos+1,pos+L,pos+((L-1)*L),pos+((L-1)*L+1)};
            cout <<"corner 0 being executed";
            }

        else if(pos==corners[1]){
            neighbors={pos-1,pos+(L-1),pos+L,pos+((L-1)*L)};
            cout <<"corner 1 being executed";
            }
        else if(pos==corners[2]){
            neighbors={pos-(L-1)*L,pos-L,pos-(L+1),pos+1};
            cout <<"corner 2 being executed";
            }
        else if(pos==corners[3]){
            neighbors={pos-(L-1)*L-1,pos-(L-1)*L,pos-L,pos-1};
            cout <<"corner 3 being executed";
            }

        else if(std::find(edge_top.begin(), edge_top.end(),pos) !=edge_top.end()){
            if(pos%2==0)
            {neighbors={pos-1, pos+1, pos+L, pos+(L-1)*L-1, pos+(L-1)*L, pos+(L-1)*L+1};}
            else
            {neighbors={pos-1, pos+1, pos+L-1, pos+L, pos+L+1, pos+(L-1)*L};}
            }

        else if(std::find(edge_bottom.begin(), edge_bottom.end(),pos) !=edge_bottom.end()){
            if(pos%2==0)
            {neighbors={pos-(L-1)*L, pos-L-1, pos-L, pos-L+1, pos-1, pos+1};}
            else
            {neighbors={pos-(L-1)*L-1, pos-(L-1)*L, pos-(L-1)*L+1, pos-L, pos-1, pos+1};}
            }

        else if(std::find(edge_left.begin(), edge_left.end(),pos) !=edge_left.end()){
            neighbors={pos-L, pos-L+1, pos+1, pos+L};
            }

        else if(std::find(edge_right.begin(), edge_right.end(),pos) !=edge_right.end()){
            neighbors={pos-L, pos-1, pos+L-1,pos+L};
            }

        else{
            if(pos%2==0)
            {neighbors={pos-L-1, pos-L, pos-L+1, pos-1, pos+1, pos+L};}
            else
            {neighbors={pos-L, pos-1, pos+1, pos+L-1, pos+L, pos+L+1};}
            }
            return neighbors;
        }
};

#endif // lattice_H
