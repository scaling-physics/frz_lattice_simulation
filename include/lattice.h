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
        std::set<int> edge_top;
        std::set<int> edge_bottom;
        std::set<int> edge_left;
        std::set<int> edge_right;

        for (int j= 0; j<4; j++)
        {cout<<corners[j] << '\n';}


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

        else if(edge_top.find(pos) !=edge_top.end()){
        neighbors={0,0,0,0,0,0};
        }

        else if(edge_bottom.find(pos) !=edge_bottom.end()){
        neighbors={0,0,0,0,0,0};
        }

        else if(edge_left.find(pos) !=edge_left.end()){
        neighbors={0,0,0,0,0,0};
        }

        else if(edge_right.find(pos) !=edge_right.end()){
        neighbors={0,0,0,0,0,0};
        }

        else{
        neighbors={0,0,0,0,0,0};
        }




return neighbors;

}

};
