# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <array>
#include <vector>
#include "lattice.h"
# include <chrono>
# include <random>
#include <algorithm>
Lattice lattice;
//using namespace std;

//generate random number


unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);



// fct takes position of diffuse particles attempts to move one at random,
//if move goes to a lattice point without neighbours then move is accepted
//if move attempts to move next to a particle --> return particle combination for step 3
// RETURN, updated diffuse_pos and grid, create array with attempted binding

int diffuse(std::vector<int> &diffuse_pos,std::array<short,16> &grid)
{
//randomely choose diffuse_pos


    float rand = unidist(gen)*diffuse_pos.size();
    std::cout<<"rand"<<rand<<"\n";
    int particle_pos = diffuse_pos[std::floor(rand)] ;
    std::cout<<"particle"<<particle_pos<<"\n";
//choose random direction
    std::vector<int> neighbors_dif;
    neighbors_dif=lattice.get_neighbors(particle_pos);


    int new_pos =neighbors_dif[ std::floor(neighbors_dif.size()*(rand - std::floor(rand)))];
    std::cout<<"new_pos"<<new_pos<<"\n";
//check if new hex is occupied
    if(grid[new_pos]==0)// accept move
    {
        grid[particle_pos]=0;
        grid[new_pos]=1;
        diffuse_pos[std::floor(rand)]=new_pos;
        //NEED TO UPDATE DIFFUSE_pos!!!!!!
        return new_pos;

    }

    else{
    return particle_pos;}

}

// find empty hexes without neighbours
// randomely create a particle with k_on/
//RETURN updated grid and diffuse position

//Do I need to check how full the grid is? Set concentration?delete particles randomely?
//int create_particle(vector<int> &diffuse_pos,array<short,Lsq> &grid)
//{
//int const k_on = 0.02;
//
//int new_particle_pos;
//
//return new_particle pos;
//}




//calculate energy change by binding/unbinding a particle
//
void binding_attempt(std::vector<int> &bound_pos,std::array<short,16> &grid, int const &alpha, int const &J)
{
    int no_bound=bound_pos.size();




}








int main()
{
    const int MC_steps = 5; // number of Monte Carlo Steps
    int MC_counter = 0;

//constants for reaction:
    int const alpha=1;
    int const J=1;

// Input and Output arrays


    int site;
    int new_particle_site;





    while(MC_counter<MC_steps)
    {

        std::cout << MC_counter <<'\n';
//step 1: Move diffusive particles

        site=diffuse(lattice.diffuse_pos,lattice.grid);
        std::cout<<"site"<<site<<"\n";
//if move is rejected check if it binds or not




//step 2: Check and create a new particles

//new_particle_site=create_particle(diffuse_pos,grid);

//step 3: bind and unbind
//BINDING
//first take site and check if particle moved adjacent to particle and perform binding_attempt

//second take new_particle_site: if zero skip, if non-zero check if adjacent to particle and perform binding_attempt


//UNBINDING
//perform unbinding_attempt for a random bound particle


        MC_counter++;
    }
    for (auto iter = lattice.grid.begin(); iter !=lattice.grid.end(); ++iter) {
        std::cout << *iter << "\n"<<' ';
        }
    return 0;
};
