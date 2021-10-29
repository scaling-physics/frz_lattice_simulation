# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <array>
#include <vector>
# include <chrono>
# include <random>
#include <algorithm>

#include "lattice.h"
Lattice lattice;
#include "Particles.h"
Particles particles;
//using namespace std;

//generate random number


unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);



// fct takes position of diffuse particles attempts to move one at random,
//if move goes to a lattice point without neighbours then move is accepted
//if move attempts to move next to a particle --> return particle combination for step 3
// RETURN, updated diffuse_pos and grid, create array with attempted binding

int diffuse(std::vector<int> &diffuse_pos,std::array<short,Nxy> &grid)
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

    else
    {
        return particle_pos;
    }

}

// find empty hexes (without neighbors)
// randomely create a particle with k_on/
//RETURN updated grid and diffuse position

//Do I need to check how full the grid is? Set concentration?delete particles randomely?
int create_particle(std::vector<int> &diffuse_pos,std::vector<int> &bound_pos,std::array<short,Nxy> &grid)
{
    std::vector<int> empty_hex(Nxy);
    for(int x = 0; x < Nxy; ++x)
    {
        empty_hex[x] = x;
    }

    std::vector<int> occupied_hex;
    occupied_hex.insert(occupied_hex.end(),bound_pos.begin(),bound_pos.end());
    occupied_hex.insert(occupied_hex.end(),diffuse_pos.begin(),diffuse_pos.end());

    for(auto it=occupied_hex.begin(); it!=occupied_hex.end(); ++it)
    {
        empty_hex.erase(std::remove(empty_hex.begin(),empty_hex.end(),*it),empty_hex.end());
    }
    float rand = unidist(gen)*(empty_hex.size());
    int new_particle_pos=std::floor(rand);
    int const k_on = 0.5;

    if((rand-std::floor(rand))<k_on)// accept move
    {
        grid[new_particle_pos]=1;
        diffuse_pos.push_back(new_particle_pos);
        return new_particle_pos;

    }
    else {return Nxy+1}
}




//calculate energy change by binding/unbinding a particle
//
void binding_attempt(std::vector<int> &bound_pos,std::array<short,Nxy> &grid, int const &alpha, int const &J)
{
//    int no_bound=bound_pos.size();




}








int main()
{
    const int MC_steps = 5; // number of Monte Carlo Steps
    int MC_counter = 0;

//constants for reaction:
//    int const alpha=1;
//    int const J=1;

// Input and Output arrays

    int site;
    int new_particle_site;

    particles.grid[1]=1;
    particles.grid[2]=1;
    particles.grid[3]=1;
    particles.grid[4]=3;
    particles.grid[5]=2;
    particles.grid[6]=3;


    particles.get_diffuse_pos();
    particles.get_bound_pos();


    while(MC_counter<MC_steps)
    {

        std::cout << MC_counter <<'\n';
//step 1: Move diffusive particles

        site=diffuse(particles.diffuse_pos,particles.grid);
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
    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
    {
        std::cout << *iter << "\n"<<' ';
    }
    return 0;
};
