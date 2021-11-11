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
#include "definitions.h"
//using namespace std;

//random number generator

unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);

int main()
{
    const int MC_steps = 1; // number of Monte Carlo Steps
    int MC_counter = 0;
    double rand;

//constants for reaction:
    double const alpha=1;
    double const J=1;

// Input and Output arrays

    int site;
    int new_particle_site;

    particles.grid[13]=1;
//    particles.grid[2]=1;
//    particles.grid[3]=1;
    particles.grid[9]=3;
    particles.grid[5]=3;
    particles.grid[10]=3;

    particles.get_diffuse_pos();
    particles.get_bound_pos();


    while(MC_counter<MC_steps)
    {
        if(MC_steps%100==0)
        {
        // select a random hex on grid:
        double rand_size = unidist(gen) * Nxy;
        int pos = rand_size;
        double rand = rand_size-pos;
        //CREATION ATTEMPT
        if(particles.is_free(pos))
        {

        }

        //DESTRUCTION ATTEMPT
        if(particles.is_diffuse(particles.get_index(pos)))
        {

        }


        }
        std::cout <<"\n counter"<< MC_counter <<'\n';

        //choose random particle
        rand = unidist(gen)*particles.positions.size();
        std::cout<<"rand"<<rand<<"\n";
        int ind=rand;
        rand=rand-ind;

        if(particles.is_diffuse(ind))
        {



//Move diffusive particles

        site=particles.diffuse(rand, ind);
        std::cout<<"site"<<site<<"\n";

//BINDING
        particles.binding_attempt(alpha,J,site,rand);


        }
        if(particles.is_bound(ind))
        {
// UNBINDING ATTEMPT
        }

//        new_particle_site=create_particle(particles);
//        std::cout<<"new particle created"<<" "<<new_particle_site<<"\n";







//UNBINDING
//perform unbinding_attempt for a random bound particle


        MC_counter++;

    }
    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
    {
        std::cout << *iter << "\n"<<' ';
    }
    for (auto iter = particles.diffuse_pos.begin(); iter !=particles.diffuse_pos.end(); ++iter)
    {
        std::cout << *iter << ","<<' ';
    }
    return 0;
};
