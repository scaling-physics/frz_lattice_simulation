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

void init(Particles &particles)
{
    particles.grid[5]=1;
    particles.grid[9]=2;
    particles.grid[10]=3;

    particles.positions.push_back(5);
    particles.positions.push_back(9);
    particles.positions.push_back(10);

    particles.orientations.push_back(2);
    particles.orientations.push_back(2);
    particles.orientations.push_back(2);
}
void print_container(const std::vector<int>& c)
{
    for (auto &i : c)
    {
        std::cout << i << " ";
    }
    std::cout << '\n';
}

int main()
{
    const int MC_steps = 3; // number of Monte Carlo Steps
    int MC_counter = 0;
    double rand;

//constants for reaction:
    double const alpha=1;
    double const J=1;

// Input and Output arrays

    int site;
    init(particles);

    while(MC_counter<MC_steps)
    {
        std::cout <<"counter"<< MC_counter <<'\n';
        if(MC_steps%1==0)
        {
            // select a random hex on grid:
            double rand_size = unidist(gen) * Nxy;
            int pos = rand_size;
            double rand = rand_size-pos;
            std::cout<<"pos "<<pos<<"\n";
            std::cout<<"rand "<<rand<<"\n";




            //CREATION ATTEMPT
            if(particles.is_free(pos))
            {
                std::cout<<"free"<<"\n";
                particles.creation_attempt(pos,rand);
                print_container(particles.positions);
            }

        }

//        if(MC_steps%1==0)
//        {
//        // select a random hex on grid:
//        double rand_size = unidist(gen) * Nxy;
//        int pos = rand_size;
//        double rand = rand_size-pos;
//        std::cout<<"pos "<<pos<<"\n";
//        std::cout<<"rand "<<rand<<"\n";
//        //DESTRUCTION ATTEMPT
//        if(particles.is_diffuse(particles.get_ind(pos)))
//        {
//            std::cout<<"diffuse "<<"\n";
//            particles.destruction_attempt(pos,rand);
//
//        }
//
//        }
        //choose random particle
//        rand = unidist(gen)*particles.positions.size();
//        std::cout<<"rand"<<rand<<"\n";
//        int ind=rand;
//        rand=rand-ind;
//
//        if(particles.is_diffuse(ind))
//        {
//
//
//
////Move diffusive particles
//
//        site=particles.diffuse(rand, ind);
//        std::cout<<"site"<<site<<"\n";
//
////BINDING
//        particles.binding_attempt(alpha,J,site,rand);
//
//
//        }
//        if(particles.is_bound(ind))
//        {
//// UNBINDING ATTEMPT
//        }
//
////        new_particle_site=create_particle(particles);
////        std::cout<<"new particle created"<<" "<<new_particle_site<<"\n";
//
//
//
//
//
//
//
////UNBINDING
////perform unbinding_attempt for a random bound particle


        MC_counter++;

    }
    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
    {
        std::cout << *iter << "\n"<<' ';
    }

    return 0;
};
