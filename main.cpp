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

//random number generator

unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);

inline static float fast_exp (float x)
{
    //if(x<-1000) return 0.0;
    volatile union
    {
        float f;
        unsigned int i;
    } cvt;

    /* exp(x) = 2^i * 2^f; i = floor (log2(e) * x), 0 <= f <= 1 */
    float t = x * 1.442695041f;
    float fi = floorf (t);
    float f = t - fi;
    int i = (int)fi;
    cvt.f = (0.3371894346f * f + 0.657636276f) * f + 1.00172476f; /* compute 2^f */
    cvt.i += (i << 23);                                          /* scale by 2^i */
    return cvt.f;
}

// find empty hexes (without neighbors?)
// randomely create a particle with k_on/
//RETURN updated grid and diffuse position

//Do I need to check how full the grid is? Set concentration?delete particles randomely?
int create_particle(Particles &particles)
{
    double const k_on = 0.5;
    std::vector<int> empty_hex(Nxy);
    for(int x = 0; x < Nxy; ++x)
    {
        empty_hex[x] = x;
    }

    std::vector<int> occupied_hex;
    occupied_hex.insert(occupied_hex.end(),particles.bound_pos.begin(),particles.bound_pos.end());
    occupied_hex.insert(occupied_hex.end(),particles.diffuse_pos.begin(),particles.diffuse_pos.end());

    for(auto it=occupied_hex.begin(); it!=occupied_hex.end(); ++it)
    {
        empty_hex.erase(std::remove(empty_hex.begin(),empty_hex.end(),*it),empty_hex.end());
    }
    double rand = unidist(gen)*(empty_hex.size());
    std::cout<<"creation rand"<<rand<<"\n";
    int new_particle_pos=empty_hex[std::floor(rand)];
    std::cout<<"creation site"<<new_particle_pos<<"\n";


    if((rand-std::floor(rand))<k_on)// create particle
    {
        particles.grid[new_particle_pos]=1;
        particles.diffuse_pos.push_back(new_particle_pos);
        return new_particle_pos;

    }
    else
    {
        return Nxy+1;
    }
}





//calculate energy change by binding/unbinding a particle
//


//somehow a negative 3 appears






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
