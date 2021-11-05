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




// fct takes position of diffuse particle: attempts to move one at random,
//if move goes to a empty lattice point is accepted
// updated diffuse_pos and grid
//RETURN (new) particle position


int diffuse(Particles &particles)
{
//randomely choose diffuse_pos


    double  rand= unidist(gen)*particles.diffuse_pos.size();
    int ind=rand;
    std::cout<<"rand"<<rand<<"\n";
    int particle_pos = particles.diffuse_pos[ind] ;
    std::cout<<"particle"<<particle_pos<<"\n";
//choose random direction
    Neighbours neighbors_dif=lattice.get_neighbors(particle_pos);

    double rand2=(rand-ind)*neighbors_dif.positions.size();
    int dir=rand2;
    int new_pos =neighbors_dif.positions[dir];
    std::cout<<"new_pos"<<new_pos<<"\n";
//check if new hex is occupied
    if(particles.grid[new_pos]==0)// accept move
    {
        particles.grid[particle_pos]=0;
        particles.grid[new_pos]=1;
        particles.diffuse_pos[ind]=new_pos;
        //NEED TO UPDATE DIFFUSE_pos!!!!!!
        return new_pos;

    }
    else
    {
        return particle_pos;
    }

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

double delta_H(int const alpha, int const J, int no_dif, std::vector<int> orientations)
{
    double delta_E=J*(-orientations.size()+(no_dif*alpha));
    return delta_E<0.0f ? 1.0 : exp(-delta_E);
}



//calculate energy change by binding/unbinding a particle
//
void binding_attempt(Particles &particles, int const alpha, int const J, int pos)
{
    double rand=unidist(gen)*3-1;
    int ori = rand;
    //std::cout<<"creation rand"<<rand<<"\n";

    std::vector<int> orientations;
    std::vector<int> slopes;
    std::vector<int> possible_binding_pos;

    int no_dif = 1;
    int i=0;

    Neighbours neighbors=lattice.get_neighbors(pos);
    int orientation_a=ori;


    for(auto it=neighbors.positions.begin(); it!=neighbors.positions.end(); ++it)
    {
        if(particles.grid[*it]==1)//if a neighbor is diffuse
        {
            rand=(rand-orientation_a)*3-1;
            orientation_a = rand;
            int slope_a = neighbors.slopes[i];
            if(((ori+slope_a) * (orientation_a+slope_a))!=0)
            {
                orientations.push_back(orientation_a);
                //slopes.push_back(slope_a);
                possible_binding_pos.push_back(*it);

                no_dif++;
            }

        }
        if(particles.grid[*it]>1)
        {
            orientation_a = particles.grid[*it]-3;
            int slope_a = neighbors.slopes[i];
            if(((ori+slope_a) * (orientation_a+slope_a))!=0)
            {
                orientations.push_back(orientation_a);
                //slope_orientations.push_back(slope_a);
                possible_binding_pos.push_back(*it);
            }
        }

        i++;
    }
    std::cout<<"pos orientation "<<ori;
    std::cout<<"\n orientation ";
    for (auto iter = orientations.begin(); iter !=orientations.end(); ++iter)
    {
        std::cout << *iter << ","<<' ';
    }
//    std::cout<<"\n slope_orientation ";
//    for (auto iter = slope_orientations.begin(); iter !=slope_orientations.end(); ++iter)
//    {
//        std::cout << *iter << ","<<' ';
//    }
    std::cout<<"\n possible_binding_pos ";
    for (auto iter = possible_binding_pos.begin(); iter !=possible_binding_pos.end(); ++iter)
    {
        std::cout << *iter << ","<<' ';
    }
    std::cout<<"\n no of diff "<< no_dif;
    std::cout<<"\n";



    if(orientations.size()>0)
    {
        double delta = delta_H(alpha,J,no_dif,orientations);
        std::cout<<"delta "<< delta << "\n";
        rand = rand-orientation_a;

        //the binding attempt is succesful if rand < delta_E
        if(rand<delta)
        {
            std::cout<<"success \n";
            particles.grid[pos]=ori+3;
            int i=0;
            for(auto it=possible_binding_pos.begin(); it!=possible_binding_pos.end(); ++it)
            {
                particles.grid[*it]=orientations[i]+3;
                auto search_pos=std::find(particles.diffuse_pos.begin(), particles.diffuse_pos.end(), *it);
                if( *search_pos==*it)
                {
                    particles.diffuse_pos.erase(search_pos);
                }
                if(!(std::find(particles.bound_pos.begin(), particles.bound_pos.end(), *it) != particles.bound_pos.end()))
                {
                    particles.bound_pos.push_back(*it);
                }
                i++;
            }

        }
    }

}

//somehow a negative 3 appears






int main()
{
    const int MC_steps = 10; // number of Monte Carlo Steps
    int MC_counter = 0;

//constants for reaction:
    int const alpha=1;
    int const J=1;

// Input and Output arrays

    int site;
    int new_particle_site;

    particles.grid[1]=1;
    particles.grid[2]=1;
    particles.grid[3]=1;
    particles.grid[9]=3;
    particles.grid[5]=3;
    particles.grid[10]=3;

    particles.get_diffuse_pos();
    particles.get_bound_pos();


    while(MC_counter<MC_steps)
    {

        std::cout <<"\n counter"<< MC_counter <<'\n';
//step 1: Move diffusive particles

        site=diffuse(particles);
        std::cout<<"site"<<site<<"\n";





//step 2: Check and create a new particles

        new_particle_site=create_particle(particles);
        std::cout<<"new particle created"<<" "<<new_particle_site<<"\n";


//step 3: bind and unbind
//BINDING
//first take site and check if particle moved adjacent to particle and perform binding_attempt
        binding_attempt(particles,alpha,J,13);

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
