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
#include "Particles.h"
#include "definitions.h"
//using namespace std;

//random number generator

unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);






int main()
{
    const int MC_steps = 2; // number of Monte Carlo Steps
    int MC_counter = 0;
    double rand;

//constants for reaction:
    double const alpha=0.5;
    double const J=1.0;

// Input and Output arrays

    int site;

    Lattice lattice;
    Particles particles(lattice);

    while(MC_counter<MC_steps)
    {
        std::cout <<"counter"<< MC_counter <<'\n';
        if(MC_steps%1==0)
        {
            // select a random hex on grid:
            double rand_size1 = unidist(gen) * Nxy;
            int pos1 = rand_size1;
            double rand1 = rand_size1-pos1;
//            std::cout<<"pos "<<pos<<"\n";
//            std::cout<<"rand "<<rand<<"\n";




            //CREATION ATTEMPT
            if(particles.is_free(pos1))
            {
                //std::cout<<"free"<<"\n";
                particles.creation_attempt(pos1,rand1);
                print_container(particles.positions);
            }



            //DESTRUCTION ATTEMPT
            else if(particles.is_diffuse(particles.get_ind(pos1)))
            {
                std::cout<<"diffuse "<<"\n";
                particles.destruction_attempt(pos1,rand1);

            }
        }

//      choose random particle
        rand = unidist(gen)*particles.positions.size();
//        std::cout<<"rand"<<rand<<"\n";
        int ind=rand;
        std::cout<<ind<<"\n";
        rand=rand-ind;

        if(particles.is_diffuse(ind))
        {

//Move diffusive particles

            site=particles.diffuse(rand, ind);
            print_container(particles.positions);


//BINDING
            particles.binding_attempt(alpha,J,site,rand);
            print_container(particles.positions);



        }

        else if(particles.is_bound(ind))
        {
// UNBINDING ATTEMPT
            particles.unbinding_attempt(alpha,J,ind,rand);
            print_container(particles.positions);

        }


        MC_counter++;
        std::cout << particles.positions.size() << '\n';

    }
    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
    {
        //std::cout << *iter << "\n"<<' ';
    }

    return 0;
};
