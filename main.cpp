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







int main()
{
    const int MC_steps = pow(10,5); // number of Monte Carlo Steps
    int MC_counter = 0;
    double rand;
    double rand_size;

//constants for reaction:
    double const alpha=0.0;
    double const J=1;

    Lattice lattice;
    Particles particles(lattice);

    std::ofstream MyFile("grid_rand_J_1_alpha_0.txt");
    MyFile << "Nx "  << Nx << ", Ny "<<Ny<<"\n";

    while(MC_counter<MC_steps)
    {
        std::cout <<"counter"<< MC_counter <<'\n';
        if(MC_counter%10==0)
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
                particles.attempt_creation(pos1,rand1);
//                print_container(particles.positions);
            }



            //DESTRUCTION ATTEMPT
            else if(particles.is_diffuse(particles.get_ind(pos1)))
            {
                particles.attempt_destruction(pos1,rand1);
            }
        }

//      choose random particle
        rand = unidist(gen);
        rand_size = rand*particles.positions.size();
        int ind=rand_size;
        rand=rand_size-ind;
        assert(rand<0);
        if(particles.is_diffuse(particles.get_pos(ind)))
        {

//Move diffusive particles

            particles.attempt_diffusion(rand, ind);
//            print_container(particles.positions);


//BINDING
            particles.attempt_binding(alpha,J,ind,rand);
//            print_container(particles.positions);



        }

        else if(particles.is_bound(particles.get_pos(ind)))
        {
// UNBINDING ATTEMPT
            particles.attempt_unbinding(alpha,J,ind,rand);
//            print_container(particles.positions);

        }

        if (MyFile.is_open())
        {
            for(int count = 0; count < Nxy; count ++)
            {
                MyFile << particles.grid[count] << " " ;
            }
            MyFile << "\n";
        }
        else std::cout << "Unable to open file";
        MC_counter++;
        std::cout << particles.positions.size() << '\n';

//            if(MC_counter%100==0)
//            {
//                std::cout<<"counter "<<MC_counter<<"\n";
//                std::cout<<
//            }

    }
//    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
//    {
//        //std::cout << *iter << "\n"<<' ';
//    }
    MyFile.close();
    return 0;
};
