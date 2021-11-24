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

//void label(int ind, int i,std::vector<int> &labels)
//{
//    labels[ind]=i;
//    std::vector<int> neigh=get_bonded_neighbours(ind, get_ori(ind))
//                           for(inds : neigh )
//    {
//        if(labels[inds]==0)
//        {
//            label(inds,i,labels);
//        }
//    }
//}



    int main()
    {
        const int MC_steps = 20000; // number of Monte Carlo Steps
        int MC_counter = 0;
        double rand;

//constants for reaction:
        double const alpha=0.0;
        double const J=-1000;

        Lattice lattice;
        Particles particles(lattice);

        std::ofstream MyFile("grid_J_-1000_alpha_0.txt");
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
                    particles.creation_attempt(pos1,rand1);
//                print_container(particles.positions);
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
            rand=rand-ind;

            if(particles.is_diffuse(particles.get_pos(ind)))
            {

//Move diffusive particles

                particles.diffuse(rand, ind);
//            print_container(particles.positions);


//BINDING
                particles.binding_attempt(alpha,J,ind,rand);
//            print_container(particles.positions);



            }

            else if(particles.is_bound(particles.get_pos(ind)))
            {
// UNBINDING ATTEMPT
                particles.unbinding_attempt(alpha,J,ind,rand);
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

            if(MC_counter%100==0)
            {
                std::cout<<"counter "<<MC_counter<<"\n";
                std::cout<<
            }

        }
//    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
//    {
//        //std::cout << *iter << "\n"<<' ';
//    }
        MyFile.close();
        return 0;
    };
