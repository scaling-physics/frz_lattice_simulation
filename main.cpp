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
    const int MC_steps = pow(10,8); // number of Monte Carlo Steps
//    const int MC_steps =500;
    int MC_counter = 0;
    double rand;
    double rand_size;

//constants for reaction:
    double const alpha=0.05;
    double const J=1;

    Lattice lattice;
    Particles particles(lattice);

//    std::ofstream MyFile("grid_ori_pos.txt");
//    std::ofstream File_Labels("labels.txt");
    std::ofstream File_grid("grid_J_large_half.txt");
//    MyFile << "Nx "  << Nx << ", Ny "<<Ny<<"\n";
    File_grid << "Nx "  << Nx << ", Ny "<<Ny<<"\n";
    while(MC_counter<MC_steps)
    {
//        std::cout <<"counter"<< MC_counter <<'\n';
//        if(MC_counter%10==0)
//        {
//            // select a random hex on grid:
//            double rand_size1 = unidist(gen) * Nxy;
//            int pos1 = rand_size1;
//            double rand1 = rand_size1-pos1;
////            std::cout<<"pos "<<pos<<"\n";
////            std::cout<<"rand "<<rand<<"\n";
//
//
//            //CREATION ATTEMPT
//            if(particles.is_free(pos1))
//            {
//                //std::cout<<"free"<<"\n";
//                particles.attempt_creation(pos1,rand1);
////                print_container(particles.positions);
//            }
//
//            //DESTRUCTION ATTEMPT
//            else if(particles.is_diffuse(particles.get_ind(pos1)))
//            {
//                particles.attempt_destruction(pos1,rand1);
//            }
//        }

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

//        if (MC_counter>1 && MC_counter%1000==0)
//        {
//            std::vector<int> labels(particles.positions.size(),0);
//            int i =1;
//            for (unsigned int index=0; index < particles.positions.size(); index++)
//            {
//                if(particles.is_bound(particles.get_pos(index)) && labels[index]==0)
//                {
//                    particles.label(index,i,labels);
//                    i++;
//                }
//            }
//
//            std::cout<<labels.size()<<" "<<particles.positions.size()<<'\n';
//            print_container(labels);
//            particles.print_labels(File_Labels,labels);
//
//        }


//        particles.print(MyFile);

//        if (MyFile.is_open())
//        {
//            for(int count = 0; count < Nxy; count ++)
//            {
//                MyFile << particles.grid[count] << " " ;
//            }
//            MyFile << "\n";
//        }
//        else std::cout << "Unable to open file";
        MC_counter++;
//        std::cout << particles.positions.size() << '\n';

//            if(MC_counter%100==0)
//            {
//                std::cout<<"counter "<<MC_counter<<"\n";
//                std::cout<<
//            }

    }
    particles.print_grid(File_grid);
//    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
//    {
//        //std::cout << *iter << "\n"<<' ';
//    }
    File_grid.close();
    return 0;
};
