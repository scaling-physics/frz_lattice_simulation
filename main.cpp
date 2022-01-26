# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <array>
#include <vector>
# include <chrono>
# include <random>
#include <algorithm>
#include <ranges>

# define NDEBUG //comment out to turn on assert.
# include <assert.h>	// for assert()

#include "lattice.h"
#include "Particles.h"

//random number generator

int main(int argc,char *argv[])
{

    double J,alpha,rate;
    int slurm_index;
    if(argc==3)
    {
        titration_concentration_frzb= atof(argv[1]);
        slurm_index = atof(argv[2]);
    }
    else
    {
        titration_concentration_frzb=20;
        slurm_index = 5;
        std::cout << "Using default parameters." << '\n';
    }

    alpha=0.0;
    J=2.6;
    density = 0.2;
    rate=0.5;

    const long int MC_steps = pow(10,6); // number of Monte Carlo Steps
//    const int MC_steps =500;
    long int MC_counter = 0;
//    long double rand;
    double rand;
//    long double rand_size;
    double rand_size;

    Lattice lattice;
    Particles particles(lattice);
///////////////////////////
    std::ostringstream fn;
    fn << "weird_FrzB_" << J << "_" << alpha << "_" <<titration_concentration_frzb<<"_"<< slurm_index << ".txt";//k_un << "_" << k << ".txt";
    std::ofstream out;
    out.open(fn.str());

    std::ostringstream fn2;
    fn2 << "weird_FrzB_labels_" << J << "_" << alpha<<"_"  <<titration_concentration_frzb<<"_" << slurm_index << ".txt";// k_un << "_" << k << ".txt";
    std::ofstream out2;
    out2.open(fn2.str());

    std::ostringstream fn3;
    fn3 << "weird_FrzB_long_flags_" << J << "_" << alpha<<"_"  <<titration_concentration_frzb<<"_" << slurm_index << ".txt";// k_un << "_" << k << ".txt";
    std::ofstream out3;
    out3.open(fn3.str());

    out << Nx << '\t' << Ny << '\t' <<'\n';
    out2 << Nx << '\t' << Ny << '\t' <<'\n';
    out3 << Nx << '\t' << Ny << '\t' <<'\n';
///////////////////////////



    std::ofstream MyFile("counters.txt");
//    std::ofstream File_Labels("labels.txt");
//    std::ofstream File_grid("grid_J_large_half_0.txt");
////    MyFile << "Nx "  << Nx << ", Ny "<<Ny<<"\n";
//    File_grid << "Nx "  << Nx << ", Ny "<<Ny<<"\n";
    while(MC_counter<MC_steps)
    {
//        std::cout<<"MC_counter "<< MC_counter<<'\n';

        for(unsigned int attempted_moves=0; attempted_moves<particles.particles.size(); attempted_moves++)
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

//      choose random particle for FrzB binding/unbinding
            rand = unidist(gen);
            rand_size = rand*particles.particles.size();
            int ind=rand_size;
            assert(rand<1);
           //BINDING OF FrzB
            particles.binding_FrzB(ind, FrzB_num, rand);

            rand = unidist(gen);
            rand_size = rand*particles.particles.size();
            ind=rand_size;
            //UNBINDING OF FrezB
            particles.unbinding_FrzB(ind, FrzB_num,rate, rand);



//      choose random particle for interaction and movement
            rand = unidist(gen);
            rand_size = rand*particles.particles.size();
            ind=rand_size;
            rand=rand_size-ind;
            assert(rand<1);
            if(particles.is_diffuse(particles.get_pos(ind)))
            {
//                std::cout<<" diffuse pos "<<particles.get_pos(ind)<<"\n";


//Move diffusive particles
                particles.attempt_diffusion(ind, rand);
//            print_container(particles.positions);


//BINDING
                particles.attempt_binding(alpha,J,ind,rand);
//            print_container(particles.positions);



            }

            else if(particles.is_bound(particles.get_pos(ind)))
            {
//            std::cout<<" bound pos "<<particles.get_pos(ind)<<"\n";
// UNBINDING ATTEMPT
                particles.attempt_unbinding(alpha,J,ind,rand);
//            print_container(particles.positions);

            }





        }
        if(MC_counter%20000==0)
        {

//                std::cout<<"Number of bonds: ";
            std::vector<int> labels(particles.particles.size(),0);
            int label_i = 1;
            for (unsigned int label_index=0; label_index < particles.particles.size(); label_index++)
            {
                if(particles.is_bound(particles.get_pos(label_index)) && labels[label_index]==0)
                {
                    int num_bonds=0;
                    particles.label(label_index,label_i,labels,num_bonds);
                    label_i++;
//                        std::cout << num_bonds << '\t';
                    out << num_bonds<< '\t';
                }
//                int check_num=0;
//                for(unsigned int check=0; check<labels.size(); check++)
//                {
//                    if(labels[check]==label_i-1)
//                    {
//                        check_num++;
//                    }
//                    //assert(check_num>1);
//                    if(check_num==1){
//                    std::cout<<label_i-1;
//
//                    }
//                }
            }



            //std::cout<<labels.size()<<" "<<particles.positions.size()<<'\n';
            std::cout<<"Number of clusters: " << std::ranges::max(labels) << '\n';
//            print_container(labels);
            out << '\n';
            particles.print_labels(out2,labels);
            particles.print(out2);
            particles.print_Frz(out3);
        }
        MC_counter++;
    }
//
//    std::vector<int> labels(particles.positions.size(),0);
//    int label_i = 1;
//    for (unsigned int label_index=0; label_index < particles.positions.size(); label_index++)
//    {
//        if(particles.is_bound(particles.get_pos(label_index)) && labels[label_index]==0)
//        {
//            particles.label(label_index,label_i,labels);
//            label_i++;
//        }
//    }
//
//    std::cout<<labels.size()<<" "<<particles.positions.size()<<'\n';
//    print_container(labels);
//    particles.print_labels(out2,labels);
//    particles.print(out2);
//
//    particles.print_grid(out);
//    for (auto iter = particles.grid.begin(); iter !=particles.grid.end(); ++iter)
//    {
//        //std::cout << *iter << "\n"<<' ';
//    }
    std::cout<<"diffuse attempt "<<diffuse_attempt<<"\n";
    std::cout<<"diffuse succ "<<diffuse_succ<<"\n";
    double dif_rat = static_cast<double>(diffuse_succ)/diffuse_attempt;
    std::cout<<"diffuse ratio "<<dif_rat<<"\n";
    std::cout<<"binding attempt "<<binding_attempt<<"\n";
    std::cout<<"binding succ "<<binding_succ<<"\n";
    double bind_rat = static_cast<double>(binding_succ)/binding_attempt;
    std::cout<<"binding ratio "<<bind_rat<<"\n";
    std::cout<<"unbinding attempt "<<unbinding_attempt<<"\n";
    std::cout<<"unbindind succ "<<unbinding_succ<<"\n";
    double unbind_rat = static_cast<double>(unbinding_succ)/unbinding_attempt;
    std::cout<<"unbindind ratio "<< unbind_rat<<"\n";
    double ratio_dif_steps = static_cast<double>(diffuse_attempt)/MC_steps;
    double ratio_bind_steps = static_cast<double>(binding_attempt)/MC_steps;
    double ratio_unbind_steps = static_cast<double>(unbinding_attempt)/MC_steps;
    std::cout<<"ratios\n "<<"dif/MC_steps "<< ratio_dif_steps<<"\n";
    std::cout<<"bind/MC_steps "<< ratio_bind_steps<<"\n";
    std::cout<<"unbind/MC_steps "<< ratio_unbind_steps<<"\n";

    MyFile<<"diffuse attempt "<<diffuse_attempt<<"\n";
    MyFile<<"diffuse succ "<<diffuse_succ<<"\n";
    MyFile<<"diffuse ratio "<<dif_rat<<"\n";
    MyFile<<"binding attempt "<<binding_attempt<<"\n";
    MyFile<<"binding succ "<<binding_succ<<"\n";
    MyFile<<"binding ratio "<<bind_rat<<"\n";
    MyFile<<"unbinding attempt "<<unbinding_attempt<<"\n";
    MyFile<<"unbindind succ "<<unbinding_succ<<"\n";
    MyFile<<"unbindind ratio "<< unbind_rat<<"\n";
    MyFile<<"ratios\n "<<"dif/MC_steps "<< ratio_dif_steps<<"\n";
    MyFile<<"bind/MC_steps "<< ratio_bind_steps<<"\n";
    MyFile<<"unbind/MC_steps "<< ratio_unbind_steps<<"\n";

    return 0;
};
