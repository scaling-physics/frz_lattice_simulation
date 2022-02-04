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

//# define NDEBUG //comment out to turn on assert.
# include <assert.h>	// for assert()

#include "lattice.h"
#include "Particles.h"

//random number generator

int main(int argc,char *argv[])
{
    double J,alpha,rate;
    int slurm_index;
    int FrzB_num;
    if(argc==6)
    {
        J= atof(argv[1]);
        alpha=atof(argv[2]);
        FrzB_num=atof(argv[3]);
        rate=atof(argv[4]);
        slurm_index = atof(argv[5]);
    }
    else
    {
        J= 4;
        alpha=0.5;
        FrzB_num=210;
        rate=0.03;
        slurm_index = 5;
        std::cout << "Using default parameters." << '\n';
    }
    titration_concentration_frzb= 1;
    density = 0.2;

    const long int MC_steps = 1*pow(10,6); // number of Monte Carlo Steps
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
    fn << "aaFrzB_J_" << J << "_alpha_" << alpha << "_FrzB_" <<FrzB_num<<"_off_" <<rate<<"_"<< slurm_index << ".txt";//k_un << "_" << k << ".txt";
    std::ofstream out;
    out.open(fn.str());

    std::ostringstream fn2;
    fn2 << "aaFrzB_labels_J_" << J << "_alpha_" << alpha << "_FrzB_" <<FrzB_num<<"_off_" <<rate<<"_"<< slurm_index << ".txt";// k_un << "_" << k << ".txt";
    std::ofstream out2;
    out2.open(fn2.str());

    std::ostringstream fn3;
    fn3 << "aaFrzB_flags_J_" << J << "_alpha_" << alpha << "_FrzB_" <<FrzB_num<<"_off_" <<rate<<"_"<< slurm_index << ".txt";// k_un << "_" << k << ".txt";
    std::ofstream out3;
    out3.open(fn3.str());

    out << Nx << '\t' << Ny << '\t' <<'\n';
    out2 << Nx << '\t' << Ny << '\t' <<'\n';
    out3 << Nx << '\t' << Ny << '\t' <<'\n';
///////////////////////////
//    int c=0;
//    int d=0;
//    for(unsigned int i=0; i<particles.particles.size();i++)
//    {
//        int pos = particles.get_pos(i);
//        std::vector<int> free_sites(particles.get_free_flag_sites(pos));
//
//        std::cout<<pos<<"\t"<<free_sites.size()<<"\n";
//
//
//        if(free_sites.size()==2){c++;}
//        else if(free_sites.size()==1){d++;}
//
//    }
//
//    std::vector<int> free_sites(particles.get_free_flag_sites(258));
//    std::cout<<"c= \t"<<c<<"\n";
//    std::cout<<"d= \t"<<d<<"\n";
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

////////////////// FrzB Un-/Binding Section //////////////////////////

//      choose random particle for FrzB binding/unbinding
            rand = unidist(gen);
            rand_size = rand*particles.particles.size();
            int ind=rand_size;
            rand = rand_size-ind;
            assert(rand<1);
            //BINDING OF FrzB
            particles.binding_FrzB(ind, FrzB_num, rand);

            rand = unidist(gen);
            rand_size = rand*particles.particles.size();
            ind=rand_size;
            rand = rand_size-ind;
            //UNBINDING OF FrezB
            particles.unbinding_FrzB(ind, FrzB_num,rate, rand);
//////////////////////////////////////////////////////////////////////

//////////////////// Particle Movement, Un-/Binding to Clusters etc ///

//      choose random particle for interaction and movement
            rand = unidist(gen);
            rand_size = rand*particles.particles.size();
            ind=rand_size;
            rand=rand_size-ind;
            assert(rand<1);
            if(particles.is_diffuse(particles.get_pos(ind)))
            {
//Move diffusive particles
                particles.attempt_diffusion(ind, rand);

//BINDING
                particles.attempt_binding(alpha,J,ind,rand);
            }

            else if(particles.is_bound(particles.get_pos(ind)))
            {
// UNBINDING ATTEMPT
                particles.attempt_unbinding(alpha,J,ind,rand);
            }
//////////////////////////////////////////////////////////////////////////
        }

////////////////////////// Labelling of clusters, Printout////////////////
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
            }



            //std::cout<<labels.size()<<" "<<particles.positions.size()<<'\n';
            std::cout<<"Number of clusters: " << std::ranges::max(labels) << '\n';
//            print_container(labels);
            out << '\n';
            particles.print_labels(out2,labels);
            particles.print(out2);
            particles.print_Frz(out3);
        }
//////////////////////////////////////////////////////////////////////////

        MC_counter++;
    }



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
