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

    double J,alpha;
    int slurm_index,exponent;
    if(argc==5)
    {
        J= atof(argv[1]);
        alpha = atof(argv[2]);
        exponent= atof(argv[3]);
        slurm_index = atof(argv[4]);
    }
    else
    {
        alpha=0.1;
        J=2.6;
        exponent=6;
        slurm_index = 21;
        std::cout << "Using default parameters." << '\n';
    }
    density = 0.1;
    const int MC_steps = 5*pow(10,6); // number of Monte Carlo Steps
//    const int MC_steps =500;
    int MC_counter = 0;
//    long double rand;
    double rand;
//    long double rand_size;
    double rand_size;
    double gamma = 1.8*pow(10,-exponent);
    double const  k_off=5*pow(10,-5);
    //double const k_on=k_off/55;
    double const k_on=2*pow(10,-(exponent+1));


    Lattice lattice;
    Particles particles(lattice);
///////////////////////////
//    std::ostringstream fn;
//    fn << slurm_index << "growth_dFrzB_bondspercluster_J_" << J << "_alpha_" << alpha << ".txt";//k_un << "_" << k << ".txt";
//    std::ofstream out;
//    out.open(fn.str());
//    std::ostringstream fn2;
//    fn2 << slurm_index << "growth_dFrzB_labels_J_" << J << "_alpha_" << alpha << ".txt";// k_un << "_" << k << ".txt";
//    std::ofstream out2;
//    out2.open(fn2.str());
    std::ostringstream fn;
    fn << slurm_index << "_growth_dFrzB_bondspercluster_J_" << J << "_alpha_" << alpha << "_exp_" << exponent << ".txt";//k_un << "_" << k << ".txt";
    std::ofstream out;
    out.open(fn.str());
    std::ostringstream fn2;
    fn2 << slurm_index << "_growth_dFrzB_labels_J_" << J << "_alpha_" << alpha << "_exp_" << exponent << ".txt";// k_un << "_" << k << ".txt";
    std::ofstream out2;
    out2.open(fn2.str());

//    out << Nx << '\t' << Ny << '\t' <<'\n';
//    out2 << Nx << '\t' << Ny << '\t' <<'\n';
///////////////////////////



    std::ofstream MyFile("counters.txt");
//    std::ofstream File_Labels("labels.txt");
//    std::ofstream File_grid("grid_J_large_half_0.txt");
////    MyFile << "Nx "  << Nx << ", Ny "<<Ny<<"\n";
//    File_grid << "Nx "  << Nx << ", Ny "<<Ny<<"\n";
    while(MC_counter<MC_steps)
    {


        for(unsigned int attempted_moves=0; attempted_moves<particles.particles.size(); attempted_moves++)
        {
//        std::cout<<attempted_moves<<"\n";
//            for(unsigned int i =0; i<particles.positions.size(); i++)
//            {
//                int p=particles.get_pos(i);
//                std::cout<<"pos \t"<<particles.positions[i]<<"\t"<<particles.particles[i]->pos<<"\t ori"<<particles.grid[p]<<"\t"<<particles.particles[i]->ori<<"\n";
//            }
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
            rand_size = rand*particles.particles.size();
            int ind=rand_size;
//            std::cout<<"ind \t"<<ind<<"\n";
            rand=rand_size-ind;
            assert(rand<1);
            if(particles.is_diffuse(particles.get_pos(ind)))
            {

//Move diffusive particles
                if(rand<0)
                {

                }
                particles.attempt_diffusion(ind, rand);
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

        }
//
        if(Nx<300)
        {

        if(rand<gamma*Nx){particles.cell_growth(rand,Nx);}}

//        particles.attempt_destruction(rand);
//
        if(particles.particles.size()<1501){particles.attempt_creation(k_on);}



        if(MC_counter%10000==0)
        {
            std::cout<<Nx<<"num of particles \t"<<particles.particles.size()<<'\n';

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
        }
        MC_counter++;
    }
//    std::cout<<"START"<<"\n";
//    for(int i =0; i<particles.positions.size(); i++)
//    {
//        int p=particles.get_pos(i);
//        std::cout<<"pos \t"<<particles.positions[i]<<"\t"<<particles.particles[i]->pos<<"\t ori"<<particles.grid[p]<<"\t"<<particles.particles[i]->ori<<"\n";
//    }
////
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
