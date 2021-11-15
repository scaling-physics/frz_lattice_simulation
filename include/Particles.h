#ifndef PARTICLES_H
#define PARTICLES_H

# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
#include <array>
#include <algorithm>
#include "lattice.h"

class Particles

{
private:

public:
    std::array<short,Nxy> grid{0}; //0 if empty, ind+1 if occupied
    std::vector<int> positions; //stores position of particles on the grid
    std::vector<int> orientations; // stores the orientations of particles in positions, 0 if diffuse, 1,2,3 if bound

    //std::vector<int> diffuse_pos; //reference to position of diffuse particles
    //std::vector<int> bound_pos; //reference to position of bound particles

    inline bool is_free(int pos)
    {
        return grid[pos]==0;
    }

    inline bool is_diffuse(int ind)
    {
        return orientations[ind]==0;
    }

    inline bool is_bound(int ind)
    {
        return orientations[ind]>0;
    }

    inline int get_ind(int pos)
    {
        return grid[pos]-1;
    }
//COUNTING THE BOUND PARTICLES OF A GIVEN POSITION
    int count_bound_neighbors(int pos, int orientation)
    {
        Neighbours n=lattice.get_neighbors(pos);
        int count_bound = 0;
        int i=0;
        int ori = orientation -2;

        for(auto it=n.positions.begin(); it!=n.positions.end(); ++it)
        {
            if(is_bound(get_ind(*it)))
            {
                int orientation_a = orientations[get_ind(*it)]-2;
                int slope_a = n.slopes[i];
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
                    count_bound++;
                }
            }
            i++;
        }
        std::cout<<count_bound<<"\n";
    }

//CREATION ATTEMPT
    void creation_attempt(int pos, double rand)
    {
        double const k_on=0.5;
        if(rand<=k_on)
        {
            std::cout<<"particle created"<<"\n";
            positions.push_back(pos);
            orientations.push_back(0);
            grid[pos]=positions.size();
        }
    }
//DESTRUCTION ATTEMPT
    void destruction_attempt(int pos,double rand)
    {
        double const  k_off=0.5;
        if(rand<k_off)
        {
            std::cout<<"particle destroyed"<<"\n";
            int ind = get_ind(pos);
            grid[pos]=0;
            positions.erase(positions.begin()+ind);
            orientations.erase(orientations.begin()+ind);
        }
        else{}
    }
//DIFFUSE PARTICLES
    int diffuse(double &rand, int ind)
    {
        // get diffuse particle
        int particle_pos = positions[ind] ;
        std::cout<<"particle"<<particle_pos<<"\n";
        //choose random direction
        Neighbours neighbors_dif=lattice.get_neighbors(particle_pos);

        double rand_size=rand*neighbors_dif.positions.size();
        int dir=rand;
        rand = rand_size-dir;
        int new_pos = neighbors_dif.positions[dir];
        std::cout<<"new_pos"<<new_pos<<"\n";
        //check if new hex is occupied
        if(is_free(new_pos))// accept move
        {
            grid[particle_pos]=0;
            grid[new_pos]=ind+1;
            positions[ind]=new_pos;
            return new_pos;

        }
        else
        {
            return particle_pos;
        }
    }

//DELTA_H
    double delta_H(double const alpha, double const J, int no_dif, std::vector<int> ori_neighbors)
    {
        std::cout<<"\n no_dif "<<no_dif<<"\n ori_size "<< ori_neighbors.size()<<"\n";
        double delta_E=(-1*J)*(ori_neighbors.size()-(no_dif*alpha));
        return delta_E<0.0f ? 1.0 : exp(-delta_E);
    }

//BINDING ATTEMPT
    void binding_attempt(int const alpha, int const J, int pos, double &rand)
    {
        double rand_size = rand*3;
        int ori = rand_size-1;
        rand=rand_size-ori;
        //std::cout<<"creation rand"<<rand<<"\n";

        std::vector<int> ori_neighbors;
        std::vector<int> slopes;
        std::vector<int> possible_binding_pos;

        int no_dif = 1;
        int no_bonds = 0;
        int i=0;

        Neighbours neighbors=lattice.get_neighbors(pos);
        int orientation_a=ori;


        for(auto it=neighbors.positions.begin(); it!=neighbors.positions.end(); ++it)
        {
            if(is_diffuse(get_ind(*it)))//if a neighbor is diffuse
            {
                rand_size=(rand)*3;
                orientation_a = rand_size-1;
                rand = rand_size-orientation_a;
                int slope_a = neighbors.slopes[i];
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
                    no_bonds++;
                    ori_neighbors.push_back(orientation_a);
                    //slopes.push_back(slope_a);
                    possible_binding_pos.push_back(*it);

                    no_dif++;
                }

            }
            if(is_bound(get_ind(*it)))
            {
                orientation_a = orientations[get_ind(*it)]-2;
                int slope_a = neighbors.slopes[i];
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
                    no_bonds++;
                    ori_neighbors.push_back(orientation_a);
                    //slope_orientations.push_back(slope_a);
                    possible_binding_pos.push_back(*it);
                }
            }

            i++;
        }
        std::cout<<"pos orientation "<<ori;
        std::cout<<"\n orientation ";
        for (auto iter = ori_neighbors.begin(); iter !=ori_neighbors.end(); ++iter)
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
            std::cout<<"\n rand "<<rand;
            //the binding attempt is succesful if rand < delta_E
            if(rand<delta)
            {
                std::cout<<"success \n";
                orientations[get_ind(pos)]=ori+2;
                int i=0;
                for(auto it=possible_binding_pos.begin(); it!=possible_binding_pos.end(); ++it)
                {
                    orientations[get_ind(*it)]=ori_neighbors[i]+2;
                    i++;
                }

            }
        }

    }

//UNBINDING ATTEMPT
    void unbinding_attempt(int const alpha, int const J, int ind,double &rand)
    {
        // get bound particle
        int particle_pos = positions[ind];
        int ori = orientations[ind];
        Neighbours neighbors=lattice.get_neighbors(particle_pos);
        std::vector<int> ori_neighbors;
        std::vector<int> possible_binding_pos;
        int i=0;
        for(auto it=neighbors.positions.begin(); it!=neighbors.positions.end(); ++it)
        {
            if(is_bound(get_ind(*it)))
            {
                int orientation_a = orientations[get_ind(*it)]-2;
                int slope_a = neighbors.slopes[i];
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
                    ori_neighbors.push_back(orientation_a);
                    //slope_orientations.push_back(slope_a);
                    possible_binding_pos.push_back(*it);
                }
            }
            i++;
        }
    }

};


#endif // PARTICLES_H
