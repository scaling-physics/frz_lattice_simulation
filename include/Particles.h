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

    inline int get_index(int pos)
    {
        return grid[pos]-1;
    }











//    void get_diffuse_pos()
//    {
//        for (int j=0; j<Nxy; j++)
//        {
//            if (grid[j]==1)
//            {
//                if(!(std::find(diffuse_pos.begin(), diffuse_pos.end(), j) != diffuse_pos.end()))
//                {
//                    diffuse_pos.push_back(j);
//                }
//            }
//        }
//        //sort(diffuse_pos.begin(),diffuse_pos.end());
//    }
//
//    void get_bound_pos()
//    {
//        for (int j=0; j<Nxy; j++)
//        {
//            if (grid[j]>1)
//            {
//                if(!(std::find(bound_pos.begin(), bound_pos.end(), j) != bound_pos.end()))
//                {
//                    bound_pos.push_back(j);
//                }
//            }
//        }
//        sort(bound_pos.begin(),bound_pos.end());
//    }

//DIFFUSE PARTICLES
    int diffuse(double &rand)
    {
//randomely choose diffuse_pos
        double rand_size = rand*positions.size();
        int ind=rand_size;
        rand = rand_size-ind;
        int particle_pos = positions[ind] ;
        std::cout<<"particle"<<particle_pos<<"\n";
//choose random direction
        Neighbours neighbors_dif=lattice.get_neighbors(particle_pos);

        rand_size=rand*neighbors_dif.positions.size();
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

    double delta_H(double const alpha, double const J, int no_dif, std::vector<int> ori_neighbors)
    {
        std::cout<<"\n no_dif "<<no_dif<<"\n ori_size "<< ori_neighbors.size()<<"\n";
        double delta_E=(-1*J)*(ori_neighbors.size()-(no_dif*alpha));
        return delta_E<0.0f ? 1.0 : exp(-delta_E);
    }


    void binding_attempt(int const alpha, int const J, int pos,double &rand)
    {
        double rand_size = rand*3;
        int ori = rand_size-1;
        rand=rand_size-ori;
        //std::cout<<"creation rand"<<rand<<"\n";

        std::vector<int> ori_neighbors;
        std::vector<int> slopes;
        std::vector<int> possible_binding_pos;

        int no_dif = 1;
        int i=0;

        Neighbours neighbors=lattice.get_neighbors(pos);
        int orientation_a=ori;


        for(auto it=neighbors.positions.begin(); it!=neighbors.positions.end(); ++it)
        {
            if(is_diffuse(get_index(*it)))//if a neighbor is diffuse
            {
                rand_size=(rand)*3;
                orientation_a = rand_size-1;
                rand = rand_size-orientation_a;
                int slope_a = neighbors.slopes[i];
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
                    ori_neighbors.push_back(orientation_a);
                    //slopes.push_back(slope_a);
                    possible_binding_pos.push_back(*it);

                    no_dif++;
                }

            }
            if(is_bound(get_index(*it)))
            {
                orientation_a = orientations[get_index(*it)]-2;
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
                orientations[get_index(pos)]=ori+2;
                int i=0;
                for(auto it=possible_binding_pos.begin(); it!=possible_binding_pos.end(); ++it)
                {
                    orientations[get_index(*it)]=ori_neighbors[i]+2;
                    i++;
                }

            }
        }

    }

    void unbinding_attempt(int const alpha, int const J, int pos,double &rand)
    {

    }

};


#endif // PARTICLES_H
