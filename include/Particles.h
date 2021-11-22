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



void print_container(const std::vector<int>& c)
{
    for (auto &i : c)
    {
        std::cout << i << " ";
    }
    std::cout << '\n';
}

struct Diffuse_Neighbors
{
    std::vector<int> diffuse_neighbors;
    std::vector<int> orientations;
    std::vector<int> slopes;
};

struct Bound_Neighbors
{
    std::vector<int> bound_neighbors;
    std::vector<int> orientations;
    std::vector<int> slopes;
};

struct Interactions
{
    std::vector<int> possible_interaction_pos;
    std::vector<int> orientations;
    int num_diffuse;
    int num_bonds;
};



class Particles

{
private:

public:
    std::array<short,Nxy> grid; //0 if empty, 1 if diffuse, {2,3,4} if bound as ori
    std::vector<int> positions; //stores position of particles on the grid
    Lattice &lattice;


    //std::vector<int> diffuse_pos; //reference to position of diffuse particles
    //std::vector<int> bound_pos; //reference to position of bound particles


    Particles(Lattice &lattice):grid{0},lattice(lattice)
    {
        for(int i=0; i<Nxy; i++)
        {
            if(i%3==0)
            {
                grid[i]=1;
                positions.emplace_back(i);
            }
        }
//        grid[5]=3;
//        grid[9]=3;
//        grid[10]=3;
//        grid[13]=1;
//
//        positions.emplace_back(5);
//        positions.emplace_back(9);
//        positions.emplace_back(10);
//        positions.emplace_back(13);
    }

    inline int get_pos(int ind)
    {
        return positions[ind];
    }
    inline bool is_free(int pos)
    {
        return grid[pos]==0;
    }
    inline bool is_diffuse(int pos)
    {
        return grid[pos]==1;
    }


    inline bool is_bound(int pos)
    {
        return grid[pos]>1;
    }

    inline int get_ind(int pos)
    {
//    NOT WORKING PROPERLY
        return pos;

    }
    inline int get_orientation(int pos)
    {
        int ori=grid[pos];
        return ori==1 ? ori : ori-3;
    }
    inline void set_orientation(int pos, int ori)
    {
        grid[pos]=ori;
    }
    inline void set_pos(int old_pos,int new_pos, int ind)
    {
        grid[old_pos]=0;
        grid[new_pos]=1;
        positions[ind]=new_pos;
    }
    inline bool is_interaction_allowed(int ori_1, int ori_2, int slope)
    {
        return ((ori_1 + slope)*(ori_2+slope))!=0;
    }

//CREATION ATTEMPT
    void creation_attempt(int pos, double rand)
    {
        double const k_on=0.5;
        if(rand<=k_on)
        {
            std::cout<<"particle created"<<"\n";
            positions.emplace_back(pos);
            grid[pos]=1;
        }
    }
//DESTRUCTION ATTEMPT
    void destruction_attempt(int pos,double rand)
    {
        double const  k_off=0.5;
        if(rand<k_off)
        {
            std::cout<<"particle destroyed"<<"\n";
            grid[pos]=0;
            positions.erase(std::find(begin(positions),end(positions),pos));
        }
    }
//DIFFUSE PARTICLES
    void diffuse(double &rand, int ind)
    {
        // get diffuse particle
        int particle_pos = positions[ind] ;

        //choose random direction
        Neighbours neighbors_dif=lattice.get_neighbors(particle_pos);

        double rand_size=rand*neighbors_dif.positions.size();
        int dir=rand_size;
        rand = rand_size-dir;
        int new_pos = neighbors_dif.positions[dir];

        //check if new hex is occupied
        if(is_free(new_pos))// accept move
        {
//            std::cout<<"particle"<<particle_pos<<"\n";
//            std::cout<<"new_pos"<<new_pos<<"\n";
            set_pos(particle_pos,new_pos, ind);
        }
    }
//COUNTING THE BOUND PARTICLES OF A GIVEN POSITION
    int count_bound_neighbors(int pos, int ori)
    {
        Neighbours n=lattice.get_neighbors(pos);
        int count_bound = 0;

        for(unsigned int i=0; i<n.positions.size(); i++)
        {
            if(!is_free(n.positions[i]))
            {

                int ind=get_ind(n.positions[i]);
                if(is_bound(ind))
                {
                    int orientation_a = get_orientation(ind)-2;
                    int slope_a = n.slopes[i];
                    if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                    {
                        count_bound++;
                    }
                }

            }
            i++;

        }
//        std::cout<<count_bound<<"\n";
        return count_bound;
    }
//GET DIFFUSE NEIGHBORS, RANDOMELY ASSIGN A ORIENTATION & GET THEIR SLOPES
    Diffuse_Neighbors get_diffuse_neighbors(int ind, double &rand)
    {
        Diffuse_Neighbors d_n;
        int pos = get_pos(ind);
        Neighbours n=lattice.get_neighbors(pos);
        std::vector<int> diffuse_neighbors;
        double rand_size;
        int ori;
        for(unsigned int i=0; i<n.positions.size(); i++)
        {
            if(is_diffuse(n.positions[i]))
            {
                d_n.diffuse_neighbors.emplace_back(n.positions[i]);
                d_n.slopes.emplace_back(n.slopes[i]);
                rand_size = rand*3;
                ori = rand_size;
                rand=rand_size-ori;
                d_n.orientations.emplace_back(ori-1);
            }
        }
        return d_n;
    }
//GET BOUND NEIGHBORS, THEIR ORIENTATIONS AND THEIR SLOPES
    Bound_Neighbors get_bound_neighbors(int ind)
    {
        Bound_Neighbors b_n;
        int pos = get_pos(ind);
        Neighbours n=lattice.get_neighbors(pos);
        for(unsigned int i=0; i<n.positions.size(); i++)
        {
            if(is_bound(n.positions[i]))
            {
                b_n.bound_neighbors.emplace_back(n.positions[i]);
                b_n.orientations.emplace_back(get_orientation(n.positions[i]));
                b_n.slopes.emplace_back(n.slopes[i]);
            }
        }
        return b_n;
    }

//DELTA_H
    double delta_H(double const alpha, double const J, int num_dif, int num_bonds, double a)
    {
        //std::cout<<"num_dif "<<num_dif<<"\n num_bonds "<< num_bonds<<"\n";
        double delta_E=(a*J)*(num_bonds-(num_dif*alpha));
        return delta_E<0.0f ? 1.0 : exp(-delta_E);
    }
//BINDING ATTEMPT
    void binding_attempt(int const alpha, int const J, int ind, double &rand)
    {
        int pos = get_pos(ind);
        Interactions interactions;
        Diffuse_Neighbors d_n=get_diffuse_neighbors(ind,rand);
        Bound_Neighbors b_n=get_bound_neighbors(ind);

        double rand_size = rand*3;
        int rand_int = rand_size;
        int ori = rand_int-1;
        rand=rand_size-rand_int;


        interactions.num_bonds=0;
        interactions.num_diffuse=1;
        Neighbours neighbors=lattice.get_neighbors(pos);
        if(d_n.diffuse_neighbors.size()>0)
        {
            for(unsigned int i=0; i<d_n.diffuse_neighbors.size(); i++)
            {
                if(is_interaction_allowed(ori,d_n.orientations[i],d_n.slopes[i]))
                {
                    interactions.possible_interaction_pos.emplace_back(d_n.diffuse_neighbors[i]);
                    interactions.orientations.emplace_back(d_n.orientations[i]);
                    interactions.num_bonds++;
                    interactions.num_diffuse++;
                    interactions.num_bonds=interactions.num_bonds+count_bound_neighbors(d_n.diffuse_neighbors[i],d_n.orientations[i]);

                }
            }
        }
        if(b_n.bound_neighbors.size()>0)
        {
            for(unsigned int i=0; i<b_n.bound_neighbors.size(); i++)
            {
                if(is_interaction_allowed(ori,b_n.orientations[i],b_n.slopes[i]))
                {
                    interactions.orientations.emplace_back(b_n.orientations[i]);
                    interactions.possible_interaction_pos.emplace_back(b_n.bound_neighbors[i]);
                    interactions.num_bonds++;
                }
            }
        }
//            else if(is_diffuse(ind))//if a neighbor is diffuse
//            {
//                rand_size=(rand)*3;
//                rand_int=rand_size;
//                orientation_a = rand_int-1;
//                rand = rand_size-rand_int;
//                int slope_a = neighbors.slopes[i];
////                std::cout<<"neighbor_dif "<<*it<<"\n";
////                std::cout<<"ori "<<orientation_a<<"\n";
//                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
//                {
//                    num_bonds++;
//                    ori_neighbors.emplace_back(orientation_a);
//                    possible_binding_pos.emplace_back(*it);
//
//                    num_dif++;
//                    //check if diffuse neighbor has bound neighbors
//                    int num_neigh_bound = count_bound_neighbors(*it,orientation_a);
//                    num_bonds=num_bonds+num_neigh_bound;
//                }
//
//            }
//            else if(is_bound(ind))
//            {
//                orientation_a = get_orientation(ind)-2;
//                int slope_a = neighbors.slopes[i];
////                std::cout<<"neighbor_bound "<<*it<<"\n";
////                std::cout<<"ori "<<orientation_a<<"\n";
//                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
//                {
//                    num_bonds++;
//                    ori_neighbors.emplace_back(orientation_a);
//                    possible_binding_pos.emplace_back(*it);
//                }
//            }
//
//            i++;
//        }
        if(interactions.num_bonds>0)
        {
            double delta = delta_H(alpha,J,interactions.num_diffuse,interactions.num_bonds, -1);
            //std::cout<<"delta "<< delta << "\n";
            //std::cout<<"\n rand "<<rand;
            //the binding attempt is succesful if rand < delta_E
            if(rand<delta)
            {
                std::cout<<"binding"<<"\n";
                set_orientation(pos,ori+3);
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],interactions.orientations[i]+3);
                }

            }
        }

    }

//UNBINDING ATTEMPT
    void unbinding_attempt(int const alpha, int const J, int ind,double &rand)
    {
        // get bound particle
        int particle_pos = get_pos(ind);
        int ori = get_orientation(particle_pos);
        Neighbours neighbors=lattice.get_neighbors(particle_pos);

        Interactions interactions;
        interactions.num_bonds=0;
        interactions.num_diffuse=1;
        Bound_Neighbors b_n=get_bound_neighbors(ind);


        if(b_n.bound_neighbors.size()>0)
        {
            for (unsigned int i=0; i<b_n.bound_neighbors.size(); i++)
            {
                if(is_interaction_allowed(ori,b_n.orientations[i],b_n.slopes[i]))
                {
                    interactions.orientations.emplace_back(b_n.orientations[i]);
                    interactions.num_bonds++;

                    int bound_neigh_of_neigh = count_bound_neighbors(b_n.bound_neighbors[i],b_n.orientations[i]);
                    if(bound_neigh_of_neigh==1)
                    {
                        interactions.possible_interaction_pos.emplace_back(b_n.bound_neighbors[i]);
                        interactions.num_diffuse++;
                    }

                }
            }
        }

        double delta_E =delta_H(alpha,J,interactions.num_diffuse,interactions.num_bonds, 1);

        if (rand<delta_E)
        {
            std::cout<<"unbinding"<<"\n";
            set_orientation(particle_pos,1);
            if(interactions.possible_interaction_pos.size()>0)
            {
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],1);
                }
            }

        }
    }

};

#endif // PARTICLES_H
