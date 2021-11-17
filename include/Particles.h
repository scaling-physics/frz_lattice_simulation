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

        grid[5]=3;
        grid[9]=3;
        grid[10]=3;
        grid[13]=1;

        positions.emplace_back(5);
        positions.emplace_back(9);
        positions.emplace_back(10);
        positions.emplace_back(13);
    }

    inline int get_pos(int ind)
    {
        return positions[ind];
    }
    inline bool is_free(int ind)
    {
        return grid[get_pos(ind)]==0;
    }
    inline bool is_diffuse(int ind)
    {
        return grid[get_pos(ind)]==1;
    }


    inline bool is_bound(int ind)
    {
        return grid[get_pos(ind)]>1;
    }

    inline int get_ind(int pos)
    {
//    NOT WORKING PROPERLY
        return pos;

    }
    inline int get_orientation(int ind)
    {
        int ori=grid[get_pos(ind)];
        return ori==1 ? ori : ori-3;
    }
    inline void set_orientation(int ind, int ori)
    {
        grid[get_pos(ind)]=ori;
    }
    inline void set_pos(int old_pos,int new_pos, int ind)
    {
        grid[old_pos]=0;
        grid[new_pos]=1;
        positions[ind]=new_pos;
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
    int diffuse(double &rand, int ind)
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
            std::cout<<"particle"<<particle_pos<<"\n";
            std::cout<<"new_pos"<<new_pos<<"\n";
            set_pos(particle_pos,new_pos, ind);
            return new_pos;
        }
        else
        {
            return particle_pos;
        }
    }
//COUNTING THE BOUND PARTICLES OF A GIVEN POSITION
    int count_bound_neighbors(int pos, int ori)
    {
        Neighbours n=lattice.get_neighbors(pos);
        int count_bound = 0;

        for(int i=0; i<n.positions.size(); i++)
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

//DELTA_H
    double delta_H(double const alpha, double const J, int no_dif, int no_bonds, double a)
    {
        //std::cout<<"no_dif "<<no_dif<<"\n no_bonds "<< no_bonds<<"\n";
        double delta_E=(a*J)*(no_bonds-(no_dif*alpha));
        return delta_E<0.0f ? 1.0 : exp(-delta_E);
    }
//BINDING ATTEMPT
    void binding_attempt(int const alpha, int const J, int pos, double &rand)
    {
        double rand_size = rand*3;
        int rand_int = rand_size;
        int ori = rand_int-1;
        rand=rand_size-rand_int;
        std::cout<<"particle"<<pos<<"\n";
        std::cout<<"ori "<<ori<<"\n";

        std::vector<int> ori_neighbors;
        std::vector<int> slopes;
        std::vector<int> possible_binding_pos;

        int no_dif = 1;
        int no_bonds = 0;
        int i=0;

        Neighbours neighbors=lattice.get_neighbors(pos);
        int orientation_a;


        for(auto it=neighbors.positions.begin(); it!=neighbors.positions.end(); ++it)
        {
            int ind = get_ind(*it);
            if(ind<0) {}
            else if(is_diffuse(ind))//if a neighbor is diffuse
            {
                rand_size=(rand)*3;
                rand_int=rand_size;
                orientation_a = rand_int-1;
                rand = rand_size-rand_int;
                int slope_a = neighbors.slopes[i];
                std::cout<<"neighbor_dif "<<*it<<"\n";
                std::cout<<"ori "<<orientation_a<<"\n";
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
                    no_bonds++;
                    ori_neighbors.emplace_back(orientation_a);
                    possible_binding_pos.emplace_back(*it);

                    no_dif++;
                    //check if diffuse neighbor has bound neighbors
                    int no_neigh_bound = count_bound_neighbors(*it,orientation_a);
                    no_bonds=no_bonds+no_neigh_bound;
                }

            }
            else if(is_bound(ind))
            {
                orientation_a = get_orientation(ind)-2;
                int slope_a = neighbors.slopes[i];
                std::cout<<"neighbor_bound "<<*it<<"\n";
                std::cout<<"ori "<<orientation_a<<"\n";
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
                    no_bonds++;
                    ori_neighbors.emplace_back(orientation_a);
                    possible_binding_pos.emplace_back(*it);
                }
            }

            i++;
        }
        if(no_bonds>0)
        {
            double delta = delta_H(alpha,J,no_dif,no_bonds, -1);
            //std::cout<<"delta "<< delta << "\n";
            //std::cout<<"\n rand "<<rand;
            //the binding attempt is succesful if rand < delta_E
            if(rand<delta)
            {
                set_orientation(get_ind(pos),ori+2);
                for(int i=0;i<possible_binding_pos.size(); i++)
                {
                    set_orientation(get_ind(possible_binding_pos[i]),ori_neighbors[i]+2);
                }

            }
        }

    }

//UNBINDING ATTEMPT
    void unbinding_attempt(int const alpha, int const J, int ind,double &rand)
    {
        // get bound particle
        int particle_pos = positions[ind];
        int ori = get_orientation(ind)-2;
        Neighbours neighbors=lattice.get_neighbors(particle_pos);
        std::vector<int> ori_neighbors;
        std::vector<int> possible_unbinding_pos;
        int i=0;
        int no_bonds=0;
        int no_dif=1;
        //print_container(neighbors.positions);
        for(auto it=neighbors.positions.begin(); it!=neighbors.positions.end(); ++it)
        {
            int ind1 = get_ind(*it);
            if(ind1>=0 && is_bound(ind1))
            {
                int orientation_a = get_orientation(ind1)-2;
                int slope_a = neighbors.slopes[i];
                if(((ori+slope_a) * (orientation_a+slope_a))!=0)
                {
//                    ori_neighbors.emplace_back(orientation_a);
                    no_bonds++;
                    int no_neigh_bound = count_bound_neighbors(*it,orientation_a);

                    if (no_neigh_bound==1)
                    {
                        no_dif++;
                        possible_unbinding_pos.emplace_back(*it);
                    }
                }
            }
            i++;
        }

        double delta_E =delta_H(alpha,J,no_dif,no_bonds, 1);

        if (rand<delta_E)
        {
            //std::cout<<"unbinding"<<"\n";
            set_orientation(ind,0);
            for(int i=0;i<possible_unbinding_pos.size();i++)
            {
                set_orientation(get_ind(possible_unbinding_pos[i]),0);
            }

        }
    }

};


#endif // PARTICLES_H
