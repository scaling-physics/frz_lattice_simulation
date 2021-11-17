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
    std::array<short,Nxy> grid; //0 if empty, ind+1 if occupied
    std::vector<int> positions; //stores position of particles on the grid
    std::vector<int> orientations; // stores the orientations of particles in positions, 0 if diffuse, 1,2,3 if bound
    Lattice &lattice;

    //std::vector<int> diffuse_pos; //reference to position of diffuse particles
    //std::vector<int> bound_pos; //reference to position of bound particles


    Particles(Lattice &lattice):grid{0},lattice(lattice)
    {

        grid[5]=1;
        grid[9]=2;
        grid[10]=3;
        grid[13]=4;

        positions.emplace_back(5);
        positions.emplace_back(9);
        positions.emplace_back(10);
        positions.emplace_back(13);

        orientations.emplace_back(2);
        orientations.emplace_back(2);
        orientations.emplace_back(2);
        orientations.emplace_back(0);
    }


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
        int ind=grid[pos]-1;
        return ind;
    }
    inline int get_orientation(int ind)
    {
        return orientations[ind];
    }
    inline void change_orientation(int ind, int ori)
    {
        orientations[ind]=ori;
    }
    inline void change_pos(int old_pos,int new_pos, int ind)
    {
        grid[old_pos]=0;
        grid[new_pos]=ind+1;
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

//CREATION ATTEMPT
    void creation_attempt(int pos, double rand)
    {
        double const k_on=0.5;
        if(rand<=k_on)
        {
            std::cout<<"particle created"<<"\n";
            positions.emplace_back(pos);
            orientations.emplace_back(0);
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
                std::cout<<"binding "<< '\t' << ori <<'\t' <<orientations[get_ind(pos)] << '\n';
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
            orientations[ind]=0;
            for(auto it=possible_unbinding_pos.begin(); it!=possible_unbinding_pos.end(); ++it)
            {
                orientations[get_ind(*it)]=0;
                i++;
            }

        }
    }

};


#endif // PARTICLES_H
