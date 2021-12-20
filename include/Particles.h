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

int unbinding_attempt=0;
int unbinding_succ=0;
int binding_attempt=0;
int binding_succ=0;
int diffuse_attempt=0;
int diffuse_succ=0;

double const k_on=0.05;
double const  k_off=0.5;

double density;
unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);

void print_container(const std::vector<int>& c)
{
    for (auto &i : c)
    {
        std::cout << i << " ";
    }
    std::cout << '\n';
}


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
    std::array<short,Nxy> grid; //0 if empty, 1 if diffuse, {2,3,4} if bound as ori

public:

    std::vector<int> positions; //stores position of particles on the grid
    Lattice &lattice;



    Particles(Lattice &lattice):grid{0},lattice(lattice)
    {
        int initial_num = Nxy*density;
        positions.resize(initial_num);

        for (int i=0; i<initial_num; i++)
        {
//            long double random=unidist(gen);
            double random=unidist(gen);
            int pos = random*Nxy;
            if(grid[pos]==0)
            {
                grid[pos]=1;
                positions[i]=pos;
            }
            else
            {
                i--;
            }
        }
    }

    inline int get_pos(const int ind) const
    {
        return positions[ind];
    }
    inline bool is_free(const int pos) const
    {
        return grid[pos]==0;
    }
    inline bool is_diffuse(const int pos) const
    {
        return grid[pos]==1;
    }
    inline bool is_bound(const int pos) const
    {
        return grid[pos]>1;
    }
    inline int get_orientation(const int pos) const
    {
        int ori=grid[pos];
        assert(is_bound(pos));
        return ori-3;
    }
    inline void set_orientation(const int pos, const int ori)
    {
        grid[pos]=ori;
    }
    inline void set_pos(const int old_pos,const int new_pos, const int ind)
    {
        grid[old_pos]=0;
        grid[new_pos]=1;
        positions[ind]=new_pos;
    }
    inline bool is_interaction_allowed(const int ori_1, const int ori_2, const int slope) const
    {
        return (ori_1 != ori_2 && ((ori_1 + slope)*(ori_2+slope))!=0);
    }


//CREATION ATTEMPT
    void attempt_creation(const int pos, const double &rand)
    {
        if(rand<=k_on)
        {
            std::cout<<"particle created"<<"\n";
            positions.emplace_back(pos);
            grid[pos]=1;
        }
    }


//DESTRUCTION ATTEMPT
    void attempt_destruction(const int pos,const double &rand)
    {
        if(rand<k_off)
        {
//            std::cout<<"particle destroyed"<<"\n";
            grid[pos]=0;
            positions.erase(std::find(begin(positions),end(positions),pos));
        }
    }


//DIFFUSE PARTICLES
    void attempt_diffusion(const int ind, double &rand)
    {
        // get diffuse particle
        diffuse_attempt++;
        int particle_pos = get_pos(ind);

        //choose random direction
        std::vector<Neighbour> n(lattice.get_neighbors2(particle_pos));

//        long double rand_size=rand*n.size();
        double rand_size=rand*n.size();
        int dir=rand_size;
        rand = rand_size-dir;
        int new_pos = n[dir].position;

        //check if new hex is occupied
        if(is_free(new_pos))// accept move
        {
            set_pos(particle_pos,new_pos, ind);
            diffuse_succ++;
        }
    }


//COUNTING THE BOUND PARTICLES OF A GIVEN POSITION
    int count_interacting_neighbors(const int pos, const int ori) const
    {
        std::vector<Neighbour> n(lattice.get_neighbors2(pos));
        int count_bound = 0;

        for(unsigned int i=0; i<n.size(); i++)
        {
            if(is_bound(n[i].position) && is_interaction_allowed(ori,get_orientation(n[i].position),n[i].slope))
            {
                count_bound++;
            }

        }
        return count_bound;
    }


//DELTA_H
    double delta_H(double const alpha, double const J, int num_dif, int num_bonds, double a) const
    {
        //std::cout<<"num_dif "<<num_dif<<"\n num_bonds "<< num_bonds<<"\n";
        double delta_E=(a*J)*(num_bonds-(num_dif*alpha));
        return delta_E<0.0f ? 1.0 : exp(-delta_E);
    }


//BINDING ATTEMPT
    void attempt_binding(double const alpha, double const J, int ind, double &rand)
    {
        int pos = get_pos(ind);
        double rand_size = rand*3;
        int rand_int = rand_size;
        int ori = rand_int-1;
        rand=rand_size-rand_int;

        Interactions interactions;
        interactions.num_bonds=0;
        interactions.num_diffuse=1;


        std::vector<Neighbour> n(lattice.get_neighbors2(pos));

        auto _is_diffuse = [this](const Neighbour &n)
        {
            return is_diffuse(n.position);
        };
        auto _get_compatible_diffuse_neighbours = [&interactions,this,&rand,ori](auto n)
        {
            double rand_size = rand*3;
            int ori2 = rand_size;
            rand=rand_size-ori2;
            ori2--;
            if(is_interaction_allowed(ori,ori2,n.slope))
            {
                interactions.possible_interaction_pos.emplace_back(n.position);
                interactions.orientations.emplace_back(ori2);
                interactions.num_bonds++;
                interactions.num_diffuse++;
                int bound_neigh_of_neigh = count_interacting_neighbors(n.position,ori);
                if(bound_neigh_of_neigh>1)
                {
                    interactions.num_bonds=interactions.num_bonds+(bound_neigh_of_neigh-1);
                }
            }
        };

        auto s= n | std::views::filter(_is_diffuse);
        std::ranges::for_each(s,_get_compatible_diffuse_neighbours);


        auto _is_bound = [this](const Neighbour &n)
        {
            return is_bound(n.position);
        };
        auto _is_allowed = [this,ori](const Neighbour &n)
        {
            return is_interaction_allowed(ori,get_orientation(n.position),n.slope);
        };
        auto _get_affected_particles = [&interactions,this](auto n)
        {
            interactions.orientations.emplace_back(get_orientation(n.position));
            interactions.possible_interaction_pos.emplace_back(n.position);
            interactions.num_bonds++;
        };

        auto s2= n | std::views::filter(_is_bound) | std::views::filter(_is_allowed);
        std::ranges::for_each(s2,_get_affected_particles);



        if(interactions.num_bonds>0)
        {
            binding_attempt++;
            double delta = delta_H(alpha,J,interactions.num_diffuse,interactions.num_bonds, -1);
            //std::cout<<"delta "<< delta << "\n";
            //std::cout<<"\n rand "<<rand;
            //the binding attempt is succesful if rand < delta_E
            if(rand<delta)
            {
//                std::cout<<"binding"<<"\n";
                binding_succ++;
                set_orientation(pos,ori+3);
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],interactions.orientations[i]+3);

                }
            }
        }
    }


//UNBINDING ATTEMPT
    void attempt_unbinding(double const alpha, double const J, int ind,double &rand)
    {
        // get bound particle
        unbinding_attempt++;
        int particle_pos = get_pos(ind);
        int ori = get_orientation(particle_pos);

        Interactions interactions;
        interactions.num_bonds=0;
        interactions.num_diffuse=1;

        int pos = get_pos(ind);
        std::vector<Neighbour> n(lattice.get_neighbors2(pos));

        auto _is_bound = [this](const Neighbour &n)
        {
            return is_bound(n.position);
        };
        auto _is_allowed = [this,ori](const Neighbour &n)
        {
            return is_interaction_allowed(ori,get_orientation(n.position),n.slope);
        };
        auto _get_affected_particles = [&interactions,this](auto n)
        {
            interactions.num_bonds++;
            int bound_neigh_of_neigh = count_interacting_neighbors(n.position,get_orientation(n.position));
            if(bound_neigh_of_neigh==1)
            {
                interactions.possible_interaction_pos.emplace_back(n.position);
                interactions.num_diffuse++;
            }
        };

        auto s= n | std::views::filter(_is_bound) | std::views::filter(_is_allowed);
        std::ranges::for_each(s,_get_affected_particles);



        double delta_E =delta_H(alpha,J,interactions.num_diffuse,interactions.num_bonds, 1);

        if (rand<delta_E)
        {
//            std::cout<<"unbinding"<<"\n";
            set_orientation(particle_pos,1);
            unbinding_succ++;
            if(interactions.possible_interaction_pos.size()>0)
            {
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],1);
                }
            }
        }
    }


    void label(int ind, int i,std::vector<int> &labels,int &num_bonds) const
    {
        labels[ind]=i;

        int pos = get_pos(ind);
        int ori = get_orientation(pos);
        std::vector<Neighbour> n(lattice.get_neighbors2(pos));

        auto _is_bound = [this](const Neighbour &n){return is_bound(n.position);};

        auto _is_allowed = [this,ori](const Neighbour &n){return is_interaction_allowed(ori,get_orientation(n.position),n.slope);};

        auto _is_labelled = [this,&labels](const Neighbour &n)
        {
            auto it= (std::find(begin(positions),end(positions),n.position));
            if(it!= positions.end())
            {
                int inds = it-positions.begin();
                return labels[inds]!=0;
            }
        };


        auto _label = [this,&labels,i,pos,&num_bonds](const Neighbour &n)
        {
            auto it= (std::find(begin(positions),end(positions),n.position));
            if(it!= positions.end())
            {
                int inds = it-positions.begin();
                if(labels[inds]==0)
                {
                    label(inds,i,labels,num_bonds);
                }
            }
        };



        auto s1= n | std::views::filter(_is_bound) | std::views::filter(_is_allowed) | std::views::filter(_is_labelled);
        auto cnt = std::ranges::distance(s1);
        num_bonds=num_bonds+cnt;

        auto s2= n | std::views::filter(_is_bound) | std::views::filter(_is_allowed);
        std::ranges::for_each(s2,_label);
    }

    std::vector<int> orientations_vec() const
    {
        std::vector<int> orientations_vector;

        for(unsigned int index=0; index<positions.size(); index++)
        {
            orientations_vector.emplace_back(grid[get_pos(index)]);
        }

        return orientations_vector;
    }

    void print(std::ofstream &out) const
    {
        std::stringstream buffer;
        std::vector<int> ori_vec=orientations_vec();

//        buffer << t << '\t';
        for(auto const &x: positions)
        {
            buffer << x<< '\t';
        }
        buffer << '\n';
        out << buffer.str();
        std::stringstream buffer1;
        for(auto const &y: ori_vec)
        {
            buffer1 << y<< '\t';
        }
        buffer1 << '\n';
        out << buffer1.str();
    }

    void print_labels(std::ofstream &out,std::vector<int> &labels) const
    {
        std::stringstream buffer;
        for(auto &x:labels)
        {
            buffer << x<< '\t';
        }
        buffer << '\n';
        out << buffer.str();
    }

    void print_grid(std::ofstream &out) const
    {
        std::stringstream buffer;
        for(auto &x:grid)
        {
            buffer << x<< '\t';
        }
        buffer << '\n';
        out << buffer.str();
    }



};



#endif // PARTICLES_H
