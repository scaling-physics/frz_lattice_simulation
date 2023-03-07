#ifndef PARTICLES_H
#define PARTICLES_H

# include <fstream>
# include <memory>
# include <cmath>
# include <iostream>	// cout, etc.
# include <sstream>      // std::stringstrea
# include <vector>
# include <array>
# include <algorithm>
# include "lattice.h"

int unbinding_attempt=0;
int unbinding_succ=0;
int binding_attempt=0;
int binding_succ=0;
int diffuse_attempt=0;
int diffuse_succ=0;

double const  k_off=5*pow(10,-4);
double const k_on=k_off/55;

double density;
double titration_concentration_frzb;



unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);
std::default_random_engine generator (seed);

void print_container(const std::vector<int>& c)
{
    for (auto &i : c)
    {
        std::cout << i << " ";
    }
    std::cout << '\n';
}

const bool Frz_Interaction_Lookuptable[6*32]=
{
    true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,
    true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,
    true,true,true,true,false,false,false,false,true,true,true,true,false,false,false,false,
    true,true,false,false,true,true,false,false,true,true,false,false,true,true,false,false,
    true,true,false,false,true,true,false,false,true,true,false,false,true,true,false,false,
    true,true,true,true,false,false,false,false,true,true,true,true,false,false,false,false,
    true,true,false,false,true,true,false,false,true,true,false,false,true,true,false,false,
    true,true,true,true,false,false,false,false,true,true,true,true,false,false,false,false,
    true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false,
    true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,
    true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,
    true,true,true,true,true,true,true,true,false,false,false,false,false,false,false,false
};
//00-15:   slope 0 ori -1,0 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//16-31:   slope 0 ori 0,-1 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//32-47:   slope 1 ori -1,1 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//48-63:   slope 1 ori 1,-1 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//64-79:   slope 2 ori 0,1  Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//80-95:   slope 2 ori 1,0  Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//96-111:  slope 3 ori -1,0 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//112-127: slope 3 ori 0,-1 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//128-143: slope 4 ori -1,1 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//144-159: slope 4 ori 1,-1 Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//160-175: slope 5 ori 0,1  Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...
//176-191: slope 5 ori 1,0  Frz_flag 0,0 0,1 0,2 0,3 1,0 1,1 1,2...

const bool Frz_occupied_site[6*12]=
{
    false,false,false,false,true,false,true,false,false,false,false,false,
    true,true,false,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,true,true,false,false,
    false,false,false,false,true,true,false,false,false,false,false,false,
    true,false,true,false,false,false,false,false,false,false,false,false,
    false,false,false,false,false,false,false,false,true,false,true,false
};

//00-11 slope 0 possible directions for Frz B to point ori order -1,0,1
//12-23 slope 1
//24-35 slope 2
//36-47 slope 3
//48-59 slope 4
//60-72 slope 5


struct Interactions
{
    std::vector<int> possible_interaction_pos;
    std::vector<int> orientations;
    int num_diffuse;
    int num_bonds;
};

struct Particle
{
    int pos=0;
    int ori=0;
    int Frz_A_B=0;
};




class Particles

{
private:

public:
    std::vector<std::weak_ptr<Particle> > grid1;
    std::vector<std::shared_ptr<Particle> > particles;

    Lattice &lattice;

//////////// Constructor

    Particles(Lattice &lattice):lattice{lattice}
    {
        int initial_num = Nxy*density;
        particles.resize(initial_num);
        grid1.resize(Nxy);

//////// Read in from file, stable initital configuration
//        std::ifstream read_pos("read_in_pos_side.txt");
//        std::ifstream read_ori("read_in_ori_side.txt");
////        std::ifstream read_frz("read_in_frz.txt");
//        int temp_pos,temp_ori;
//        int ind_pos=0;
//        int ind_ori=0;
////
//        for (int i=0; i<initial_num; i++)
//        {
//            auto p = std::make_shared<Particle>();
//            particles[i] = p;
//        }
//
//        while (read_pos>>temp_pos)
//        {
//            particles[ind_pos]->pos = temp_pos;
//            ind_pos++;
//        }
//        while (read_ori>>temp_ori)
//        {
//            particles[ind_ori]->ori = temp_ori;
//            ind_ori++;
//        }
//
//        for (int i=0; i<initial_num; i++)
//        {
//            auto p = particles[i];
//            int pos = p->pos;
//            std::weak_ptr<Particle> pw(p);
//            grid1[pos]= pw;
//        }
////          int temp_frz
////          int ind_frz=0;
//////        while (read_frz>>temp_frz)
//////        {
//////            particles[ind_frz]->Frz_A_B = temp_frz;
//////        }


//// Random initial configuration
        for (int i=0; i<initial_num; i++)
        {
            double random=unidist(gen);
            double random_size = random*Nxy;
            int pos=random_size;
            random=random_size-pos;
//            double titration_concentration_frzb=0.5;
            if(grid1[pos].expired())
            {
                auto p = std::make_shared<Particle>();
                p-> pos = pos;
                p->ori =1;

                particles[i] = p;
                std::weak_ptr<Particle> pw(p);
                grid1[pos]= pw;
            }
            else
            {
                i--;
            }

        }
    }

//////////// GETTERS, SETTERS //////////////
    inline int get_pos(const int ind) const
    {
        auto p = particles[ind];
        int pos= p->pos;
        return pos;
    }

    inline bool is_free(const int pos) const
    {
        return grid1[pos].expired();
    }

    inline bool is_diffuse(const int pos) const
    {
        int ori=0;
        if(!grid1[pos].expired())
        {
            auto p = std::shared_ptr<Particle> (grid1[pos].lock());

            ori = p-> ori;
        }
        return ori==1;
    }
    inline bool is_bound(const int pos) const
    {
        int ori=0;
        if(!grid1[pos].expired())
        {
            auto p = std::shared_ptr<Particle> (grid1[pos].lock());

            ori = p-> ori;
        }
        return ori>1;
    }
    inline int get_orientation(const int pos) const
    {
        int ori=0;
        assert(!grid1[pos].expired());
        if(!grid1[pos].expired())
        {
            auto p = std::shared_ptr<Particle> (grid1[pos].lock());

            ori = p-> ori;
        }
        return ori-3;
    }

    inline int get_flag(const int pos) const
    {
        if(!grid1[pos].expired())
        {
            auto p = std::shared_ptr<Particle> (grid1[pos].lock());

            int flag = p-> Frz_A_B;
            return flag;
        }
    }

    inline std::vector<int> get_free_flag_sites(const int pos)
    {
        std::vector<int> free_sites;
        int ori = get_orientation(pos);
        int frz = get_flag(pos);
        std::vector<Neighbour> n(lattice.get_neighbors2(pos));

        if(ori==0)
        {
            auto it1=(std::find_if(begin(n),end(n), [](Neighbour m)
            {
                return m.slope==0;
            }));
            if(it1!=n.end())
            {
                int n_ind = it1-n.begin();
                if(grid1[n[n_ind].position].expired())
                {
                    int pos_neighbor1=n[n_ind].position;
                    free_sites.emplace_back(2);
                }
            }

            auto it2=(std::find_if(begin(n),end(n), [](Neighbour m)
            {
                return m.slope==3;
            }));
            if(it2!=n.end())
            {
                int n_ind = it2-n.begin();
                if(grid1[n[n_ind].position].expired())
                {
                    int pos_neighbor2=n[n_ind].position;
                    free_sites.emplace_back(1);
                }
            }
        }
        else if (ori==1)
        {
            auto it1=(std::find_if(begin(n),end(n), [](Neighbour m)
            {
                return m.slope==2;
            }));
            if(it1!=n.end())
            {
                int n_ind = it1-n.begin();
                if(grid1[n[n_ind].position].expired())
                {
                    int pos_neighbor1=n[n_ind].position;
                    free_sites.emplace_back(1);
                }
            }

            auto it2=(std::find_if(begin(n),end(n), [](Neighbour m)
            {
                return m.slope==5;
            }));
            if(it2!=n.end())
            {
                int n_ind = it2-n.begin();
                if(grid1[n[n_ind].position].expired())
                {
                    int pos_neighbor2=n[n_ind].position;
                    free_sites.emplace_back(2);
                }
            }
        }

        else if(ori==-1)
        {
            auto it1=(std::find_if(begin(n),end(n), [](Neighbour m)
            {
                return m.slope==1;
            }));
            if(it1!=n.end())
            {
                int n_ind = it1-n.begin();
                if(grid1[n[n_ind].position].expired())
                {
                    int pos_neighbor1=n[n_ind].position;
                    free_sites.emplace_back(2);
                }
            }

            auto it2=(std::find_if(begin(n),end(n), [](Neighbour m)
            {
                return m.slope==4;
            }));
            if(it2!=n.end())
            {
                int n_ind = it2-n.begin();
                if(grid1[n[n_ind].position].expired())
                {
                    int pos_neighbor2=n[n_ind].position;
                    free_sites.emplace_back(1);
                }
            }
        }


        return free_sites;
    }

    inline void set_orientation(const int pos, const int ori)
    {
        if(!grid1[pos].expired())
        {
            auto p = std::shared_ptr<Particle> (grid1[pos].lock());

            p-> ori = ori;
        }
    }
    inline void set_pos(const int old_pos,const int new_pos, const int ind)
    {
        assert(particles[ind]->pos == old_pos);
        particles[ind] -> pos = new_pos;
        assert(particles[ind]->pos == new_pos);
        grid1[new_pos].swap(grid1[old_pos]);
    }

    inline void set_frz(const int ind, const int frz)
    {
        particles[ind] -> Frz_A_B = frz;
//        if(!grid1[pos].expired())
//        {
//            auto p = std::shared_ptr<Particle> (grid1[pos].lock());
//
//            p-> Frz_A_B = frz;
//        }
    }
    inline bool is_interaction_allowed(const int ori_self, const int ori_other, const int flag_self, const int flag_other, const int slope) const
    {
        int order_of_flags;
        if(ori_self<ori_other)
        {
            order_of_flags=0;
        }
        else
        {
            order_of_flags=16;
        }
        int key = slope*32 + order_of_flags + flag_self*4 + flag_other;
//        bool b = (ori_self != ori_other && ((ori_self + (slope%3-1))*(ori_other+(slope%3-1)))!=0 && Frz_Interaction_Lookuptable[key]);
        return (ori_self != ori_other && ((ori_self + (slope%3-1))*(ori_other+(slope%3-1)))!=0 && Frz_Interaction_Lookuptable[key]);
//        return (ori_1 != ori_2 && ((ori_1 + slope)*(ori_2+slope))!=0);
    }
////////////////////////////////////////////////////////////

//// needs fixing
//    inline bool is_occupied_by_FrzB(const int pos) const
//    {
//        std::vector<Neighbour> n(lattice.get_neighbors2(pos));
//
//        for(unsigned int i=0; i<n.size(); i++)
//        {
//            if(is_bound(n[i].position))
//            {
//                int key=n[i].slope*12+(get_orientation(n[i].position)+1)*4+get_flag(n[i].position);
//                bool b= Frz_occupied_site[key];
//                if(b)
//                {
//                    return true;
//                }
//            }
//        }
//        return false;
//    }
    inline bool is_frz_interaction_allowed(const int flag_1, const int flag_2, const int slope) const
    {
        return (flag_1!=1 && flag_2!=1);
    }

    //COUNTING THE BOUND PARTICLES OF A GIVEN POSITION
    int count_interacting_neighbors(const int pos, const int ori) const
    {
        std::vector<Neighbour> n(lattice.get_neighbors2(pos));
        int count_bound = 0;

        for(unsigned int i=0; i<n.size(); i++)
        {
            if(is_bound(n[i].position) && is_interaction_allowed(ori,get_orientation(n[i].position),get_flag(pos),get_flag(n[i].position),n[i].slope))
            {
                count_bound++;
            }

        }
        assert(count_bound>=0);
        return count_bound;
    }

////CREATION ATTEMPT
    void attempt_creation(double &rand)
    {
        int initial_num = Nxy*density;
        int free_hexes=(Nx*Ny)-particles.size();
        double k_on_1 = k_on*(initial_num-particles.size());
        std::binomial_distribution<int> bnom(free_hexes,k_on_1);
        double random=unidist(gen);

        int new_particles_created = bnom(generator);
        if(particles.size()+new_particles_created>initial_num)
        {
            new_particles_created=new_particles_created+(initial_num-particles.size()-new_particles_created);
        }
        for ( int i=0;i<new_particles_created;i++)
        {
            double random_size = random*Nxy;
            int pos=random_size;
            random=random_size-pos;
            if(grid1[pos].expired())
            {
                auto p = std::make_shared<Particle>();
                p-> pos = pos;
                p->ori =1;

                particles.emplace_back(p);
                std::weak_ptr<Particle> pw(p);
                grid1[pos]= pw;
            }
            else
            {
                i--;
            }

        }
        rand=random;
    }
//
//
////DESTRUCTION ATTEMPT
    void attempt_destruction(double &rand)
    {
        for(unsigned int ind=0; ind<particles.size(); ind++)
        {
            if(is_diffuse(get_pos(ind)))
            {
                rand = unidist(gen);
                if(rand<k_off)
                {
                    int pos = get_pos(ind);
                    assert(!grid1[pos].expired());
                    auto it = particles.begin()+ind;
                    particles.erase(it);
                    assert(grid1[pos].expired());
                }
            }

        }

    }


//CELL GROWTH
    void cell_growth(double &rand,int &Nx)
    {
        double rand_size = Nx*rand;
        int column_insert = rand_size;
        rand = rand_size-column_insert;
        Nx++;
        std::vector<std::weak_ptr<Particle> > column;
        column.resize(Ny);
        auto it = grid1.begin()+column_insert*Ny;

        grid1.insert(it, column.begin(),column.end());
        std::cout<<grid1.size();

        for(int i=(Nx*Ny-1); i>=Ny*column_insert; i--)
        {
            if(!grid1[i].expired())
            {
                auto p = std::shared_ptr<Particle> (grid1[i].lock());

                int pos = p->pos;
                int new_pos = pos+Ny;
                p->pos = new_pos;
            }
        }
        for(unsigned int i=(Nx*Ny-1); i>=Ny*column_insert; i--)
        {
            if(!grid1[i].expired() && is_bound(i))
            {
                auto p = std::shared_ptr<Particle> (grid1[i].lock());

                int pos = p->pos;
                int ori = p-> ori;
                int is_still_bound = count_interacting_neighbors(pos, ori);
                if(is_still_bound>=1)
                {
                    return;
                }
                else
                {
                    set_orientation(pos,1);
                }
            }
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
//        if(is_free(new_pos) && !is_occupied_by_FrzB(new_pos))// accept move
        if(is_free(new_pos))
        {
            set_pos(particle_pos,new_pos, ind);
            diffuse_succ++;
        }
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
        auto _get_compatible_diffuse_neighbours = [&interactions,this,&rand,ori,pos](auto n)
        {
            double rand_size = rand*3;
            int ori2 = rand_size;
            rand=rand_size-ori2;
            ori2--;
            if(is_interaction_allowed(ori,ori2,get_flag(pos),get_flag(n.position),n.slope))
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
        auto _is_allowed = [this,ori,pos](const Neighbour &n)
        {
            return is_interaction_allowed(ori,get_orientation(n.position),get_flag(pos),get_flag(n.position),n.slope);
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
            assert(interactions.possible_interaction_pos.size()>0);
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
//                std::cout<<"pos after binding "<<pos<<"\n";
//                std::cout<<"ori after binding "<<get_orientation(pos)<<"\t"<<"Frz_flag"<<get_flag(pos)<<"\n";
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],interactions.orientations[i]+3);
//                    std::cout<<"pos after binding "<<interactions.possible_interaction_pos[i]<<"\n";
//                    std::cout<<"ori after binding "<<get_orientation(interactions.possible_interaction_pos[i])<<"\t"<<"Frz_flag"<<get_flag(interactions.possible_interaction_pos[i])<<"\n";
                }

                std::vector<Neighbour> neig(lattice.get_neighbors2(pos));
                int test1=0;

                for(unsigned int i=0; i<neig.size(); i++)
                {
                    if(is_bound(neig[i].position))
                    {
                        test1++;
                    }
                }

                int test = count_interacting_neighbors(pos,ori);
                assert(test>0);
                assert(test1>0);
            }
        }
    }


//UNBINDING ATTEMPT
    void attempt_unbinding(double const alpha, double const J, int ind,double &rand)
    {
        // get bound particle
        unbinding_attempt++;

        int pos = get_pos(ind);
        int ori = get_orientation(pos);
        assert(is_bound(pos));
//        std::cout<<"pos pre unbinding "<<pos<<"\n";
//        std::cout<<"ori pre unbinding "<<ori<<"\t"<<"Frz_flag"<<get_flag(pos)<<"\n";
        std::vector<Neighbour> neig(lattice.get_neighbors2(pos));
//////        int test1=0;
//////
//////        for(unsigned int i=0; i<neig.size(); i++)
//////        {
//////            if(is_bound(neig[i].position)&& is_interaction_allowed(ori,get_orientation(neig[i].position),get_flag(pos),get_flag(neig[i].position),neig[i].slope))
//////            {
////////                        std::cout<<"pos pre unbinding "<<neig[i].position<<"\n";
////////                        std::cout<<"ori pre unbinding "<<get_orientation(neig[i].position)<<"\t"<<"Frz_flag"<<get_flag(neig[i].position)<<"\n";
//////                test1++;
//////            }
//////        }
//////        assert(test1>0);
//////        int test = count_interacting_neighbors(pos,ori);
//////        assert(test>0);

        Interactions interactions;
        interactions.num_bonds=0;
        interactions.num_diffuse=1;


        std::vector<Neighbour> n(lattice.get_neighbors2(pos));

        auto _is_bound = [this](const Neighbour &n)
        {
            return is_bound(n.position);
        };
        auto _is_allowed = [this,ori,pos](const Neighbour &n)
        {
            return is_interaction_allowed(ori,get_orientation(n.position),get_flag(pos),get_flag(n.position),n.slope);
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
            set_orientation(pos,1);
            unbinding_succ++;
//            std::cout<<"pos after unbinding "<<pos<<"\n";
//            std::cout<<"ori after unbinding "<<grid[pos].ori<<"\t"<<"Frz_flag"<<get_flag(pos)<<"\n";
            if(interactions.possible_interaction_pos.size()>0)
            {
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],1);
//                    std::cout<<"pos after unbinding "<<interactions.possible_interaction_pos[i]<<"\n";
//                    std::cout<<"ori after unbinding "<<grid[interactions.possible_interaction_pos[i]].ori<<"\t"<<"Frz_flag"<<get_flag(interactions.possible_interaction_pos[i])<<"\n";
                }
            }
        }
//        else{std::cout<<"no unbinding"<<"\n";}
    }


// Binding FrzB acceptance rate 1
    void binding_FrzB(int ind, int &FrzB_num, double &rand)
    {
//        assert (rand<1 && rand>0);
        if (FrzB_num>0)
        {


        int pos=get_pos(ind);
        int frz = get_flag(pos);
//        std::cout<<"FrzB_num \t"<<FrzB_num<<"\t"<<"pos frz \t"<< pos<<"\t"<<frz<<"\t"<<"rand \t"<<rand<<"\t";
//        double rate_acceptance = FrzB_num/(double)420;//(Nxy-particles.size());
        double rate_acceptance = 1;//(Nxy-particles.size());

//        std::cout<<"rate_acc \t"<<rate_acceptance<<"\t";

        if(is_diffuse(pos))
        {
            if(frz==3)
            {
                return;
            }
            else if(frz==0 && rand<rate_acceptance)
            {
                int new_frz = unidist(gen)*2;
                FrzB_num--;
                set_frz(ind,new_frz+1);
//                std::cout<<1<<"\t";
            }
            else if(frz>0 && frz<3 && rand<rate_acceptance)
            {
                FrzB_num--;
                set_frz(ind,3);
//                std::cout<<2<<"\t";
            }
        }
        else if(is_bound(pos))
        {
            if(frz==3)
            {
                return;
            }
            else if(frz==0 && rand<rate_acceptance )
            {
                std::vector<int> free_sites(get_free_flag_sites(pos));
                if(free_sites.size()==1 && rand<rate_acceptance)
                {
                    set_frz(ind, free_sites[0]);
                    FrzB_num--;
//                    std::cout<<4<<"\t";
                }
                else if(free_sites.size()==2 && rand<rate_acceptance)
                {
                    int new_frz = unidist(gen)*2;
                    set_frz(ind,new_frz+1);
                    FrzB_num--;
//                    std::cout<<5<<"\t";
                }
            }
            else if(frz>0 && frz<3 && rand<rate_acceptance)
            {
                std::vector<int> free_sites(get_free_flag_sites(pos));
                if(free_sites.size()==2 && rand<rate_acceptance)
                {
                    set_frz(ind,3);
                    FrzB_num--;
//                    std::cout<<6<<"\t";FrzB_{run}_image_J_{J}_alpha_{alpha}_FrzB_{FrzB}_off_{rate}')
                }
            }
        }
        }
//        std::cout<<FrzB_num<<"\n";
    }
// Unbinding FrzB acceptance rate proportional to e^(-deltaE)=const
    void unbinding_FrzB(const int ind,int &FrzB_num, const double &rate, double &rand )
    {
//        return;
        int pos=get_pos(ind);
        int frz = get_flag(pos);
        if(frz==0)
        {
            return;
        }

        else if(frz>0 && frz<3 && rand<rate)
        {
            FrzB_num++;
            set_frz(ind,0);
            //options are to go from one to no FrzB
        }
        else if(frz==3 && (rate*2)>rand)
        {
            int new_frz = unidist(gen)*2;
            FrzB_num++;
            set_frz(ind,new_frz+1);
            //options are no FrzB = 0, or one FrzB = 1 or 2
        }
    }


// labeling a cluster
    void label(int ind, int i,std::vector<int> &labels,int &num_bonds) const
    {
        labels[ind]=i;

        int pos = get_pos(ind);
        int ori = get_orientation(pos);
        std::vector<Neighbour> n(lattice.get_neighbors2(pos));
//        std::cout<<n.size();
//        int test_all=0;
//        for(unsigned int test_i=0; test_i<n.size(); test_i++)
//        {
//            if(is_bound(n[test_i].position) && is_interaction_allowed(ori, get_orientation(n[test_i].position),get_flag(pos),get_flag(n[test_i].position),n[test_i].slope))
//            {
//                test_all++;
//            }
//        }
//        assert(test_all>0);

        auto _is_bound = [this](const Neighbour &n)
        {
            return is_bound(n.position);
        };

        auto _is_allowed = [this,ori,pos](const Neighbour &n)
        {
            return is_interaction_allowed(ori,get_orientation(n.position),get_flag(pos),get_flag(n.position),n.slope);
        };

        auto _is_labelled = [this,&labels](const Neighbour &n)
        {
            auto it=(std::find_if(begin(particles),end(particles), [n](std::shared_ptr<Particle> q)
            {
                return q->pos == n.position;
            }));
            if(it!= particles.end())
            {
                int inds = it-particles.begin();
                return labels[inds]!=0;
            }
            std::cout << "error";
            return false;
        };


        auto _label = [this,&labels,i,pos,&num_bonds](const Neighbour &n)
        {
            auto it=(std::find_if(begin(particles),end(particles), [n](std::shared_ptr<Particle> q)
            {
                return q->pos == n.position;
            }));
            if(it!= particles.end())
            {
                int inds = it-particles.begin();
                if(labels[inds]==0)
                {
                    label(inds,i,labels,num_bonds);
                }
            }
        };



        auto s1= n | std::views::filter(_is_bound) | std::views::filter(_is_allowed)| std::views::filter(_is_labelled);
        auto cnt = std::ranges::distance(s1);
        num_bonds=num_bonds+cnt;

        auto s2= n | std::views::filter(_is_bound) | std::views::filter(_is_allowed);
        std::ranges::for_each(s2,_label);
    }

    std::vector<int> orientations_vec() const
    {
        std::vector<int> orientations_vector;

        for(unsigned int index=0; index<particles.size(); index++)
        {
            orientations_vector.emplace_back(particles[index]->ori);
        }

        return orientations_vector;
    }

    void print(std::ofstream &out) const
    {
        std::stringstream buffer;
        std::vector<int> ori_vec=orientations_vec();

//        buffer << t << '\t';
        for(auto const &x: particles)
        {
            buffer << x-> pos << '\t';
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


    std::vector<int> Frz_vec() const
    {
        std::vector<int> Frz_vector;

        for(unsigned int index=0; index<particles.size(); index++)
        {
            Frz_vector.emplace_back(particles[index]->Frz_A_B);
        }

        return Frz_vector;
    }

    void print_Frz(std::ofstream &out) const
    {
        std::stringstream buffer;
        std::vector<int> Frz_vector=Frz_vec();

        for(auto &x:Frz_vector)
        {
            buffer << x<< '\t';
        }
        buffer << '\n';
        out << buffer.str();
    }

//    void print_grid(std::ofstream &out) const
//    {
//        std::stringstream buffer;
//        for(auto &x:grid)
//        {
//            buffer << x.ori<< '\t';
//        }
//        buffer << '\n';
//        out << buffer.str();
//    }



};



#endif // PARTICLES_H
