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


struct Interactions
{
    std::vector<int> possible_interaction_pos;
    std::vector<int> orientations;
    int num_diffuse;
    int num_bonds;
    int num_aa_bonds;
    int num_ab_bonds;
};

struct Particle
{
    int pos=0;
    int ori=0;
    int Frz_A_B=0;
    std::vector<int> partners = {};
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
                p-> ori =1;

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
        return (ori_self != ori_other && ((ori_self + (slope%3-1))*(ori_other+(slope%3-1)))!=0 && (flag_other*flag_self==0));
    }
////////////////////////////////////////////////////////////

    //COUNTING THE BOUND PARTICLES OF A GIVEN POSITION
//    std::vector<int> count_interacting_neighbors(const int pos, const int ori) const
//    {
//        std::vector<Neighbour> n(lattice.get_neighbors2(pos));
//        std::vector<int> count_bound = {0,0};
//
//        for(unsigned int i=0; i<n.size(); i++)
//        {
//            if(is_bound(n[i].position) && is_interaction_allowed(ori,get_orientation(n[i].position),get_flag(pos),get_flag(n[i].position),n[i].slope))
//            {
//                if(get_flag(pos)==0 && get_flag(n[i].position==0))
//                    count_bound[0]++;
//                else
//                    count_bound[1]=1;
//            }
//        }
////        assert(count_bound>=0);
//        return count_bound;
//    }


//        COUNTING THE BOUND PARTICLES OF A GIVEN POSITION
    int count_interacting_neighbors(const int pos, const int ori) const
    {
        std::vector<Neighbour> n(lattice.get_neighbors2(pos));
        int count_bound = 0;

        for(unsigned int i=0; i<n.size(); i++)
        {
            if(is_bound(n[i].position) && is_interaction_allowed(ori,get_orientation(n[i].position),get_flag(pos),get_flag(n[i].position),n[i].slope))
            {
//                if(get_flag(pos)==0 && get_flag(n[i].position)==0)
                    count_bound++;
            }
        }
//        assert(count_bound>=0);
        return count_bound;
    }

    int count_interacting_neighbors_ab(const int pos, const int ori) const
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
//        assert(count_bound>=0);
        return count_bound;
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

    int count_flags() const
    {   int count=0;
        for(unsigned int index=0; index<particles.size(); index++)
            {
                if(particles[index]->Frz_A_B>0)
                count++;
            }

        return count;
    }



//DELTA_H
    double delta_H(double const alpha, double const J, int num_dif, int num_aa_bonds, int num_ab_bonds, double a, double beta) const
    {
//        if(num_ab_bonds>0)
//            std::cout<<"num_dif "<<num_dif<<"\t num_aa_bonds "<< num_aa_bonds<<"\t num_ab_bonds"<<num_ab_bonds<<"\n";
        double delta_E=a*(J*(num_aa_bonds-(num_dif*alpha)) + beta*num_ab_bonds);
//        double delta_E=(a*J)*(num_bonds-(num_dif*alpha));
        return delta_E<0.0f ? 1.0 : exp(-delta_E);
    }


//BINDING ATTEMPT
    void attempt_binding(double const alpha, double const J,double const beta, int ind, double &rand)
    {
        int pos = get_pos(ind);
        double rand_size = rand*3;
        int rand_int = rand_size;
        int ori = rand_int-1;
        rand=rand_size-rand_int;

        Interactions interactions;
        interactions.num_bonds=0;
        interactions.num_aa_bonds=0;
        interactions.num_ab_bonds=0;
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

//            std::cout<<get_flag(pos)<<'\n';
            if(is_interaction_allowed(ori,ori2,get_flag(pos),get_flag(n.position),n.slope))
            {
                interactions.possible_interaction_pos.emplace_back(n.position);
                interactions.orientations.emplace_back(ori2);
                interactions.num_diffuse++;
//                interactions.num_bonds++;


                if((get_flag(pos)==0 && get_flag(n.position)==0))              //check if aa is interacting or ab is interacting
                {
                    interactions.num_bonds++;
                }
                if((get_flag(pos)==1 && get_flag(n.position)==0))
                {
                    interactions.num_ab_bonds=1;
                }
                if((get_flag(pos)==0 && get_flag(n.position)==1))
                {
                    interactions.num_ab_bonds++;
                }

                int bound_neigh_of_neigh = count_interacting_neighbors(n.position,ori);
                if(bound_neigh_of_neigh>1)
                {
                    interactions.num_bonds=interactions.num_bonds+(bound_neigh_of_neigh-1);
                }
//                else
//                {
//                    interactions.num_ab_bonds++;
//                }
//
//                int bound_neigh_of_neigh = count_interacting_neighbors(n.position,ori);
//                int bound_neigh_of_neigh2 = count_interacting_neighbors_ab(n.position,ori);
//                if(bound_neigh_of_neigh > 1)
//                {
//                    interactions.num_aa_bonds=interactions.num_aa_bonds+(bound_neigh_of_neigh-1);
//
//                }
//                if(bound_neigh_of_neigh2 > 1)
//                {
//                    interactions.num_ab_bonds = interactions.num_ab_bonds+(bound_neigh_of_neigh2-1);
//                }

            }
         };

        auto s = n | std::views::filter(_is_diffuse);
        std::ranges::for_each(s,_get_compatible_diffuse_neighbours);


        auto _is_bound = [this](const Neighbour &n)
        {
            return is_bound(n.position);
        };
        auto _is_allowed = [this,ori,pos](const Neighbour &n)
        {
            return is_interaction_allowed(ori,get_orientation(n.position),get_flag(pos),get_flag(n.position),n.slope);
        };
        auto _get_affected_particles = [&interactions,this,pos](auto n)
        {
            interactions.orientations.emplace_back(get_orientation(n.position));
            interactions.possible_interaction_pos.emplace_back(n.position);
            interactions.num_bonds++;

            if((get_flag(pos)==0 && get_flag(n.position)==0))              //check if aa is interacting or ab is interacting
            {
                interactions.num_bonds++;
            }
            if((get_flag(pos)==1 && get_flag(n.position)==0))
            {
                interactions.num_ab_bonds=1;
            }
            if((get_flag(pos)==0 && get_flag(n.position)==1))
            {
                interactions.num_ab_bonds++;
            }

        };
//        std::cout<<interactions.num_aa_bonds<<'\n';

        auto s2= n | std::views::filter(_is_bound) | std::views::filter(_is_allowed);
        std::ranges::for_each(s2,_get_affected_particles);



        if((interactions.num_bonds+interactions.num_ab_bonds)>0)
        {
            assert(interactions.possible_interaction_pos.size()>0);
            binding_attempt++;
            double delta = delta_H(alpha,J,interactions.num_diffuse,interactions.num_bonds,interactions.num_ab_bonds, -1,beta);
//            std::cout<<"delta "<< delta << "\n";
            //std::cout<<"\n rand "<<rand;
            //the binding attempt is succesful if rand < delta_E
            if(rand<delta)
            {
//                std::cout<<"binding"<<"\n";
                binding_succ++;
                set_orientation(pos,ori+3);

//                std::cout<<"pos after binding "<<pos<<"\n";
//                std::cout<<"ori after binding "<<get_orientation(pos)<<"\t"<<"Frz_flag"<<get_flag(pos)<<"\n

//                if(get_flag(pos)==1)
//                {
//                    set_orientation(interactions.possible_interaction_pos[0],interactions.orientations[i]+3);
//                }
//                else
//                {
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],interactions.orientations[i]+3);
//                    std::cout<<"pos after binding "<<interactions.possible_interaction_pos[i]<<"\n";
//                    std::cout<<"ori after binding "<<get_orientation(interactions.possible_interaction_pos[i])<<"\t"<<"Frz_flag"<<get_flag(interactions.possible_interaction_pos[i])<<"\n";
                }
//                }
//

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

        if(get_flag(pos)>0)
            std::cout<<get_flag(pos)<<'\t'<<interactions.num_ab_bonds<<'\n';
    }


//UNBINDING ATTEMPT
    void attempt_unbinding(double const alpha, double const J,double const beta, int ind,double &rand)
    {
        // get bound particle
        unbinding_attempt++;

        int pos = get_pos(ind);
        int ori = get_orientation(pos);
        assert(is_bound(pos));
//        std::cout<<"pos pre unbinding "<<pos<<"\n";
//        std::cout<<"ori pre unbinding "<<ori<<"\t"<<"Frz_flag"<<get_flag(pos)<<"\n";
        std::vector<Neighbour> neig(lattice.get_neighbors2(pos));

        Interactions interactions;
        interactions.num_bonds=0;
        interactions.num_aa_bonds=0;
        interactions.num_ab_bonds=0;
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
        auto _get_affected_particles = [&interactions,this,pos](auto n)
        {
//            interactions.num_bonds++;
            int bound_neigh_of_neigh = count_interacting_neighbors(n.position,get_orientation(n.position));
            if(bound_neigh_of_neigh==1)
            {
                interactions.possible_interaction_pos.emplace_back(n.position);
                interactions.num_diffuse++;
            }

            if((get_flag(pos)==0 && get_flag(n.position)==0))              //check if aa is interacting or ab is interacting
            {
                interactions.num_bonds++;
            }
            if((get_flag(pos)==1 && get_flag(n.position)==0))
            {
                interactions.num_ab_bonds=1;
            }
            if((get_flag(pos)==0 && get_flag(n.position)==1))
            {
                interactions.num_ab_bonds++;
            }

//            int bound_neigh_of_neigh = count_interacting_neighbors(n.position,get_orientation(n.position)) + count_interacting_neighbors_ab(n.position,get_orientation(n.position));
//            if(bound_neigh_of_neigh==1)
//            {
//                interactions.possible_interaction_pos.emplace_back(n.position);
//                interactions.num_diffuse++;
//            }
        };

        auto s = n | std::views::filter(_is_bound) | std::views::filter(_is_allowed);
        std::ranges::for_each(s,_get_affected_particles);

        double delta_E = delta_H(alpha,J,interactions.num_diffuse,interactions.num_bonds,interactions.num_ab_bonds, 1,beta);

        if (rand<delta_E)
        {
//            std::cout<<"unbinding"<<"\n";
            set_orientation(pos,1);
            unbinding_succ++;
//            std::cout<<"pos after unbinding "<<pos<<"\n";
//            std::cout<<"ori after unbinding "<<grid[pos].ori<<"\t"<<"Frz_flag"<<get_flag(pos)<<"\n";

            if(interactions.possible_interaction_pos.size()>0)
            {
//
//                if(get_flag(pos)==1)
//                {
//                    rand = unidist(gen);
//                    rand_size = rand*interactions.possible_interaction_pos.size();
//                    int ind=rand_size;
//                    set_orientation(interactions.possible_interaction_pos[ind],interactions.orientations[ind]+3);
//                }
//                else
//                {
                for(unsigned int i=0; i<interactions.possible_interaction_pos.size(); i++)
                {
                    set_orientation(interactions.possible_interaction_pos[i],1);
//                    std::cout<<"pos after unbinding "<<interactions.possible_interaction_pos[i]<<"\n";
//                    std::cout<<"ori after unbinding "<<grid[interactions.possible_interaction_pos[i]].ori<<"\t"<<"Frz_flag"<<get_flag(interactions.possible_interaction_pos[i])<<"\n";
                }
//                }
            }
        }
//        if(get_flag(pos)>0)
//            std::cout<<get_flag(pos)<<'\t'<<interactions.num_ab_bonds<<'\n';
//        else{std::cout<<"no unbinding"<<"\n";}
    }


// Binding FrzB acceptance rate 1
    void binding_FrzB(int ind, int &FrzB_num)
    {
//        assert (rand<1 && rand>0);
        if (FrzB_num>0)
        {
        int pos=get_pos(ind);
        int frz = get_flag(pos);

        if(frz>0)
        {
            return;
        }
        else
        {
            FrzB_num--;
            set_frz(ind,1);
        }
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
//            std::cout<<x<<'\t';
            buffer <<x<< '\t';
        }
//        std::cout<<'\n';
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
