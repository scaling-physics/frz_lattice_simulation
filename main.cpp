 # include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <array>
#include <vector>
using namespace std;

// fct takes position of diffuse particles attempts to move one at random,
//if move goes to a lattice point without neighbours then move is accepted
//if move attempts to move next to a particle --> return particle combination for step 3
// RETURN, updated diffuse_pos and grid, create array with attempted binding
void diffuse(vector diffuse_pos,array grid){


}

// find empty hexes without neighbours
// randomely create a particle with k_on/
//RETURN updated grid and diffuse position

//Do I need to check how full the grid is? Set concentration?delete particles randomely?
void create_particle(vector diffuse_pos, vector bound_pos, array grid){
int const k_on = 0.02;


}




//calculate energy change by binding/unbinding a particle
//
void energy_change(vector bound_pos, array grid){
int const alpha=1;
int const J=1;
int no_bound==bound_pos.size()

float nucleation_term = alpha*J*no_bound


}








int main(){
const int MC_steps = 10; // number of Monte Carlo Steps
int MC_counter = 0;

while (MC_counter<MC_Steps){

cout << n <<'\n';
//step 1: Move diffusive particles

diffuse(a,b);
//if move is rejected check if it binds or not




//step 2: Check and create a new particles

create_particle(a,c,b);

//step 3: bind and unbind
//first loop through the array for from step 1 to check if diffusive particles bind/form a new cluster or not





MC_counter++;




}
return 0;
};
