# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <array>
#include <vector>
using namespace std;

// fct takes position of diffuse particles attempts to move them,
//if move goes to a lattice point without neighbours then move is accepted
//if move attempts to move to next to a particle --> return particle combination for step 3
// RETURN, updated diffuse_pos and grid, create array holding attempted binding
int diffuse(vector diffuse_pos,array grid){

return 0;
}

// find empty hexes
// randomely create particles with k_on
//RETURN updated grid and diffuse position
int create_particle(vector diffuse_pos, vector bound_pos, array grid){

return 0;
}

//calculate energy change by binding/unbinding a particle
//
int energy_change(vector bound_pos, array grid){

return 0;
}

int main(){
const int N = 10; // number of Monte Carlo Steps
int n = 0;

while (n<N){

cout << n <<'\n';
//step 1: Move diffusive particles

diffuse(a,b);

//step 2: Check and create new particles


//step 3: bind and unbind

n++;




}
return 0;
};
