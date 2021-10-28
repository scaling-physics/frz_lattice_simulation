# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
#include <array>
#include <algorithm>
# include <chrono>
# include <random>
    const int Nx=4;
    const int Ny=4;
//generate random number
unsigned long int seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 gen(seed);
std::uniform_real_distribution<double> unidist(0.0,1.0);


using namespace std;
vector<int> v[10];

// Function for inserting elements
// in array of vectors
void insertionInArrayOfVectors()
{

    for (int i = 0; i < 5; i++) {

        // Inserting elements at every
        // row i using push_back()
        // function in vector
        for (int j = i + 1; j < 5; j++) {
            v[i].push_back(j);
        }
    }
}

// Function to print elements in array
// of vectors
void printElements()
{

    // Traversing of vectors v to print
    // elements stored in it
    for (int i = 0; i < 5; i++) {

        cout << "Elements at index "
             << i << ": ";

        // Displaying element at each column,
        // begin() is the starting iterator,
        // end() is the ending iterator
        for (auto it = v[i].begin();
             it != v[i].end(); it++) {

            // (*it) is used to get the
            // value at iterator is
            // pointing
            cout << *it << ' ';
        }
        cout << endl;
    }
}

// Function to illustrate array
// of vectors
void arrayOfVectors()
{
    // Inserting elements in array
    // of vectors
    insertionInArrayOfVectors();

    // Print elements stored in array
    // of vectors
    printElements();
}

std::vector<int> get_neighbors_a(int const &pos, int const &L)
        {std::vector<int> neighbors;

        //corners
        std::array<int,4> corners={0,L-1, (L-1)*L,L*L-1};
        //edges
        std::vector<int> edge_top(L-2);
            for (int i=0;i<L-2;i++){edge_top[i]=i+1;}
        std::vector<int> edge_bottom(L-2);
            for (int i=0;i<L-2;i++){edge_bottom[i]=(L-1)*L+i+1;}
        std::vector<int> edge_left(L-2);
            for (int i=0;i<L-2;i++){edge_left[i]=(i+1)*L;}
        std::vector<int> edge_right(L-2);
            for (int i=0;i<L-2;i++){edge_right[i]=(i+2)*L-1;}

//        for (int j= 0; j<2; j++)
//        {cout<<edge_left[j] << '\n';}
//        for (int j= 0; j<2; j++)
//        {cout<<edge_right[j] << '\n';}


        if(pos==corners[0]){
            neighbors={pos+1,pos+L,pos+((L-1)*L),pos+((L-1)*L+1)};
            cout <<"corner 0 being executed";
            }

        else if(pos==corners[1]){
            neighbors={pos-1,pos+(L-1),pos+L,pos+((L-1)*L)};
            cout <<"corner 1 being executed";
            }
        else if(pos==corners[2]){
            neighbors={pos-(L-1)*L,pos-L,pos-(L+1),pos+1};
            cout <<"corner 2 being executed";
            }
        else if(pos==corners[3]){
            neighbors={pos-(L-1)*L-1,pos-(L-1)*L,pos-L,pos-1};
            cout <<"corner 3 being executed";
            }

        else if(std::find(edge_top.begin(), edge_top.end(),pos) !=edge_top.end()){
            if(pos%2==0)
            {neighbors={pos-1, pos+1, pos+L, pos+(L-1)*L-1, pos+(L-1)*L, pos+(L-1)*L+1};}
            else
            {neighbors={pos-1, pos+1, pos+L-1, pos+L, pos+L+1, pos+(L-1)*L};}
            }

        else if(std::find(edge_bottom.begin(), edge_bottom.end(),pos) !=edge_bottom.end()){
            if(pos%2==0)
            {neighbors={pos-(L-1)*L, pos-L-1, pos-L, pos-L+1, pos-1, pos+1};}
            else
            {neighbors={pos-(L-1)*L-1, pos-(L-1)*L, pos-(L-1)*L+1, pos-L, pos-1, pos+1};}
            }

        else if(std::find(edge_left.begin(), edge_left.end(),pos) !=edge_left.end()){
            neighbors={pos-L, pos-L+1, pos+1, pos+L};
            }

        else if(std::find(edge_right.begin(), edge_right.end(),pos) !=edge_right.end()){
            neighbors={pos-L, pos-1, pos+L-1,pos+L};
            }

        else{
            if(pos%2==0)
            {neighbors={pos-L-1, pos-L, pos-L+1, pos-1, pos+1, pos+L};}
            else
            {neighbors={pos-L, pos-1, pos+1, pos+L-1, pos+L, pos+L+1};}
            }
            return neighbors;
        }

std::array<int,2> get_index(int const pos)
    {
        //transform 1d array into 2d notation
        int x=pos%Nx;
        int y=pos/Nx;

        std::array<int, 2> coord{x,y};


        return coord;
    }

    int get_single_index(int x, int y)
    {
        int index = y*Ny+x;
        return index;

    }

    bool is_in_lattice(std::array<int,2> coord)
    {
        int x=coord[0];
        int y= coord[1];
        return ( (unsigned int)x<Nx && (unsigned int)y<Ny );//negative values are looped around to big numbers by the cast to unsigned. In this way, only one inequality is required.
    }

    std::vector<int> get_neighbors(int pos)
    {
        std::vector<int> neighbors;
        std::array<int,2> coord;
        coord = get_index(pos);

        int x=coord[0];
        int y= coord[1];



        if(y%2==0)
        {
            if(x-1>=0)
            {
                neighbors.push_back(get_single_index(x-1,y));
                neighbors.push_back(get_single_index(x-1,(y+1)%Ny));
                neighbors.push_back(get_single_index(x-1,(Ny+y-1)%Ny));
            }


        if(x+1<Nx)
        {
            neighbors.push_back(get_single_index(x+1,y));
        }

        neighbors.push_back(get_single_index(x,(y+1)%Ny));
        neighbors.push_back(get_single_index(x,(Ny+y-1)%Ny));



    }

    else
    {
        if(x-1>=0)
        {
            neighbors.push_back(get_single_index(x-1,y));
        }

        if(x+1<Nx)
        {
            neighbors.push_back(get_single_index(x+1,y));
            neighbors.push_back(get_single_index(x+1,(y+1)%Ny));
            neighbors.push_back(get_single_index(x+1,(Ny+y-1)%Ny));
        }

        neighbors.push_back(get_single_index(x,(y+1)%Ny));
        neighbors.push_back(get_single_index(x,(Ny+y-1)%Ny));



    }
    sort(neighbors.begin(),neighbors.end());
    return neighbors;
}





// Driver code
int main()
{
    //arrayOfVectors();
    vector<int> neighbors;
    vector<int> neighbors_a;
    array<int,2> coord;

//    for (auto iter = x.begin(); iter != x.end(); ++iter) {
//        std::cout << *iter << ' ';
//    }
    const int L=10;
    int Lsq = L*L;
    //cout << Lsq;
    int pos=5;

    neighbors_a=get_neighbors_a(pos,L);
    for (auto iter = neighbors_a.begin(); iter != neighbors_a.end(); ++iter) {
        std::cout << *iter << "\n"<<' ';
        }

    neighbors= get_neighbors(pos);
    cout<<"neighbor"<<"\n";
    for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        std::cout << *iter << "\n"<<' ';
        }
    sort(neighbors.begin(),neighbors.end());
    cout<<"neighbor sorted"<<"\n";
    for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        std::cout << *iter << "\n"<<' ';
        }
//    cout << unidist(gen)*Lsq;
    return 0;
}
