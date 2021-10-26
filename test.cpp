# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
#include <array>
#include <algorithm>
#include <set>
# include <chrono>
# include <random>

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

std::vector<int> get_neighbors(int const &pos, int const &L)
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

        for (int j= 0; j<2; j++)
        {cout<<edge_left[j] << '\n';}
        for (int j= 0; j<2; j++)
        {cout<<edge_right[j] << '\n';}


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
// Driver code
int main()
{
    //arrayOfVectors();
    array<int,5> x{0};
    vector<int> neighbors;

//    for (auto iter = x.begin(); iter != x.end(); ++iter) {
//        std::cout << *iter << ' ';
//    }
    const int L=4;
    int Lsq = L*L;
    //cout << Lsq;
    int pos=13;

    neighbors=get_neighbors(pos,L);
    for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        std::cout << *iter << ' ';
        }


    cout << unidist(gen)*Lsq;
    return 0;
}
