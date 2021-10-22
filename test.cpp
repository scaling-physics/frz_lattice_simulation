# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
#include <array>
#include <set>
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
std::set<int> edge_top;
std::set<int> edge_bottom;
std::set<int> edge_left;
std::set<int> edge_right;

for (int j= 0; j<4; j++)
{cout<<corners[j] << '\n';}


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

else if(edge_top.find(pos) !=edge_top.end()){
neighbors={0,0,0,0,0,0};
}

else if(edge_bottom.find(pos) !=edge_bottom.end()){
neighbors={0,0,0,0,0,0};
}

else if(edge_left.find(pos) !=edge_left.end()){
neighbors={0,0,0,0,0,0};
}

else if(edge_right.find(pos) !=edge_right.end()){
neighbors={0,0,0,0,0,0};
}

else{
neighbors={0,0,0,0,0,0};
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
    int pos=15;

    neighbors=get_neighbors(pos,L);
    for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        std::cout << *iter << ' ';
        }



    return 0;
}
