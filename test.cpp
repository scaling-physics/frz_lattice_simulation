# include <fstream>
# include <cmath>
# include <iostream>	// cout, etc.
#include <sstream>      // std::stringstrea
#include <vector>
using namespace std;
vector<int> neighbors[L];

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

// Driver code
int main()
{
    arrayOfVectors();
    return 0;
}
