#include <vector>
#include <math.h>

using namespace std;

// Calculates the binomial of two numbers, n and k
float calculateBinomial(float n, float k);

// Calculates the number of North-East lattice paths to point
int calculateNELatticeNumber(vector<int> point);

// Calculates the number of Dyck paths using the Catalan numbers
float calculateCatalanNumber(float n);

// Calculates the number of Motzkin paths using the Motzkin numbers
float calculateMotzkinNumber(float n);

// Calculates the number of combinations without replacement
int calculateCombinationNumberWithoutReplacement(float n, float r);

// Calculates the number of combinations with replacement
int calculateCombinationNumberWithReplacement(float n, float r);

// Calculates the actual permutations of integers
vector<vector<int>> calculatePermutationsInt(int numberOfChoosenItems, vector<int> itemsChoosenFrom);

// Calculates the actual permutations of characters
vector<vector<char>> calculatePermutationsChar(int numberOfChoosenItems, vector<char> itemsChoosenFrom);

class LatticePath {
public:

    // The length of the lattice path
    int length;

    // A vector of vectors which give the possible steps
    vector<vector<int>> steps;

    // Constructor
    LatticePath(int length, vector<vector<int>> steps);
    
    // Destructor
    ~LatticePath();

    // Updates the length of the lattice path
    void setLength(int lengthOfPath);

    // Updates the possible steps of the lattice paths
    void setSteps(vector<vector<int>> possibleSteps);

    // Calculates all the possible Lattice paths that end at the given point. NE is asking if the only moves are north and east.
    vector<vector<vector<int>>> calculateLatticePaths(vector<int> point, bool dyckPath, bool motzkinPath);

    // Calculates all the possible Dyck paths
    vector<vector<vector<int>>> calculateDyckPaths();

    // Calculates all the possible Motzkin paths
    vector<vector<vector<int>>> calculateMotzkinPaths();

    // Finds all the sequences shorter than maxLength that contain the input path
    vector<vector<vector<int>>> calculateContainers(vector<vector<int>> pathToContain, vector<vector<vector<int>>> allPaths);

    // Finds all the sequences shorter than maxLength that avoid the input path
    vector<vector<vector<int>>> calculateAvoiders(vector<vector<int>> pathToAvoid, vector<vector<vector<int>>> allPaths);
};


