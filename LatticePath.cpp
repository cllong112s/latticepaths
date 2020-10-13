#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <unordered_map>
#include "LatticePath.h"

using namespace std;

// This is a function to calculate the binomial coefficient of two numbers, n and k.
float calculateBinomial(float n, float k) {
    // This will be returned in the end as the binomial coefficient.
    float binomial = 1;

    // This calculates which input was larger, n or k, for later purposes.
    int maxInput = max(n, k);

    // This calculates the binomial one piece at a time.
    // I split it up like this because doing factorials not dividing it out
    // created too large to numbers.
    for (int i = 2; i <= maxInput; i++) {
        if (i <= n) {
            binomial *= i;
        }
        if (i <= k) {
            binomial /= i;
        }
        if (i <= n - k) {
            binomial /= i;
        }
    }
    // Printing the calculated BC to the console.
    printf("%i binomial %i = %i", n, k, binomial);

    // Returning the calculated BC.
    return binomial;
}

// Constructor for a lattice path object.
LatticePath::LatticePath(int length, vector<vector<int>> steps) {
    // Setting the length of the path and the steps that can be taken in the integer plane.
    this -> length = length;
    this -> steps = steps;

    // Printing the input information to the console.
    cout << "LatticePath object created with length: " << length << endl;
    cout << "LatticePath object created with steps: ";
    for (int i = 0; i < steps.size(); i++) {
        cout << "{";
        cout << steps[i][0];
        cout << ",";
        cout << steps[i][1] ;
        cout << "} ";
    }
    cout << endl;
}

// Deconstructor
LatticePath::~LatticePath() {
    cout << "Deconstructor" << endl;
}

// Sets/changes the length of the lattice path.
void LatticePath::setLength(int lengthOfPath) {
    this -> length = lengthOfPath;
    printf("Length set to %i", length);
}

// Sets/changes the steps allowed for the path.
void LatticePath::setSteps(vector<vector<int>> possibleSteps) {
    this -> steps = possibleSteps;
}

// This calculates the number of North-East lattice paths to a point.
int calculateNELatticeNumber(vector<int> point) {
    int numberOfPaths;

    // Setting the X and Y values that the path must reach.
    float x = point[0];
    float y = point[1];

    // This value is calculated with a BC from above.
    numberOfPaths = calculateBinomial(x + y, x);

    // Printing out the calculated number of NE lattice paths to the console.
    printf("Number of NE lattice paths = %i \n", numberOfPaths);

    // Returning the calculated number of paths.
    return numberOfPaths;
}

// Calculates the Catalan number of a given length. This is also the number of Dyck paths of n = length / 2.
float calculateCatalanNumber(float n) {
    // The Catalan number is calculated with two pieces, one being a BC.
    float binomial = calculateBinomial(2 * n, n);

    // Printing the calculated Catalan number to the console.
    printf("Catalan number for %.1f = %.1f \n", n, binomial/(n + 1));

    // Returning the calculated Catalan number.
    return binomial / (n + 1);
}

// Calculates the Motzkin number/number of Motzkin paths of a given length.
float calculateMotzkinNumber(float n) {
    int numberOfPaths;

    // The upper limit of the sum to find the Motzkin involves the floor function.
    int floorND2 = floor(n / 2);

    // Calculating the floor of n / 2.
    for (float k = 0; k <= floorND2; k++) {
        numberOfPaths += (calculateBinomial(n, (2 * k)) * calculateCatalanNumber(k));
    }

    // Printing the Motzkin number/number of paths to the console.
    printf("Number of Motzkin paths of length %i = %i \n", n, numberOfPaths);

    // Returning the Motzkin number/number of paths.
    return numberOfPaths;
}

// This calculates the number of possible combinations without replacement given n and r.
int calculateCombinationNumberWithoutReplacement(float n, float r) {

    // The number of combinations without replacement is given by a simple BC.
    float combinationsWithoutReplacement = calculateBinomial(n ,r);

    // Printing the number of combinations without replacement to the console.
    printf("Combinations without replacement using n = %.1f, r = %.1f : %.1f \n", n, r, combinationsWithoutReplacement);

    // Returning the number of combinations without replacement.
    return combinationsWithoutReplacement;
}

// This calculates the number of combinations with replacement.
int calculateCombinationNumberWithReplacement(float n, float r) {
    // This is set to 1 because we are going to use multiplication, where 1 is the identity.
    float combinationsWithReplacement = 1;

    // Getting the max input to find the stopping iteration for the for loop.
    float maxInput = max(n + r - 1, r);
    maxInput = max(maxInput, n - 1);

    // Calculating this involves factorials, which can get very large, so I split the calculation
    // into pieces so the number would stay small.
    for (int i = 2; i <= maxInput; i++) {
        if (i <= n + r - 1) {
            combinationsWithReplacement *= i;
        }
        if (i <= r) {
            combinationsWithReplacement /= i;
        }
        if (i <= n - 1) {
            combinationsWithReplacement /= i;
        }
    }

    // Printing the number of combinations with replacement to the console.
    printf("Combinations with replacement using n = %.1f, r = %.1f : %.1f \n", n, r, combinationsWithReplacement);

    // Returning the combinations with replacement.
    return combinationsWithReplacement;
}

// I believe that I need to set these out here because I am going to be using helper functions.
// I may be wrong, please let me know if I can improve this.
int lengthOfPermutation; // Sets the length of the permutation in question.

// Permutations list for integers.
vector<vector<int>> permutationsInt;

// Permutations list for characters.
vector<vector<char>> permutationsChar;

// A vector that will hold numbers which act as variables in the helper functions.
// For example, if index 3 is 7, then the variable that is index 3 is set to 7.
// This is so I can use many variables without creating them manually.
vector<int> vars;

// This function will be called recursively to find all integer permutations.
permutationsHelperFunctionInt(int counter, int numberOfChoosenItems, vector<int> itemsChoosenFrom) {
    // For each vars variable position until the number of items in the item set - 1, push values from the itemsChoosenFrom vector
    // into the current permutations vector and then the permutations list. After each iteration the vars variable position is incremented.
    // This is part of a multi-level for where this function is called many times.
    for (vars[counter] = 0; vars[counter] < itemsChoosenFrom.size(); vars[counter]++) {
        // If the next item is equal to the number of items being picked for the permutation.
        if (counter + 1 == numberOfChoosenItems) {
            // Setting and clearing the current permutation.
            vector<int> currentPermutation;
            currentPermutation.clear();
            // Pushing the items into the current permutation vector.
            for (int i = 0; i < numberOfChoosenItems; i++) {
                currentPermutation.push_back(itemsChoosenFrom[vars[i]]);
            }
            // Pushing the currentPermutation onto the permutations list for the integers.
            permutationsInt.push_back(currentPermutation);
        } else {
            // This is the recursive call that layers the for loop. We do counter + 1 here
            // because we are wanting to go to the next vars variable position.
            permutationsHelperFunctionInt(counter + 1, numberOfChoosenItems, itemsChoosenFrom);
        }
    }
}

// This is a function like the one above but for characters. I don't believe I could use template here
// because of the other data types in the function. Let me know if this is incorrect please.
permutationsHelperFunctionChar(int counter, int numberOfChoosenItems, vector<char> itemsChoosenFrom) {
    // For each vars variable position until the number of items in the item set - 1, push values from the itemsChoosenFrom vector
    // into the current permutations vector and then the permutations list. After each iteration the vars variable position is incremented.
    // This is part of a multi-level for loop where this function is called many times.
    cout << vars.size() << endl;
    for (vars[counter] = 0; vars[counter] < itemsChoosenFrom.size(); vars[counter]++) {
        cout << vars[counter] << endl;
        // If the next item is equal to the number of items being picked for the permutation.
        if (counter + 1 == numberOfChoosenItems) {
            // Setting and clearing the current permutation.
            vector<char> currentPermutation;
            currentPermutation.clear();
            // Pushing the items into the current permutation vector.
            for (int i = 0; i < numberOfChoosenItems; i++) {
                currentPermutation.push_back(itemsChoosenFrom[vars[i]]);
            }
            // Pushing the currentPermutation onto the permutations list for the characters.
            permutationsChar.push_back(currentPermutation);
        } else {
            // This is the recursive call that layers the for loop. We do counter + 1 here
            // because we are wanting to go to the next vars variable position.
            permutationsHelperFunctionChar(counter + 1, numberOfChoosenItems, itemsChoosenFrom);
        }
    }
}

// This will calculate all the permutations of consecutive integers. This will produce permutations of length numberOfChoosenItems,
// from itemsChoosenFrom.
vector<vector<int>> calculatePermutationsInt(int numberOfChoosenItems, vector<int> itemsChoosenFrom) {
    vector<vector<int>> outputVector;
    // Defining vars to keep my variables used in the recursive for loop calls.
    vector<int> vars;
    int counter = 0;
    int numberOfPermutations;

    // Filling the vars vector with 0's.
    for (int i = 0; i < numberOfChoosenItems; i++) {
        vars.push_back(0);
    }

    // First call to the helper function, which is all will produce the permutations.
    permutationsHelperFunctionInt(counter, numberOfChoosenItems, itemsChoosenFrom);
    outputVector = permutationsInt;
    permutationsInt.clear();
    vars.clear();

    // Printing out the number of permutations to the console.
    cout << outputVector.size() << endl;

    printf("Number of integer permutations with %i items choosen from %i = %i \n", numberOfChoosenItems, itemsChoosenFrom.size(), outputVector.size());

    // Returning the vector of permutations.
    return outputVector;
}

// This calculates the permutations of characters in the same way as the integer permutations does.
vector<vector<char>> calculatePermutationsChar(int numberOfChoosenItems, vector<char> itemsChoosenFrom) {
    vector<vector<char>> outputVector;

    // Defining the vars vector. This will keep numberOfItemsChoosen variable positions.
    int counter = 0;
    int numberOfPermutations;

    // Setting each place in the vars vector to 0, because each variable starts at zero.
    for (int k = 0; k < numberOfChoosenItems; k++) {
        vars.push_back(0);
    }

    // First call to the helper function which will produce all permutations of the characters.
    permutationsHelperFunctionChar(counter, numberOfChoosenItems, itemsChoosenFrom);
    outputVector = permutationsChar;
    permutationsChar.clear();
    vars.clear();

    // Printing out the number of permutations to the console.
    printf("Number of character permutations with %i items choosen from %i = %i \n", numberOfChoosenItems, itemsChoosenFrom.size(), outputVector.size());

    // Returning the vector of all permutations.
    return outputVector;
}

// Defining the vector that will hold all lattice paths. I am doing this here because I am using a helper function.
// Please let me know if this is incorrect.
vector<vector<vector<int>>> LatticePaths;

// This is the helper function to be called to find all lattice paths with certain steps, to a certain point in the integer plane.
// dyckPath and motzkinPath tell the function if the lattice path is of a special type.
LatticePathsHelperFunction(int counter, vector<vector<int>> steps, vector<int> point, bool dyckPath, bool motzkinPath) {
    // We are using a vars vector here again because we are using a recursively called for loop to generate permutations.
    for (vars[counter] = 0; vars[counter] < steps.size(); vars[counter]++) {
        if (counter + 1 == lengthOfPermutation) { // If this iteration gives us the length of the path we need...
            vector<vector<int>> currentPermutation; // Defining the current permutation vector.
            currentPermutation.clear(); // Clearing the current permutation vector in case it has values already in it from a call beforehand.
            // Filling the current Permutation vector with the steps as defined by the vars vector.
            for (int i = 0; i < lengthOfPermutation; i++) {
                currentPermutation.push_back(steps[vars[i]]);
            }

            // We need to check if the path made it to the point or not. So we define the x and y position.
            int sumX = 0;
            int sumY = 0;
            // Now we look to see if the paths got to the point.
            for (int i = 0; i < currentPermutation.size(); i++) {
                sumX += currentPermutation[i][0];
                sumY += currentPermutation[i][1];

                // A Dyck path is a special lattice path that cannot go above the x = y line, so here we check for that.
                if (dyckPath && sumY < 0) {
                    break;
                }

                // A motzkin path cannot fall below the x-axis, here we check for that.
                if (motzkinPath && sumY < 0) {
                    break;
                }
            }
            // If the path got the point, place the current permutation into the lattice paths vector to be outputted.
            if (sumX == point[0] && sumY == point[1]) {
                LatticePaths.push_back(currentPermutation);

            }

        // If the length isn't there yet, call the function again.
        } else {
            LatticePathsHelperFunction(counter + 1, steps, point, dyckPath, motzkinPath);
        }
    }
}

// This function will find all the pattice paths using the steps provided going to the point provided.
// You can also state if the path is a special Dyck path or Motzkin path.
vector<vector<vector<int>>> LatticePath::calculateLatticePaths(vector<int> point, bool dyckPath, bool motzkinPath) {
    if (dyckPath) {
        setSteps({{1,1},{1,-1}});
    }
    if (motzkinPath) {
        setSteps({{1,1},{1,-1},{1,0}});
    }

    vector<vector<vector<int>>> outputVector;
    int counter = 0;
    lengthOfPermutation = length;

    // This sets each index in vars to a value of 0. This will be used to count for loop iterations later.
    for (int i = 0; i < length; i++) {
        vars.push_back(0);
    }

    // Calling the recursive helper function to find all lattice paths.
    LatticePathsHelperFunction(counter, steps, point, dyckPath, motzkinPath);
    outputVector = LatticePaths;

    LatticePaths.clear();
    vars.clear();

    // This prints out the number of lattice paths found.
    printf("Number of lattice paths with the given steps and length %i = %i \n", length, outputVector.size());

    // Returning the vector of all lattice paths fitting the steps, and point.
    return outputVector;
}

// This is to calculate a special lattice path called a dyck path. This just calls calculateLatticePaths with the dyckPath set to true.
// Dyck paths must be of even length.
vector<vector<vector<int>>> LatticePath::calculateDyckPaths() {
    vector<vector<vector<int>>> dyckPaths = calculateLatticePaths({length, 0}, true, false);

    // Prints the number of Dyck paths to the terminal
    printf("Number of Dyck paths with length %i = %i \n", length, dyckPaths.size());

    // Returning the calculated Dyck paths
    return dyckPaths;
}

// This calculates all the Motzkin paths of a certain length.
vector<vector<vector<int>>> LatticePath::calculateMotzkinPaths() {
    vector<vector<vector<int>>> motzkinPaths = calculateLatticePaths({length, 0}, false, true);

    // Printing out the number of Motzkin paths fitting the inputs.
    printf("Number of Motzkin paths with length %i = %i \n", length, motzkinPaths.size());

    // Returning the vector of Motzkin paths found.
    return motzkinPaths;
}

// A path A is avoided by a path B if the steps in path A are not found in order within path B.
vector<vector<vector<int>>> LatticePath::calculateAvoiders(vector<vector<int>> pathToAvoid, vector<vector<vector<int>>> allPaths) {
    vector<vector<vector<int>>> avoiders;
    // This looks through each path within allPaths and checks for avoidance to the pathToAvoid.
    for (int i = 0; i < allPaths.size(); i++) {
        int j = 0;
        int k = 0;
        int matchCounter = 0;
        while (j < allPaths[i].size()) {
            if (allPaths[i][j] == pathToAvoid[k]) {
                matchCounter++;

                k++;
            }
            j++;
        }
        // Checking to see if avoidance is met.
        if (matchCounter < pathToAvoid.size()) {
            avoiders.push_back(allPaths[i]);
        }
    }
    // Printing the number of avoiding paths to the console.
    printf("Paths avoiding the input path: %i \n", avoiders.size());

    // Returning the vector of avoiders.
    return avoiders;
}

// A path A is said to contain path B if the steps in path A are all found in order within path B.
vector<vector<vector<int>>> LatticePath::calculateContainers(vector<vector<int>> pathToContain, vector<vector<vector<int>>> allPaths) {
    vector<vector<vector<int>>> containers;
    // For each path in allPaths, check if that path is contained within the given path. This just means that the
    // pathToContain's steps are all given one of the paths in allPaths in order. For example, take the path ABCD, AC is contained with the path.
    for (int i = 0; i < allPaths.size(); i++) {
        int j = 0;
        int k = 0;
        int matchCounter = 0;
        // This is looking how many matching steps there are (in order) that are within one of the allPaths. If all the steps from pathToContain
        // are within the allPath vector then the path is said to contain the path.
        while (j < allPaths[i].size()) {
            if (allPaths[i][j] == pathToContain[k]) {
                matchCounter++;
                k++;
            }
            j++;
        }
        // Checking to see if containment was met.
        if (matchCounter >= pathToContain.size() - 1) {
            containers.push_back(allPaths[i]);
        }
    }
    printf("Paths containing the input path: %i \n", containers.size());

    // Returning the paths that contained the path inputted.
    return containers;
}

