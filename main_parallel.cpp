#include"dim3_parallel.hpp"

// Input table R(x,y), S(z,y)
myvector<pair<int, int>> R, S;

// R join S on <y>, then project and deduplicate on <x,z>
// Output table result(x,z)
myvector<pair<int, int>> result;

int main(int argc, char* argv[]) {
    gen_rand_data(R,S);
    result.clear();

    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    DIM3_parallel<int, int, int, myHasher, myHasher, myHasher, myHasher, uint32_t>::doJoinProject_parallel(R, S, result);

    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    
    printf("It takes %f seconds.\n", chrono::duration<double, std::ratio<1, 1>>(endTime - startTime).count());
    printf("|Result|=%u\n", result.size());

    return 0;
}
