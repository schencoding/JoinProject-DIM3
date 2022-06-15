#include"dim3.hpp"

// Input table R(x,y), S(z,y)
myvector<pair<int, int>> R, S;

// R join S on <y>, then project and deduplicate on <x,z>
// Output table result(x,z)
myvector<pair<int, int>> result;

void gen_rand_data(
        myvector<pair<int, int>>& R, myvector<pair<int, int>>& S,
        int nR = 1e7, int nS = 1e7,
        int MODX = 1e6, int MODY = 1e6, int MODZ = 1e6,
        int seed = 492) {
    
    R.clear();
    S.clear();

    default_random_engine random_engine(seed);
    for (int i = 0; i < nR; i++) R.push_back({ random_engine() % MODX,random_engine() % MODY });
    for (int i = 0; i < nS; i++) S.push_back({ random_engine() % MODZ,random_engine() % MODY });
}

int main(int argc, char* argv[]) {
    gen_rand_data(R,S);
    result.clear();

    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    DIM3<int, int, int, myHasher, myHasher, myHasher, myHasher, uint32_t>::doJoinProject(R, S, result);

    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    
    printf("It takes %f seconds.\n", chrono::duration<double, std::ratio<1, 1>>(endTime - startTime).count());
    printf("|Result|=%u\n", result.size());

    return 0;
}
