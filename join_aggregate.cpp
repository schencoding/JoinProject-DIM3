#include"dim3.hpp"

myvector<TUPLE<int, int, float>> R, S, result;
/*
for table R(int x, int y, float v) and S(int z, int y, float v),
we want to calculate SQL:
    SELECT R.x,S.z,avg(R.v+S.v)
    FROM R,S
    WHERE R.y=S.y
*/

int main(int argc, char* argv[]) {
    gen_rand_data(R,S);
    result.clear();

    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    DIM3JoinAggregate<int, int, int, float, AggregateFunctionAvg<float>,
        myHasher, myHasher, myHasher, myHasher, uint32_t>
        ::doJoinAggregate(R, S, result);

    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    
    printf("It takes %f seconds.\n", chrono::duration<double, std::ratio<1, 1>>(endTime - startTime).count());
    printf("|Result|=%u\n", result.size());

    return 0;
}
