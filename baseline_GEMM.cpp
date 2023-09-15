#include"dim3.hpp"
#include<mkl.h>

// Input table R(x,y), S(z,y)
myvector<pair<int, int>> R, S;

// R join S on <y>, then project and deduplicate on <x,z>
// Output table result(x,z)
myvector<pair<int, int>> result;

int main(int argc, char* argv[]) {
    // Note: We can't use the default parameters here 
    // because that's a relatively sparse dataset and 
    // would result in oversized matrices!
    gen_rand_data(R,S,1e7,1e7,1e4,1e4,1e4); 
    printf("Note that a different dataset is used here.\n");

    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    int n, k, m;
    int xmax, ymax, zmax;

    ID2VALUE<int, int> ID2V;
    myvector<pair<int, int>> Rm, Sm;

    long long OUT_J_hat = EstimateJoinCardinality<int, int, int, myHasher, uint32_t>::estimate(R, S);

    baseline_mapping(R, S, ID2V, Rm, Sm, xmax, ymax, zmax);

    n = xmax + 1;
    k = ymax + 1;
    m = zmax + 1;

    char* A = (char*)mkl_calloc(n * k, sizeof(char), 64);
    unsigned char* B = (unsigned char*)mkl_calloc(k * m, sizeof(unsigned char), 64);
    int *C = (int*)mkl_calloc(n * m, sizeof(int), 64);

    if (A == NULL || B == NULL || C == NULL) {
        printf("\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
        mkl_free(A);
        mkl_free(B);
        mkl_free(C);
        exit(-1);
    }

    for (int i = 0, _i = Rm.cur_num; i < _i; i++) {
        A[Rm.data[i].first * k + Rm.data[i].second] = 1;
    }

    for (int i = 0, _i = Sm.cur_num; i < _i; i++) {
        B[Sm.data[i].second * m + Sm.data[i].first] = 1;
    }

    MKL_INT32 oc = 0;
    cblas_gemm_s8u8s32(CblasRowMajor, CblasNoTrans, CblasNoTrans, CblasFixOffset,
        n, m, k, 1.0, A, k, 0, B, m, 0, 0.0, C, m, &oc);

    result.clear();
    result.reserve(min(OUT_J_hat,(long long)n*m));
    int ii = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (C[ii] > 0) {
                result.push_back({ ID2V.x[i],ID2V.z[j] });
            }
            ii++;
        }
    }

    chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();
    
    printf("It takes %f seconds.\n", chrono::duration<double, std::ratio<1, 1>>(endTime - startTime).count());
    printf("|Result|=%u\n", result.size());

    return 0;
}
