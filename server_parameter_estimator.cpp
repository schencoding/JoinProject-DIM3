// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// DIM3 is licensed under Mulan PSL v2.

#include <bits/stdc++.h>
#include <immintrin.h>
#include <chrono>
#include "Include/flat_hash_map.hpp"
#include "Include/myvector.hpp"
#include "Include/structures.hpp"
#include "Include/L+S_mapper.hpp"
#include "Include/myBitVector.hpp"

using namespace std;

//===============
//L1: 32KB?
//L2: 256KB?
//L3: 12MB?
const unsigned int L1cache_size = 32768; // in bytes
const unsigned int L2cache_size = 262144; // in bytes
const unsigned int L3cache_size = 12582912; // in bytes
//===============

default_random_engine random_engine(492);

void part4_func(int* p, int* q, unsigned int n, string s) {
    const int N = 1e2;
    const int nR = 1e7;
    for (unsigned int j = 0; j < nR; j++) p[j] = random_engine() % n;
    for (unsigned int j = 0; j < n; j++) q[j] = random_engine();

    chrono::high_resolution_clock::time_point t1, t2;
    chrono::duration<double, std::ratio<1, 1>> d1, d2;

    int a = 0;
    t1 = chrono::high_resolution_clock::now();
    for (int nn = 0; nn < N; nn++) {
        for (int j = 0; j < nR; j++) {
            a += q[p[j]];
        }
    }
    t2 = chrono::high_resolution_clock::now();
    d1 = t2 - t1;

    int b = random_engine.max() / 2;
    t1 = chrono::high_resolution_clock::now();
    for (int nn = 0; nn < N; nn++) {
        for (int j = 0; j < nR; j++) {
            q[p[j]] += nn;
        }
    }
    t2 = chrono::high_resolution_clock::now();
    d2 = t2 - t1;

    fprintf(stderr, "%d %u %d %d %f %f\n", a, n, q[105], q[1234], d1.count(), d2.count());
    printf("const double Para_Sz_mem_rand_read_%s = %f; // 1e-9 seconds\n", s.c_str(), d1.count());
    printf("const double Para_Sz_mem_rand_read_write_%s = %f; // 1e-9 seconds\n", s.c_str(), d2.count());
    fflush(stdout);
}

int main(int argc, char* argv[]) {
    {
        // Part1
        printf("const unsigned int __L1cache_size = %u;\n", L1cache_size);
        printf("const unsigned int __L2cache_size = %u;\n", L2cache_size);
        printf("const unsigned int __L3cache_size = %u;\n", L3cache_size);
        printf("const unsigned int __L1cache_size_by_int = __L1cache_size / sizeof(int);\n");
        printf("const unsigned int __L2cache_size_by_int = __L2cache_size / sizeof(int);\n");
        printf("const unsigned int __L3cache_size_by_int = __L3cache_size / sizeof(int);\n");
        printf("const unsigned int __cache_line_size = 64;\n");
        printf("const unsigned int __partition_para = 10000;\n");
        printf("\n");
        printf("\n");
        fflush(stdout);
    }

    {
        // Part2
        myvector<pair<int, int>> R, S;
        ska::flat_hash_set<pair<int, int>, myHasher> fset;
        int N = 1e7;
        int M = 1e2;
        int dataRangeN = 1e8;
        int dataRangeK = 1e8;
        int dataRangeM = 1e8;

        chrono::high_resolution_clock::time_point t1, t2;
        chrono::duration<double, std::ratio<1, 1>> d1, d2, d3, d4;

        R.clear();S.clear();
        default_random_engine random_engine(492);
        for (int i = 0; i < N; i++) R.push_back({ random_engine() % dataRangeN,random_engine() % dataRangeK });
        for (int i = 0; i < N; i++) S.push_back({ random_engine() % dataRangeM,random_engine() % dataRangeK });

        int* p = (int*)malloc(dataRangeK * sizeof(int));
        memset(p, -1, dataRangeK * sizeof(int));

        t1 = chrono::high_resolution_clock::now();
        for (int nn = 0; nn < M; nn++) {
            for (int i = 0; i < N; i++) {
                p[R.data[i].second] = nn;
            }
        }
        t2 = chrono::high_resolution_clock::now();
        d1 = (t2 - t1) * (1e9 / (double)M / (double)N);

        fset.reserve(3 * N);
        fset.clear();
        t1 = chrono::high_resolution_clock::now();
        for (int i = 0; i < N; i++) {
            fset.insert({ R.data[i].first,R.data[i].second });
        }
        t2 = chrono::high_resolution_clock::now();
        d2 = (t2 - t1) * (1e9 / (double)N);

        LS_mapper<int, myHasher> hashmap(N, 0, 16384 / 4);
        t1 = chrono::high_resolution_clock::now();
        for (int i = 0; i < N; i++) {
            S.data[i].first = R.data[i].first;
            S.data[i].second = hashmap.insert(R.data[i].second, p);
            //R.data[i] = { R.data[i].first,hashmap.insert(R.data[i].second) };
        }
        t2 = chrono::high_resolution_clock::now();
        d3 = (t2 - t1) * (1e9 / (double)N);

        int a = S.data[107].second;

        t1 = chrono::high_resolution_clock::now();
        for (int i = 0, j = 0; i < N; i++) {
            j++;
            S.data[j].first = hashmap.find(R.data[i].first);
            S.data[j].second = R.data[i].second;
        }
        t2 = chrono::high_resolution_clock::now();
        d4 = (t2 - t1) * (1e9 / (double)N);

        a += S.data[107].second;
        fprintf(stderr, "%d\n", a);
        fprintf(stderr, "%f\n", d1.count());

        printf("const double Para_Radix_Flat_Hash_Table_cost = %f; // 1e-9 seconds\n", d2.count());
        printf("const double Para_Radix_LSMapper_build_cost = %f; // 1e-9 seconds\n", d3.count());
        printf("const double Para_Radix_LSMapper_probe_cost = %f; // 1e-9 seconds\n", d4.count());
        printf("\n");
        fflush(stdout);
    }

    {
        // Part3
        myvector<pair<int, int>> R, S;
        int N = 1e2;
        int nR = 1e6;
        int nS = nR;
        int dataRangeN = 6e7;
        int dataRangeK = 6e7;
        int dataRangeM = 6e7;

        chrono::high_resolution_clock::time_point t1, t2;
        chrono::duration<double, std::ratio<1, 1>> d1, d2, d3, d4;

        R.clear();S.clear();
        default_random_engine random_engine(492);
        for (int i = 0; i < nR; i++) R.push_back({ random_engine() % dataRangeN,random_engine() % dataRangeK });
        for (int i = 0; i < nS; i++) S.push_back({ random_engine() % dataRangeM,random_engine() % dataRangeK });

        N = 1e3;
        int a = 0;
        t1 = chrono::high_resolution_clock::now();
        for (int nn = 0; nn < N; nn++) {
            for (int i = 0; i < nR; i++) {
                a += 3;
            }
        }
        t2 = chrono::high_resolution_clock::now();
        d1 = t2 - t1;

        N = 1e2;
        int bitvector_len = (dataRangeK / myBitVector<>::_Bitsperword) + (dataRangeK % myBitVector<>::_Bitsperword == 0 ? 0 : 1);
        assert(bitvector_len % 8 == 0);
        myBitVector<> bvR = myBitVector<>(bitvector_len);
        myBitVector<> bvS = myBitVector<>(bitvector_len);
        for (int i = 0; i < nR; i++) {
            bvR.setTrue(R.data[i].second);
        }
        for (int i = 0; i < nS; i++) {
            bvS.setTrue(S.data[i].second);
        }
        int b = 0;
        t1 = chrono::high_resolution_clock::now();
        for (int nn = 0; nn < N; nn++) {
            for (int t = 0; t + 8 < bitvector_len; t += 8) {
                __m256i ma = _mm256_loadu_si256((__m256i*)(bvR._Array + t));
                __m256i mb = _mm256_loadu_si256((__m256i*)(bvS._Array + t));
                if (_mm256_testz_si256(ma, mb) == 0) b += 3;
                else b += 1;
            }
        }
        t2 = chrono::high_resolution_clock::now();
        d2 = (t2 - t1) * (1e9 / (double)N / (double)(bitvector_len / 8));

        N = 1e2;
        int c = 0;
        t1 = chrono::high_resolution_clock::now();
        for (int nn = 0; nn < N; nn++) {
            for (int i = 0; i < nR; i++) {
                if (bvS.checkIfTrue(R.data[i].second)) c += 3;
                else c += 1;
                if (bvR.checkIfTrue(S.data[i].second)) c += 3;
                else c += 1;
            }
        }
        t2 = chrono::high_resolution_clock::now();
        d3 = (t2 - t1) * (1e9 / (double)N / (double)nR /2.);

        t1 = chrono::high_resolution_clock::now();
        bvR.clear(bitvector_len);
        for (int i = 0; i < nR; i++) bvR.setTrue(R.data[i].first);
        bvS.clear(bitvector_len);
        for (int i = 0; i < nS; i++) bvS.setTrue(R.data[i].first);
        bvR.clear(bitvector_len);
        for (int i = 0; i < nR; i++) bvR.setTrue(R.data[i].second);
        bvS.clear(bitvector_len);
        for (int i = 0; i < nS; i++) bvS.setTrue(R.data[i].second);
        t2 = chrono::high_resolution_clock::now();
        d4 = (t2 - t1) * (1e9 / (double)nR / 4.);

        fprintf(stderr, "a=%d\n", a);
        fprintf(stderr, "b=%d\n", b);
        fprintf(stderr, "c=%d\n", c);

        fprintf(stderr, "%f\n", d1.count());

        printf("const double Para_Rx_Dense = %f; // 1e-9 seconds\n", d2.count());
        printf("const double Para_Rx_Sparse = %f; // 1e-9 seconds\n", d3.count());
        printf("const double Para_Sz_Bitmap = %f; // 1e-9 seconds\n", d4.count());
        printf("\n");
        fflush(stdout);
    }

    {
        // Part4
        int* p = (int*)malloc(5e7 * sizeof(int));
        assert(p != NULL);
        int* q = (int*)malloc(5e7 * sizeof(int));
        assert(q != NULL);

        for (unsigned int j = 0; j < 50000000; j++) q[j] = random_engine();
        chrono::high_resolution_clock::time_point t1, t2;
        chrono::duration<double, std::ratio<1, 1>> d1;
        int a = 0;
        t1 = chrono::high_resolution_clock::now();
        for (int nn = 0; nn < 20; nn++) {
            for (int j = 0; j < 50000000; j++) {
                a += q[j];
            }
        }
        t2 = chrono::high_resolution_clock::now();
        d1 = t2 - t1;
        fprintf(stderr, "%d %f\n", a, d1.count());
        printf("const double Para_Sz_mem_seq_read = %f; // 1e-9 seconds\n", d1.count());

        part4_func(p, q, 1000, "1k");
        part4_func(p, q, L2cache_size / 4 * 2, "2L2");
        part4_func(p, q, L3cache_size / 4 / 2, "05L3");
        part4_func(p, q, L3cache_size / 4 * 1.5, "15L3");
        part4_func(p, q, L3cache_size / 4 * 3, "3L3");
        part4_func(p, q, L3cache_size / 4 * 10, "10L3");
        printf("\n");
    }

    return 0;
}
