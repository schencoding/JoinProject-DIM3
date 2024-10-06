// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// DIM3 is licensed under Mulan PSL v2.

#pragma once

#define THREAD_NUM 8

#include <omp.h>
#include "dim3.hpp"

template <typename _Tx, /*typename _Ty,*/ typename _Tz, typename _Tcounter = uint32_t>
struct Hybrid_solution_parallel {
public:
    static void dojoinproject(CSR<>& R_all, S_DENSE& S_dense, CSR<>& S_sparse, ID2VALUE<_Tx,_Tz>& ID2V,
                                const int nx, const int ny, const int nz,
                                myvector<pair<_Tx,_Tz>, _Tcounter>& result, size_t reserve_size/*?*/) {
        result.reserve(reserve_size);

        myvector<pair<_Tx, _Tz>, _Tcounter> local_result[THREAD_NUM];
        for (int i = 0; i < THREAD_NUM; i++) {
            local_result[i].init(result.max_num / THREAD_NUM * 2);
        }

        if (S_dense.size > 0) {
            //DenseEC
            int* origin_z_val = (int*)malloc(S_dense.size * sizeof(int));
            if (origin_z_val == NULL) {
                printf("Fail to malloc origin_z_val in SPMX.\n");
                exit(-1);
            }
            for (int i = 0; i < S_dense.size; i++) {
                origin_z_val[i] = ID2V.z[S_dense.id2v[i]];
            }

#pragma omp parallel num_threads(THREAD_NUM)
            {
                auto tid = omp_get_thread_num();
                auto my_tmp_result = local_result[tid];
                int dense_bitvector_len = S_dense.dense_bitvector_len;
                myBitVector<int> dense_x_bitvector(dense_bitvector_len);
                int* f3_threshold_cache= (int*)malloc(ny * sizeof(int));
                memset(f3_threshold_cache,-1,ny * sizeof(int));
#pragma omp for schedule(dynamic,30)
                for (int i = 0; i < nx; i++) {
                    if (R_all.JR[i + 1] - R_all.JR[i] == 0) continue;
                    int origin_x_val = ID2V.x[i];
                    bool is_dense_bitvector_initialized = false;
                    int _t=R_all.JR[i + 1] - R_all.JR[i];
                    if (f3_threshold_cache[_t]==-1) f3_threshold_cache[_t]=f3_threshold::cal_SIMD_method_threshold(_t, ny);
                    int nsthreshold = f3_threshold_cache[_t];
                    for (int j = 0; j < S_dense.size; j++) {
                        if (S_dense.id2cnt[j] > nsthreshold) {
                            //SIMD method
                            if (!is_dense_bitvector_initialized) {
                                memset(dense_x_bitvector._Array, 0, dense_bitvector_len * sizeof(int));
                                for (int j = R_all.JR[i]; j < R_all.JR[i + 1]; j++) {
                                    dense_x_bitvector.setTrue(R_all.IC[j]);
                                }
                                is_dense_bitvector_initialized = true;
                            }

                            int* z_bitvector_array = S_dense.data[j]._Array;
                            for (int t = 0; t < dense_bitvector_len; t += 8) {
                                __m256i a = _mm256_loadu_si256((__m256i*)(z_bitvector_array + t));
                                __m256i b = _mm256_loadu_si256((__m256i*)(dense_x_bitvector._Array + t));
                                if (_mm256_testz_si256(a, b) == 0) {
                                    my_tmp_result.push_back({ origin_x_val,origin_z_val[j] });
                                    break;
                                }
                            }
                        }
                        else {
                            //non-SIMD method
                            for (int t = R_all.JR[i]; t < R_all.JR[i + 1]; t++) {
                                if (S_dense.data[j].checkIfTrue(R_all.IC[t])) {
                                    my_tmp_result.push_back({ origin_x_val,origin_z_val[j] });
                                    break;
                                }
                            }
                        }
                    }
                }
                free(f3_threshold_cache);
                free(dense_x_bitvector._Array);
                local_result[tid] = my_tmp_result;
            }
            free(origin_z_val);
        }

        if (S_sparse.JR != NULL) {
#pragma omp parallel num_threads(THREAD_NUM)
            {
                auto tid = omp_get_thread_num();
                auto my_tmp_result = local_result[tid];
                int *SPAw = (int*)malloc(nz * sizeof(int));
                if (SPAw == NULL) {
                    printf("Fail to malloc SPAw in SPMX.\n");
                    exit(-1);
                }
                memset(SPAw, -1, nz * sizeof(int));

#pragma omp for schedule(dynamic,30)
                for (int cur_x = 0; cur_x < nx; cur_x++) {
                    int origin_x_val = ID2V.x[cur_x];
                    for (int i = R_all.JR[cur_x], _i = R_all.JR[cur_x + 1]; i < _i; i++) {
                        int cur_y = R_all.IC[i];
                        for (int j = S_sparse.JR[cur_y], _j = S_sparse.JR[cur_y + 1]; j < _j; j++) {
                            if (SPAw[S_sparse.IC[j]] != cur_x) {
                                SPAw[S_sparse.IC[j]] = cur_x;
                                my_tmp_result.push_back({ origin_x_val,ID2V.z[S_sparse.IC[j]] });
                            }
                        }
                    }
                }
                free(SPAw);
                local_result[tid] = my_tmp_result;
            }
        }

        {
            int result_cnt = 0;
            int thread_result_base[THREAD_NUM];
            for (int i = 0; i < THREAD_NUM; i++) {
                thread_result_base[i] = result_cnt;
                result_cnt += local_result[i].cur_num;
            }
            result.reserve(result_cnt);
#pragma omp parallel for num_threads(THREAD_NUM)
            for (int i = 0; i < THREAD_NUM; i++) {
                memcpy(result.data + thread_result_base[i], local_result[i].data, local_result[i].cur_num * sizeof(pair<_Tx, _Tz>));
            }
            result.cur_num = result_cnt;
        }
    }
};

template <typename _Tx, typename _Ty, typename _Tz,
    typename _Hx = hash<_Tx>, typename _Hy = hash<_Ty>, typename _Hz = hash<_Tz>,
    typename _Hxz = hash<pair<_Tx, _Tz>>, typename _Tcounter = uint32_t>
struct DIM3_parallel : public DIM3<_Tx,_Ty,_Tz,_Hx,_Hy,_Hz,_Hxz,_Tcounter> {
public:
    static void doJoinProject_parallel(const myvector<pair<_Tx, _Ty>, _Tcounter> _R,
        const myvector<pair<_Tz, _Ty>, _Tcounter> _S,
        myvector<pair<_Tx, _Tz>, _Tcounter>& result) {

        result.clear();
        myvector<pair<_Tx, _Ty>, _Tcounter> R,S;
        if (_R.size()>=_S.size()) {R=_R;S=_S;}
        else {R=_S;S=_R;}
        // R larger

        LL OUT_J_hat = EstimateJoinCardinality<_Tx, _Ty, _Tz, _Hy, _Tcounter>::estimate(R, S);
        // cout << "[i] OUT_J_hat= " << OUT_J_hat << endl;

        if (f1_threshold::is_classical_batter(R.size(), S.size(), OUT_J_hat)) {
            //Radix hash
            // cout << "[i] Use Radix hash" << endl;

            printf("[E] Parallel Radix_hash_join_and_project not implemented\n");
            assert(false);
            return;
        }

        //hybrid solution
        // cout << "[i] Use hybrid solution" << endl;

        int n, k, m;
        int xmax, ymax, zmax;
        S_DENSE S_dense;
        CSR<> S_sparse;
        CSR<> R_all;

        ID2VALUE<_Tx, _Tz> ID2V;
        myvector<pair<int, int>> Rm, Sm;

        //in RDBMS when we know the |X||Y||Z|, calll mapping
        mapping_step<_Tx, _Ty, _Tz, _Hx, _Hy, _Hz>::mapping(R, S, ID2V, Rm, Sm, xmax, ymax, zmax);
        
        // cout << "[i] |Rm|=" << Rm.size() << " |Sm|=" << Sm.size() << endl;

        if (Rm.size() == 0 || Sm.size() == 0)  return;

        n = xmax + 1;
        k = ymax + 1;
        m = zmax + 1;

        // cout << "[i] n=" << n << " k=" << k << " m=" << m << endl;

        DIM3<_Tx,_Ty,_Tz,_Hx,_Hy,_Hz,_Hxz,_Tcounter>::divide_R(Rm, n, R_all);
        DIM3<_Tx,_Ty,_Tz,_Hx,_Hy,_Hz,_Hxz,_Tcounter>::divide_S(Sm, n, k, m, R_all, S_dense, S_sparse);

        free(Rm.data);
        free(Sm.data);

        Hybrid_solution_parallel<_Tx, _Tz, _Tcounter>::dojoinproject(R_all, S_dense, S_sparse, ID2V, n, k, m, result, min(OUT_J_hat*3, (LL)n * m));

        // cout << "[i] |S_dense|=" << S_dense.size << " (" << (double)S_dense.size / m * 100. << "%)" << endl;

        for (int i = 0, _i = S_dense.size; i < _i; i++) {
            free(S_dense.data[i]._Array);
        }
        free(S_dense.id2cnt);
        free(S_dense.id2v);
        free(S_sparse.JR);
        free(S_sparse.IC);
        free(R_all.JR);
        free(R_all.IC);

    }
};
