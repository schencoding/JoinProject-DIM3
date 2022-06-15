#pragma once

#include <bits/stdc++.h>
#include <unordered_map>
#include <immintrin.h>
#include <chrono>
#include "Include/DataLoader.hpp" //To be rm
#include "Include/flat_hash_map.hpp"
#include "Include/structures.hpp"
#include "Include/L+S_mapper.hpp"
#include "Include/myBitVector.hpp"

#include "machine_parameters.hpp"

using namespace std;

//===============================================
/*
* Threshold functions 
* f1,f2,f3
*/

struct f1_threshold {
    static bool is_classical_batter(int nR, int nS, LL OUT_J) {
        double f1 = (2 * nR + nS) * Para_Radix_LSMapper_build_cost
            + nS * Para_Radix_LSMapper_probe_cost
            + OUT_J * cal_mem_rand_read_write_cost(min(nR, nS))
            - OUT_J * Para_Radix_Flat_Hash_Table_cost;
        // cout << "f1=" << f1 << endl;
        return f1 > 0;
    }
};

struct f2_threshold {
private:
    int* R_counter_array;
    int R_counter_max;
    double base_sparse_hat;
    double base_dense_hat;
    double cur_para_mem_rand_read_ny, cur_para_mem_rand_read_nz, cur_para_mem_rand_read_nR;
    double cur_para_mem_rand_read_and_wirte_nz;

    double TMP_Scache[105];
    int TMP_Sintervel;

public:
    f2_threshold(int* R_JR, int nx, int ny, int nz) {
        cur_para_mem_rand_read_ny = cal_mem_rand_read_cost(ny);
        cur_para_mem_rand_read_nR = cal_mem_rand_read_cost(R_JR[nx]);
        cur_para_mem_rand_read_and_wirte_nz = cal_mem_rand_read_write_cost(nz);
        double dnR = (double)R_JR[nx] / 1e3;
        double dRy = (double)R_JR[nx] / (double)ny;

        base_sparse_hat = 0;
        base_sparse_hat += 2 * (double)nx / 1e3 * Para_Sz_mem_seq_read / nz;
        base_sparse_hat += dnR * Para_Sz_mem_seq_read / nz;
        base_sparse_hat += 2 * dnR * cur_para_mem_rand_read_ny / nz;

        base_dense_hat = dnR * Para_Sz_Bitmap / nz;

        R_counter_array = (int*)malloc(200 * sizeof(int)); //assume 1e20 as max
        if (R_counter_array == NULL) {
            printf("Fail to malloc R_counter_array in f2_threshold.\n");
            exit(-1);
        }
        assert(R_counter_array != NULL);
        memset(R_counter_array, 0, 200 * sizeof(int));

        R_counter_max = 0;
        for (int i = 0; i < nx; i++) {
            int t = val2idx(R_JR[i + 1] - R_JR[i]);
            R_counter_max = max(R_counter_max, t);
            R_counter_array[t]++;
        }

        assert(R_counter_max < 200);
        int* p = (int*)realloc(R_counter_array, (R_counter_max + 1) * sizeof(int));
        if (p == NULL) {
            printf("Fail to realloc R_counter_array f2_threshold.\n");
            free(R_counter_array);
            exit(-1);
        }
        R_counter_array = p;

        for (int i = 0; i < 105; i++) TMP_Scache[i] = -1;
        TMP_Sintervel = max(ny / 100, 1);
    }

    bool is_DenseEC_better(int Sm_cnt, int ny, int outj) {
        int t = Sm_cnt / TMP_Sintervel;
        if (TMP_Scache[t] == -1) TMP_Scache[t] = cal_dense_hat(t * TMP_Sintervel + TMP_Sintervel / 2, ny);
        return TMP_Scache[t] < cal_sparse_hat(outj);
    }

private:
    inline int val2idx(int val) {
        if (val <= 0) return 0;
        int _t = log10(val);
        return _t * 10 + (val / (int)pow(10, _t));
    }

    double cal_sparse_hat(double dout_j) {
        return base_sparse_hat + dout_j / 1e3 * (Para_Sz_mem_seq_read + cur_para_mem_rand_read_and_wirte_nz); // 1e-6 seconds
    }

    double cal_dense_hat(int Sm_cnt, int ny) {
        double p1 = min(Sm_cnt / (double)ny, 1.);
        double _p = (double)Sm_cnt / ((double)ny * (double)ny);
        ny = (ny - 1) / 256 + 1;
        double ret = 0.;
        for (int i = 0, t = 0, tp = 1; i <= R_counter_max; i++, t += tp) {
            double p256 = 1 - pow(1 - min(_p * t, 1.), 256);
            double hat256 = p256 < 1e-5 ? 9e99 : (1 - pow(1 - p256, ny)) / p256 * Para_Rx_Dense;
            double hat1 = (1 - pow(1 - p1, t)) / p1 * Para_Rx_Sparse;
            ret += min(hat256, hat1) * R_counter_array[i];
            if (i != 0 && i % 10 == 0) {
                t -= tp / 2;
                tp *= 10;
                t -= tp / 2;
            }
        }
        return base_dense_hat + (Sm_cnt * Para_Sz_Bitmap + ret) * 1e-3; // 1e-6 seconds
    }
};

struct f3_threshold {
public:
    static int cal_SIMD_method_threshold(int nr, int k) {
        int low = 1;
        int hig = k;
        while (hig - low > 5) {
            int mid = (low + hig) >> 1;
            if (cal_dense_hat(nr, mid, k) < cal_sparse_hat(nr, mid, k)) hig = mid;
            else low = mid;
        }
        for (int i = low; i <= hig; i++) {
            if (cal_dense_hat(nr, i, k) <= cal_sparse_hat(nr, i, k)) return i;
        }
        return hig;
    }

    static bool is_dense_better(int nr, int ns, int k) {
        double dense_hat = cal_dense_hat(nr, ns, k);
        double sparse_hat = cal_sparse_hat(nr, ns, k);
        return dense_hat < sparse_hat;
    }

private:
    static double cal_dense_t_hat(int nr, int ns, int k) {
        double p = (double)nr / (double)k * (double)ns / (double)k;
        p = 1 - pow(1 - p, 256);
        k = (k - 1) / 256 + 1;
        return (1 - pow(1 - p, k)) / p;
    }

    static double cal_dense_hat(int nr, int ns, int k) {
        return cal_dense_t_hat(nr, ns, k) * 1e-3 * Para_Rx_Dense; // 1e-6 seconds
    }

    static double cal_sparse_t_hat(int nr, int ns, int k) {
        double p = (double)ns / (double)k;
        return (1 - pow(1 - p, nr)) / p;
    }

    static double cal_sparse_hat(int nr, int ns, int k) {
        return cal_sparse_t_hat(nr, ns, k) * 1e-3 * Para_Rx_Sparse; // 1e-6 seconds
    }
};
//===============================================


template <typename _Tx, typename _Ty, typename _Tz, typename _H = hash<_Ty>, typename _Tcounter = uint32_t>
struct EstimateJoinCardinality {
public:
    static LL estimate(const myvector<pair<_Tx, _Ty>, _Tcounter> R, const myvector<pair<_Tz, _Ty>, _Tcounter> S) {
        assert(R.size() >= S.size());

        const int sample_bits = cal_sample_bits(S.size());
        const int sample_mask = (1 << sample_bits) - 1;
        _H hash;
        ska::flat_hash_map<_Ty, _Tcounter, _H> hash_map;
        hash_map.reserve(base_sample_size);
        hash_map.set_hash_shr_bits(sample_bits);
        for (int i = 0, _i = S.size(); i < _i; i++) {
            if ((hash(S.data[i].second) & sample_mask) == 0) {
                hash_map[S.data[i].second]++;
            }
        }

        LL result = 0;
        for (int i = 0, _i = R.size(); i < _i; i++) {
            if ((hash(R.data[i].second) & sample_mask) == 0) {
                result += hash_map[R.data[i].second];
            }
        }

        if (result == 0) {
            hash_map.clear();
            hash_map.set_hash_shr_bits(0);
            for (int i = 0; i < 1000; i++) hash_map[R.data[i].second]++;
            for (int i = 0; i < 1000; i++) result += hash_map[S.data[i].second];
            return result * ((double)R.size() / 1e3) * ((double)S.size() / 1e3);
        }

        return result << sample_bits;
    }

private:
    static const int base_sample_size = 10000;
    static int cal_sample_bits(const int N) {
        int sample_size = N, sample_bits = 0;
        while (sample_size > base_sample_size) {
            sample_bits++;
            sample_size /= 2;
        }
        return sample_bits;
    }
};


template <typename _Tx, typename _Ty, typename _Tz, typename _H, typename _Hpair, typename _Tcounter = uint32_t>
struct Radix_hash_join_and_project {
private:
    ska::flat_hash_set<pair<_Tx, _Tz>, _Hpair> project_set;
    myvector<myvector<pair<_Tx, _Ty>>> hashBuckets;
    const int __cur_partition_para = __L2cache_size / max(sizeof(pair<_Tz, _Ty>), sizeof(pair<_Tx, _Ty>)) / 2;
    const int hash_bucket_bits = (int)ceil(log2(__L2cache_size / 2 / __cache_line_size));
    const int hash_bucket_size = 1 << hash_bucket_bits;
    const int hash_bucket_mask = hash_bucket_size - 1;
    int partition_bits, partition_size;
    _Tcounter nR, nS;

public:
    void dojoinproject(const myvector<pair<_Tx, _Ty>, _Tcounter> R,
                            const myvector<pair<_Tz, _Ty>, _Tcounter> S,
                            const uint32_t OUT_hat,
                            myvector<pair<_Tx, _Tz>, _Tcounter>& result) {
        
        _Tcounter* hist, * psumR, * psumS;
        nR = R.size();
        nS = S.size();
        project_set.reserve(OUT_hat);
        partition_bits = (int)ceil(log2(max((int)max(nR, nS), (__cur_partition_para << 1)) / __cur_partition_para));//TODO
        assert(hash_bucket_bits < 30);
        assert(partition_bits < 30);
        partition_size = 1 << partition_bits;

        hist = (_Tcounter*)malloc((partition_size + 1) * sizeof(_Tcounter));
        psumR = (_Tcounter*)malloc((partition_size + 1) * sizeof(_Tcounter));
        psumS = (_Tcounter*)malloc((partition_size + 1) * sizeof(_Tcounter));
        if (hist == NULL || psumR == NULL || psumS == NULL) {
            printf("Fail to malloc hist/psumR/psumS in mapping, it seens that they need %fG memory.\n",
                3 * (LL)(partition_size + 1) * sizeof(_Tcounter) / 1073741824.0);
            exit(-1);
        }

        pair<_Tx, _Ty>* pR = partition(R.data, hist, psumR, nR);
        pair<_Tz, _Ty>* pS = partition(S.data, hist, psumS, nS);


        hashBuckets.fixSize(hash_bucket_size);
        unsigned int init_size = max(1U, __cache_line_size / (unsigned)sizeof(pair<_Tx, _Ty >));
        for (int i = 0; i < hash_bucket_size; i++) {
            hashBuckets.data[i].init(init_size);
        }

        for (int i = 0; i < partition_size; i++) {
            join(pR, pS, psumR[i], psumR[i + 1], psumS[i], psumS[i + 1]);
        }

        free(hist);
        free(psumR);
        free(psumS);
        free(pR);
        free(pS);
        
        result.reserve(project_set.size());
        for (auto i : project_set) {
            result.push_back(i);
        }
    }

private:
    template <typename _Ta, typename _Tb>
    pair<_Ta, _Tb>* partition(const pair<_Ta, _Tb>* V,
                                _Tcounter* hist, _Tcounter* psum,
                                const _Tcounter nV) {

        _H hash;
        const int partition_mask = partition_size - 1;

        memset(hist, 0, partition_size * sizeof(_Tcounter));
        for (_Tcounter i = 0; i < nV; i++) {
            hist[hash(V[i].second) & partition_mask]++;
        }

        psum[0] = 0;
        for (int i = 1; i < partition_size; i++) {
            psum[i] = psum[i - 1] + hist[i - 1];
        }
        psum[partition_size] = psum[partition_size - 1] + hist[partition_size - 1];
        assert(psum[partition_size] == nV);
        memcpy(hist, psum, (partition_size + 1) * sizeof(_Tcounter));

        pair<_Ta, _Tb>* V_p = (pair<_Ta, _Tb>*)malloc(nV * sizeof(pair<_Ta, _Tb>));
        if (V_p == NULL) {
            printf("Fail to malloc V_p in partition_internal_second.\n");
            exit(-1);
        }

        for (_Tcounter i = 0; i < nV; i++) {
            V_p[hist[hash(V[i].second) & partition_mask]++] = V[i];
        }
        return V_p;
    }

    void join(const pair<_Tx, _Ty>* R, const pair<_Tz, _Ty>* S,
        _Tcounter Rbegin, _Tcounter Rend, _Tcounter Sbegin, _Tcounter Send) {
        _H hash;

        for (int i = 0; i < hashBuckets.size(); i++) {
            hashBuckets.data[i].clear();
        }

        for (int i = Rbegin; i < Rend; i++) {
            hashBuckets.data[(hash(R[i].second) >> partition_bits) & hash_bucket_mask].push_back(R[i]);
        }

        for (int i = Sbegin; i < Send; i++) {
            auto& v = hashBuckets.data[(hash(S[i].second) >> partition_bits)& hash_bucket_mask];
            auto v_cnt = v.size();
            for (int j = 0; j < v_cnt; j++) {
                if (S[i].second == v.data[j].second) {
                    project_set.insert({ v.data[j].first,S[i].first });
                }
            }
        }
    }
};


template <typename _Tx, typename _Ty, typename _Tz,
    typename _Hx = hash<_Tx>, typename _Hy = hash<_Ty>, typename _Hz = hash<_Tz>,
    typename _Tcounter = uint32_t>
struct mapping_step {
public:
    static void mapping(const myvector<pair<_Tx, _Ty>, _Tcounter> R,
        const myvector<pair<_Tz, _Ty>, _Tcounter> S,
        ID2VALUE<_Tx, /*_Ty,*/ _Tz>& ID2V,
        myvector<pair<int, int>, _Tcounter>& Rout,
        myvector<pair<int, int>, _Tcounter>& Sout,
        int& xmax, int& ymax, int& zmax) {

        const unsigned int partition_target_size = max((unsigned int)(4 * __L2cache_size / sizeof(pair<_Tx, _Ty>)), __partition_para);
        const int partition_bits = (int)ceil(log2(max(max(R.size(), S.size()), (partition_target_size << 1)) / partition_target_size));//TODO
        assert(partition_bits < 30);
        const int partition_size = 1 << partition_bits;
        const _Tcounter nR = R.size();
        _Tcounter nS = S.size();
        _Tcounter* hist, * psumR, * psumS;
        pair<_Tx, int>* R_t;
        pair<_Tz, int>* S_t;
        _Hx hasherx;
        _Hy hashery;
        _Hz hasherz;

        hist = (_Tcounter*)malloc((partition_size + 1) * sizeof(_Tcounter));
        psumR = (_Tcounter*)malloc((partition_size + 1) * sizeof(_Tcounter));
        psumS = (_Tcounter*)malloc((partition_size + 1) * sizeof(_Tcounter));
        if (hist == NULL || psumR == NULL || psumS == NULL) {
            printf("Fail to malloc hist/psumR/psumS in mapping, it seens that they need %fG memory.\n",
                3 * (LL)(partition_size + 1) * sizeof(_Tcounter) / 1073741824.0);
            exit(-1);
        }

        //mapping R.y & S.y
        pair<_Tx, _Ty>* pR = partition_y(R.data, hist, psumR, nR, partition_size, hashery);
        pair<_Tz, _Ty>* pS = partition_y(S.data, hist, psumS, nS, partition_size, hashery);
        mapping_y(pR, pS, psumR, psumS, nR, nS, partition_bits, hashery, /*ID2V.y,*/ ymax, R_t, S_t);

        if (nR == 0 || nS == 0) {
            Rout.cur_num = Sout.cur_num = 0;
            return;
        }

        //mapping R.x
        partition_xz(R_t, hist, psumR, nR, partition_size, hasherx);
        mapping_xz(R_t, psumR, nR, partition_bits, hasherx, ID2V.x, xmax, Rout.data);
        Rout.cur_num = Rout.max_num = nR;

        //mapping S.z
        partition_xz(S_t, hist, psumS, nS, partition_size, hasherz);
        mapping_xz(S_t, psumS, nS, partition_bits, hasherz, ID2V.z, zmax, Sout.data);
        Sout.cur_num = Sout.max_num = nS;

        free(pR);
        free(pS);
        free(R_t);
        free(S_t);
        free(hist);
        free(psumR);
        free(psumS);
    }

    static void mapping(const myvector<pair<int, int>, _Tcounter> R,
        const myvector<pair<int, int>, _Tcounter> S,
        ID2VALUE<int, /*_Ty,*/ int>& ID2V,
        myvector<pair<int, int>, _Tcounter>& Rout,
        myvector<pair<int, int>, _Tcounter>& Sout,
        int& xmax, int& ymax, int& zmax,
        int hintxyz) {

        if (hintxyz * 2 * 4 + 16384 < __L3cache_size_by_int) {

            Rout.reserve(R.size());
            Sout.reserve(S.size());
            ID2V.x = (int*)malloc(hintxyz * sizeof(int));
            ID2V.z = (int*)malloc(hintxyz * sizeof(int));
            ska::flat_hash_map<int, int, myHasher> m;
            m.reserve(hintxyz);

            int cnt = -1;
            for (int i = 0, _i = R.size(); i < _i; i++) {
                auto tmp = m.find(R.data[i].second);
                if (tmp == m.end()) {
                    cnt++;
                    m[R.data[i].second] = cnt;
                    Rout.data[i].second = cnt;
                }
                else {
                    Rout.data[i].second = tmp->second;
                }
            }
            int t = 0;
            for (int i = 0, _i = S.size(); i < _i; i++) {
                auto tmp = m.find(S.data[i].second);
                if (tmp != m.end()) {
                    Sout.data[t] = { S.data[i].first,tmp->second };
                    t++;
                }
            }
            Rout.size() = R.size();
            Sout.size() = t;
            ymax = cnt;

            cnt = -1;
            m.clear();
            for (int i = 0, _i = Rout.size(); i < _i; i++) {
                auto tmp = m.find(R.data[i].first);
                if (tmp == m.end()) {
                    cnt++;
                    m[R.data[i].first] = cnt;
                    Rout.data[i].first = cnt;
                }
                else {
                    Rout.data[i].first = tmp->second;
                }
            }
            xmax = cnt;
            for (auto& i : m) ID2V.x[i.second] = i.first;

            cnt = -1;
            m.clear();
            for (int i = 0, _i = Sout.size(); i < _i; i++) {
                auto tmp = m.find(Sout.data[i].first);
                if (tmp == m.end()) {
                    cnt++;
                    m[Sout.data[i].first] = cnt;
                    Sout.data[i].first = cnt;
                }
                else {
                    Sout.data[i].first = tmp->second;
                }
            }
            zmax = cnt;
            for (auto& i : m) ID2V.z[i.second] = i.first;

            //ID2V.x = (int*)realloc(ID2V.x, (xmax + 1) * sizeof(int));
            //ID2V.z = (int*)realloc(ID2V.z, (zmax + 1) * sizeof(int));

            return;
        }
    
        mapping(R, S, ID2V, Rout, Sout, xmax, ymax, zmax);
    }

private:
    template <typename _Ta, typename _Tb, typename _H>
    static void partition_xz(pair<_Ta, _Tb>* V,
        _Tcounter* hist, _Tcounter* psum,
        const _Tcounter nV, const int partition_size, const _H hasher) {

        const int partition_mask = partition_size - 1;

        memset(hist, 0, partition_size * sizeof(_Tcounter));
        for (_Tcounter i = 0; i < nV; i++) {
            hist[hasher(V[i].first) & partition_mask]++;
        }

        psum[0] = 0;
        for (int i = 1; i < partition_size; i++) {
            psum[i] = psum[i - 1] + hist[i - 1];
        }
        psum[partition_size] = psum[partition_size - 1] + hist[partition_size - 1];
        assert(psum[partition_size] == nV);
        memcpy(hist, psum, (partition_size + 1) * sizeof(_Tcounter));

        pair<_Ta, _Tb>* V_p = (pair<_Ta, _Tb>*)malloc(nV * sizeof(pair<_Ta, _Tb>));
        if (V_p == NULL) {
            printf("Fail to malloc V_p in partition_xz.\n");
            exit(-1);
        }

        for (_Tcounter i = 0; i < nV; i++) {
            V_p[hist[hasher(V[i].first) & partition_mask]++] = V[i];
        }

        memcpy(V, V_p, nV * sizeof(pair<_Ta, _Tb>));
        free(V_p);
    }

    template <typename _Ta, typename _Tb,typename _H >
    static pair<_Ta, _Tb>* partition_y(const pair<_Ta, _Tb>* V,
        _Tcounter* hist, _Tcounter* psum,
        const _Tcounter nV, const int partition_size, const _H hasher) {

        const int partition_mask = partition_size - 1;

        memset(hist, 0, partition_size * sizeof(_Tcounter));
        for (_Tcounter i = 0; i < nV; i++) {
            hist[hasher(V[i].second) & partition_mask]++;
        }

        psum[0] = 0;
        for (int i = 1; i < partition_size; i++) {
            psum[i] = psum[i - 1] + hist[i - 1];
        }
        psum[partition_size] = psum[partition_size - 1] + hist[partition_size - 1];
        assert(psum[partition_size] == nV);
        memcpy(hist, psum, (partition_size + 1) * sizeof(_Tcounter));

        pair<_Ta, _Tb>* V_p = (pair<_Ta, _Tb>*)malloc(nV * sizeof(pair<_Ta, _Tb>));
        if (V_p == NULL) {
            printf("Fail to malloc V_p in partition_y.\n");
            exit(-1);
        }

        for (_Tcounter i = 0; i < nV; i++) {
            V_p[hist[hasher(V[i].second) & partition_mask]++] = V[i];
        }

        return V_p;
    }

    template <typename _Ta, typename _H>
    static void mapping_xz(const pair<_Ta, int>* V, const _Tcounter* psum,
        const _Tcounter nV, const int partition_bits, const _H hasher,
        _Ta*& id2v, int& the_max, pair<int, int>*& Vout) {

        const int partition_size = 1 << partition_bits;
        LS_mapper<_Ta, _H> hashmap(psum[1] - psum[0], partition_bits, 16384 / 4);
        id2v = (_Ta*)malloc(nV * sizeof(_Ta));
        Vout = (pair<int, int>*)malloc(nV * sizeof(pair<int, int>));
        if (id2v == NULL || Vout == NULL) {
            printf("Fail to malloc id2v/Vout in mapping_xz.\n");
            exit(-1);
        }

        for (int i = 0; i < partition_size; i++) {
            hashmap.clear();
            for (_Tcounter j = psum[i], _j = psum[i + 1]; j < _j; j++) {
                Vout[j].first = hashmap.insert(V[j].first, id2v);
                Vout[j].second = V[j].second;
                //Vout[j] = { hashmap.insert(V[j].first , id2v) ,V[j].second };
            }
        }
        id2v = (_Ta*)realloc(id2v, hashmap.size() * sizeof(_Ta));
        if (id2v == NULL) {
            printf("Fail to realloc id2v in mapping_xz.\n");
            exit(-1);
        }
        the_max = hashmap.size() - 1;
    }

    template <typename _Ta, typename _Tb, typename _Tc, typename _H>
    static void mapping_y(const pair<_Ta, _Tb>* V1, const pair<_Tc, _Tb>* V2,
        const _Tcounter* psum1, const _Tcounter* psum2,
        const _Tcounter n1, _Tcounter& n2,
        const int partition_bits, const _H hasher,
        /*IMDBstring*& id2v,*/ int& the_max,
        pair<_Ta, int>*& V1out, pair<_Tc, int>*& V2out) {

        const int partition_size = 1 << partition_bits;
        LS_mapper<_Tb, _H> hashmap(psum1[1] - psum1[0], partition_bits, 16384 / 4);
        _Tcounter v2cnt = 0;
        //id2v = (IMDBstring*)malloc(n1 * sizeof(IMDBstring));
        V1out = (pair<_Ta, int>*)malloc(n1 * sizeof(pair<_Ta, int>));
        V2out = (pair<_Tc, int>*)malloc(n2 * sizeof(pair<_Tc, int>));
        if (/*id2v == NULL || */V1out == NULL || V2out == NULL) {
            printf("Fail to malloc V1out/V2out in mapping_y.\n");
            exit(-1);
        }

        for (int i = 0; i < partition_size; i++) {
            hashmap.clear();
            for (_Tcounter j = psum1[i], _j = psum1[i + 1]; j < _j; j++) {
                V1out[j].first = V1[j].first;
                V1out[j].second = hashmap.insert(V1[j].second /*, id2v*/);
                //V1out[j] = { V1[j].first,hashmap.insert(V1[j].second /*, id2v*/) };
            }

            for (_Tcounter j = psum2[i], _j = psum2[i + 1]; j < _j; j++) {
                int t = hashmap.find(V2[j].second);
                if (t != -1) {
                    V2out[v2cnt].first = V2[j].first;
                    V2out[v2cnt].second = t;
                    v2cnt++;
                    //V2out[v2cnt++] = { V2[j].first,t };
                }
            }
        }

        //id2v = (IMDBstring*)realloc(id2v, hashmap.size() * sizeof(IMDBstring));
        //if (id2v == NULL) exit(-1);
        the_max = hashmap.size() - 1;
        n2 = v2cnt;
    }

};


template <typename _Tx, /*typename _Ty,*/ typename _Tz, typename _Tcounter = uint32_t>
struct Hybrid_solution {
public:
    static void dojoinproject(CSR<>& R_all, S_DENSE& S_dense, CSR<>& S_sparse, ID2VALUE<_Tx,_Tz>& ID2V,
                                const int nx, const int ny, const int nz,
                                myvector<pair<_Tx,_Tz>, _Tcounter>& result, size_t reserve_size/*?*/) {
        result.reserve(reserve_size);

        if (S_dense.size > 0) {
            //DenseEC
            int dense_bitvector_len = S_dense.dense_bitvector_len;
            myBitVector<int> dense_x_bitvector(dense_bitvector_len);
            int* origin_z_val = (int*)malloc(S_dense.size * sizeof(int));
            int* f3_threshold_cache= (int*)malloc(ny * sizeof(int));
            if (origin_z_val == NULL || f3_threshold_cache == NULL) {
                printf("Fail to malloc origin_z_val/f3_threshold_cache in SPMX.\n");
                exit(-1);
            }
            for (int i = 0; i < S_dense.size; i++) {
                origin_z_val[i] = ID2V.z[S_dense.id2v[i]];
            }
            memset(f3_threshold_cache,-1,ny * sizeof(int));
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
                                result.push_back({ origin_x_val,origin_z_val[j] });
                                break;
                            }
                        }
                    }
                    else {
                        //non-SIMD method
                        for (int t = R_all.JR[i]; t < R_all.JR[i + 1]; t++) {
                            if (S_dense.data[j].checkIfTrue(R_all.IC[t])) {
                                result.push_back({ origin_x_val,origin_z_val[j] });
                                break;
                            }
                        }
                    }
                }
            }
            free(origin_z_val);
            free(f3_threshold_cache);
            free(dense_x_bitvector._Array);
        }

        if (S_sparse.JR != NULL) {
            int *SPAw = (int*)malloc(nz * sizeof(int));
            if (SPAw == NULL) {
                printf("Fail to malloc SPAw in SPMX.\n");
                exit(-1);
            }
            memset(SPAw, -1, nz * sizeof(int));

            for (int cur_x = 0; cur_x < nx; cur_x++) {
                int origin_x_val = ID2V.x[cur_x];
                for (int i = R_all.JR[cur_x], _i = R_all.JR[cur_x + 1]; i < _i; i++) {
                    int cur_y = R_all.IC[i];
                    for (int j = S_sparse.JR[cur_y], _j = S_sparse.JR[cur_y + 1]; j < _j; j++) {
                        if (SPAw[S_sparse.IC[j]] != cur_x) {
                            SPAw[S_sparse.IC[j]] = cur_x;
                            result.push_back({ origin_x_val,ID2V.z[S_sparse.IC[j]] });
                        }
                    }
                }
            }
            free(SPAw);
        }
    }
};


//===============================================
/*
* DIM3
*/
template <typename _Tx, typename _Ty, typename _Tz,
    typename _Hx = hash<_Tx>, typename _Hy = hash<_Ty>, typename _Hz = hash<_Tz>,
    typename _Hxz = hash<pair<_Tx, _Tz>>, typename _Tcounter = uint32_t>
struct DIM3 {
public:
    static void doJoinProject(const myvector<pair<_Tx, _Ty>, _Tcounter> _R,
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

            Radix_hash_join_and_project<_Tx, _Ty, _Tz, _Hy, _Hxz> classical_solution;
            classical_solution.dojoinproject(R, S, (uint32_t)min(OUT_J_hat, (LL)1e8), result);

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

        divide_R(Rm, n, R_all);
        divide_S(Sm, n, k, m, R_all, S_dense, S_sparse);

        free(Rm.data);
        free(Sm.data);

        Hybrid_solution<_Tx, _Tz, _Tcounter>::dojoinproject(R_all, S_dense, S_sparse, ID2V, n, k, m, result, min(OUT_J_hat*3, (LL)n * m));

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

private:
    static void divide_R(const myvector<pair<int, int>> V, const int nx, CSR<>& R_all) {
        int nV = V.size();
        R_all.JR = (int*)malloc((nx + 1) * sizeof(int));
        int* JRt = (int*)malloc((nx + 1) * sizeof(int));
        R_all.IC = (int*)malloc(nV * sizeof(int));
        if (R_all.JR == NULL || JRt == NULL || R_all.IC == NULL) {
            printf("Fail to malloc R_all.JR/JRt/R_all.IC in divide_R.\n");
            exit(-1);
        }
        memset(JRt, 0, (nx + 1) * sizeof(int));

        for (int i = 0; i < nV; i++) {
            JRt[V.data[i].first]++;
        }

        R_all.JR[0] = 0;
        for (int i = 1; i < nx; i++) {
            R_all.JR[i] = R_all.JR[i - 1] + JRt[i - 1];
        }
        R_all.JR[nx] = R_all.JR[nx - 1] + JRt[nx - 1];
        assert(R_all.JR[nx] == nV);
        memcpy(JRt, R_all.JR, (nx + 1) * sizeof(int));

        for (int i = 0; i < nV; i++) {
            R_all.IC[JRt[V.data[i].first]++] = V.data[i].second;
        }

        free(JRt);
    }

    static void divide_S(const myvector<pair<int, int>> V,
                            const int nx, const int ny, const int nz,
                            const CSR<> R_all,
                            S_DENSE& S_dense, CSR<>& S_sparse) {
        const int nV = V.size();
        int* zCNT = NULL, * z2denseid = NULL, * JRt = NULL;
        zCNT = (int*)malloc(nz * sizeof(int));
        z2denseid = (int*)malloc(nz * sizeof(int));
        if (zCNT == NULL || z2denseid == NULL) {
            printf("Fail to malloc zCNT/z2denseid in divide_S.\n");
            exit(-1);
        }
        memset(zCNT, 0, nz * sizeof(int));
        memset(z2denseid, -1, nz * sizeof(int));

        for (int i = 0; i < nV; i++) {
            zCNT[V.data[i].first]++;
        }

        int* Ry_counter = (int*)malloc(ny * sizeof(int));
        memset(Ry_counter, 0, ny * sizeof(int));
        for (int i = 0, _t = R_all.JR[nx]; i < _t; i++) {
            Ry_counter[R_all.IC[i]]++;
        }
        int* OUTJ_CNT = (int*)malloc(nz * sizeof(int));
        memset(OUTJ_CNT, 0, nz * sizeof(int));
        for (int i = 0; i < nV; i++) {
            OUTJ_CNT[V.data[i].first] += Ry_counter[V.data[i].second];
        }

        int densecnt = 0;
        assert(R_all.JR != NULL);
        f2_threshold f2 = f2_threshold(R_all.JR, nx, ny, nz);
        for (int i = 0; i < nz; i++) {
            if (f2.is_DenseEC_better(zCNT[i], ny, OUTJ_CNT[i])) {
                z2denseid[i] = densecnt;
                densecnt++;
            }
        }

        S_dense.size = densecnt;
        S_dense.dense_bitvector_len = (ny / myBitVector<>::_Bitsperword) + (ny % myBitVector<>::_Bitsperword == 0 ? 0 : 1);
        S_dense.dense_bitvector_len = (S_dense.dense_bitvector_len % 8 == 0 ?
            S_dense.dense_bitvector_len : (S_dense.dense_bitvector_len / 8 + 1) * 8);
        assert(S_dense.dense_bitvector_len % 8 == 0);
        
        // cout << "densecnt=" << densecnt << endl;
        
        assert(densecnt <= nz);
        if (densecnt == nz) {
            S_dense.id2v = (int*)malloc(densecnt * sizeof(int));
            S_dense.id2cnt = (int*)malloc(densecnt * sizeof(int));
            if (S_dense.id2v == NULL || S_dense.id2cnt == NULL) {
                printf("Fail to malloc S_dense.id2v in divide_S.\n");
                exit(-1);
            }
            for (int i = 0; i < densecnt; i++) {
                S_dense.id2v[i] = i;
                S_dense.id2cnt[i] = zCNT[i/*S_dense.id2v[i]==i*/];
            }

            S_dense.data = (myBitVector<>*)malloc(densecnt * sizeof(myBitVector<>));
            if (S_dense.data == NULL) {
                printf("Fail to malloc S_dense.data in divide_S.\n");
                exit(-1);
            }
            for (int i = 0; i < densecnt; i++) {
                S_dense.data[i] = myBitVector<>(S_dense.dense_bitvector_len);
            }
            for (int i = 0; i < nV; i++) {
                S_dense.data[V.data[i].first].setTrue(V.data[i].second);
            }

            S_sparse.JR = NULL;
            S_sparse.IC = NULL;
        }
        else if (densecnt > 0) {
            S_dense.id2v = (int*)malloc(densecnt * sizeof(int));
            S_dense.id2cnt = (int*)malloc(densecnt * sizeof(int));
            if (S_dense.id2v == NULL || S_dense.id2cnt == NULL) {
                printf("Fail to malloc S_dense.id2v/id2cnt in divide_S.\n");
                exit(-1);
            }
            for (int i = 0; i < nz; i++) {
                if (z2denseid[i] != -1) {
                    S_dense.id2v[z2denseid[i]] = i;
                    S_dense.id2cnt[z2denseid[i]] = zCNT[i];
                }
            }

            S_dense.data = (myBitVector<>*)malloc(densecnt * sizeof(myBitVector<>));
            if (S_dense.data == NULL) {
                printf("Fail to malloc S_dense.data in divide_S.\n");
                exit(-1);
            }
            for (int i = 0; i < densecnt; i++) {
                S_dense.data[i] = myBitVector<>(S_dense.dense_bitvector_len);
            }
            S_sparse.JR = (int*)malloc((ny + 1) * sizeof(int));
            JRt = (int*)malloc((ny + 1) * sizeof(int));
            if (S_sparse.JR == NULL || JRt == NULL) {
                printf("Fail to malloc S_sparse.JR/JRt in divide_S.\n");
                exit(-1);
            }
            memset(JRt, 0, (ny + 1) * sizeof(int));
            int nsparse = 0;
            for (int i = 0; i < nV; i++) {
                if (z2denseid[V.data[i].first] != -1) {
                    //dense
                    S_dense.data[z2denseid[V.data[i].first]].setTrue(V.data[i].second);
                }
                else {
                    //sparse
                    JRt[V.data[i].second]++;
                    V.data[nsparse] = V.data[i];
                    nsparse++;
                }
            }

            S_sparse.JR[0] = 0;
            for (int i = 1; i < ny; i++) {
                S_sparse.JR[i] = S_sparse.JR[i - 1] + JRt[i - 1];
            }
            S_sparse.JR[ny] = S_sparse.JR[ny - 1] + JRt[ny - 1];
            assert(S_sparse.JR[ny] == nsparse);
            memcpy(JRt, S_sparse.JR, (ny + 1) * sizeof(int));

            S_sparse.IC = (int*)malloc(nsparse * sizeof(int));
            if (S_sparse.IC == NULL) {
                printf("Fail to malloc S_sparse.IC in divide_S.\n");
                exit(-1);
            }

            for (int i = 0; i < nsparse; i++) {
                S_sparse.IC[JRt[V.data[i].second]++] = V.data[i].first;
            }
        }
        else {
            S_dense.id2v = NULL;
            S_dense.id2cnt = NULL;
            S_dense.data = NULL;

            S_sparse.JR = (int*)malloc((ny + 1) * sizeof(int));
            S_sparse.IC = (int*)malloc(nV * sizeof(int));
            JRt = (int*)malloc((ny + 1) * sizeof(int));
            if (S_sparse.JR == NULL || S_sparse.IC == NULL || JRt == NULL) {
                printf("Fail to malloc S_sparse.JR/JRt/S_sparse.JR in divide_S.\n");
                exit(-1);
            }
            memset(JRt, 0, (ny + 1) * sizeof(int));

            for (int i = 0; i < nV; i++) {
                JRt[V.data[i].second]++;
            }

            S_sparse.JR[0] = 0;
            for (int i = 1; i < ny; i++) {
                S_sparse.JR[i] = S_sparse.JR[i - 1] + JRt[i - 1];
            }
            S_sparse.JR[ny] = S_sparse.JR[ny - 1] + JRt[ny - 1];
            assert(S_sparse.JR[ny] == nV);
            memcpy(JRt, S_sparse.JR, (ny + 1) * sizeof(int));

            for (int i = 0; i < nV; i++) {
                S_sparse.IC[JRt[V.data[i].second]++] = V.data[i].first;
            }
        }

        if (JRt) free(JRt);
        free(zCNT);
        free(z2denseid);
    }
};

