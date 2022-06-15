#pragma once

#include <bits/stdc++.h>
#include "flat_hash_map.hpp"

using namespace std;

template <typename  T = int, typename H = hash<T>, typename _Tcounter = int>
struct LS_mapper {
private:
    H hasher;
    ska::flat_hash_map<T, _Tcounter, H> cache_hash_map;//TODO: hasher
    pair<T, _Tcounter>* large_map;
    _Tcounter large_map_size;
    _Tcounter large_map_mask;
    _Tcounter cnt = 0;
    int hash_shr_bits = 0;
public:
    LS_mapper(_Tcounter n, int _the_shr_bits = 0, int cache_n = 16384) {
        large_map_size = 1 << ((_Tcounter)log2(n) + 2);
        large_map_mask = large_map_size - 1;
        large_map_size += 3; //linear probe 4
        large_map = (pair<T, _Tcounter>*)malloc(large_map_size * sizeof(pair<T, _Tcounter>));
        if (large_map == NULL) exit(-1);
        memset(large_map, -1, large_map_size * sizeof(pair<T, _Tcounter>));
        cache_hash_map.reserve(cache_n);
        cache_hash_map.set_hash_shr_bits(_the_shr_bits);
        hash_shr_bits = _the_shr_bits;
    }
    ~LS_mapper() {
        free(large_map);
    }

    void clear() {
        memset(large_map, -1, large_map_size * sizeof(pair<T, _Tcounter>));
        cache_hash_map.clear();
    }

    void reset_cnt() {
        cnt = 0;
    }

    void __set_cnt(_Tcounter val) {
        cnt = val;
    }

    _Tcounter size() {
        return cnt;
    }

    int cache_map_bucket_count() {
        return cache_hash_map.bucket_count();
    }

    _Tcounter insert(T origin_val, T* recorder) {
        _Tcounter ret = insert(origin_val);
        recorder[ret] = origin_val;
        return ret;
    }

    _Tcounter insert(T origin_val) {
        _Tcounter base = (hasher(origin_val) >> hash_shr_bits) & large_map_mask;

        auto& lm0 = large_map[base];
        if (lm0.first == origin_val && lm0.second != -1) {
            return lm0.second;
        }
        else if (lm0.second == -1) {
            lm0 = { origin_val,cnt };
            return cnt++;
        }

        base++;
        auto& lm1 = large_map[base];
        if (lm1.first == origin_val && lm1.second != -1) {
            return lm1.second;
        }
        else if (lm1.second == -1) {
            lm1 = { origin_val,cnt };
            return cnt++;
        }

        base++;
        auto& lm2 = large_map[base];
        if (lm2.first == origin_val && lm2.second != -1) {
            return lm2.second;
        }
        else if (lm2.second == -1) {
            lm2 = { origin_val,cnt };
            return cnt++;
        }

        base++;
        auto& lm3 = large_map[base];
        if (lm3.first == origin_val && lm3.second != -1) {
            return lm3.second;
        }
        else if (lm3.second == -1) {
            lm3 = { origin_val,cnt };
            return cnt++;
        }

        auto tmp = cache_hash_map.find(origin_val);
        if (tmp == cache_hash_map.end()) {
            cache_hash_map[origin_val] = cnt;
            return cnt++;
        }
        else {
            return tmp->second;
        }
    }

    _Tcounter find(T origin_val) {
        _Tcounter base = (hasher(origin_val) >> hash_shr_bits) & large_map_mask;

        auto& lm0 = large_map[base];
        if (lm0.first == origin_val) {
            return lm0.second;
        }
        else if (lm0.second == -1) {
            return -1;
        }

        base++;
        auto& lm1 = large_map[base];
        if (lm1.first == origin_val) {
            return lm1.second;
        }
        else if (lm1.second == -1) {
            return -1;
        }

        base++;
        auto& lm2 = large_map[base];
        if (lm2.first == origin_val) {
            return lm2.second;
        }
        else if (lm2.second == -1) {
            return -1;
        }

        base++;
        auto& lm3 = large_map[base];
        if (lm3.first == origin_val) {
            return lm3.second;
        }
        else if (lm3.second == -1) {
            return -1;
        }

        auto tmp = cache_hash_map.find(origin_val);
        if (tmp != cache_hash_map.end()) {
            return tmp->second;
        }

        return -1;
    }
};