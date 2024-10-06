// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// DIM3 is licensed under Mulan PSL v2.

#pragma once
#include"myBitVector.hpp"

#define FORCE_INLINE inline __attribute__((always_inline))
#define ROTL32(x,y)	generic_rotl32(x,y)
#define ROTR32(x,y)	generic_rotr32(x,y)

using LL = long long;
using uLL = unsigned long long;

struct S_DENSE {
    myBitVector<>* data;
	int* id2v;
    int* id2cnt;
	int size;
    int dense_bitvector_len;
};

template <typename _Tv = int>
struct CSR {
	int* JR, * IC;
    _Tv* VAL;
};

template <typename _Tx, /*typename _Ty,*/ typename _Tz>
struct ID2VALUE
{
    _Tx* x;
    //_Ty* y;
    _Tz* z;
};


static FORCE_INLINE
uint32_t generic_rotr32(const uint32_t x, const unsigned bits)
{
    const unsigned n = bits % 32;
    return (x >> n) | (x << (32 - n));
}

static FORCE_INLINE
uint32_t generic_rotl32(const uint32_t x, const unsigned bits)
{
    const unsigned n = bits % 32;
    return (x << n) | (x >> (32 - n));
}



class IntPairHash
{
public:
    FORCE_INLINE uint32_t operator()(const pair<int, int>& key) const {
        uint32_t h1 = key.first;
        uint32_t k1 = key.second;

        k1 *= 0xcc9e2d51;
        k1 = ROTL32(k1, 15);
        k1 *= 0x1b873593;
        h1 ^= k1;
        h1 = ROTL32(h1, 13);
        h1 = h1 * 5 + 0xe6546b64;

        return h1;
    }
};

template <unsigned len>
struct IntArray {
    int data[len];
    FORCE_INLINE IntArray<len>& operator=(const IntArray<len>& new_val) {
        memcpy(this->data, &new_val.data, len * sizeof(int));
        return *this;
    }
    FORCE_INLINE bool operator==(const IntArray<len>& another) const {
        return memcmp(this->data, another.data, len * sizeof(int)) == 0;
    }
    FORCE_INLINE bool operator<(const IntArray<len>& another) const {
        return memcmp(this->data, another.data, len * sizeof(int)) == -1;
    }
    FORCE_INLINE int operator[](const int i) const {
        return this->data[i];
    }
    FORCE_INLINE int& operator[](const int i) {
        return this->data[i];
    }
};

struct myHasher {
    FORCE_INLINE uint32_t operator()(const int key) const {
        return key;
    }
    template <unsigned len>
    FORCE_INLINE uint32_t operator()(const IntArray<len> key) const {
        uint32_t h = cal_hash(key, 0);
        h = fmix32(h);
        return h;
    }
    template <typename _Ta, typename _Tb>
    FORCE_INLINE uint32_t operator()(const pair<_Ta, _Tb> key) const {
        uint32_t h = cal_hash(key.first, 0);
        h = cal_hash(key.second, h);
        h = fmix32(h);
        return h;
    }

private:
    FORCE_INLINE uint32_t cal_hash(const int key, uint32_t h) const {
        h ^= ROTL32(key, 15);
        h = ROTL32(h, 13);
        return h;
    }
    template <unsigned len>
    FORCE_INLINE uint32_t cal_hash(const IntArray<len> key, uint32_t h) const {
        for (int i = 0; i<len; i++)
        {
            h ^= ROTL32(key[i], 15);
            h = ROTL32(h, 13);
        }
        return h;
    }
    FORCE_INLINE uint32_t fmix32(uint32_t h) const
    {
        h ^= h >> 16;
        h *= 0x85ebca6b;
        h ^= h >> 13;
        h *= 0xc2b2ae35;
        h ^= h >> 16;
        return h;
    }
};
