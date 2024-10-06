// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// DIM3 is licensed under Mulan PSL v2.

#pragma once

template<typename _T = int>
struct myBitVector
{
public:
    static constexpr int _Bitsperword = CHAR_BIT * sizeof(_T);
    static constexpr int _Posshr = 3 + (sizeof(_T) == 4 ? 2 : (sizeof(_T) == 2 ? 1 : 0));//TODO:Only support 1/2/4 byte(s)
    static constexpr int _Posmask = _Bitsperword - 1;
    //static_assert(sizeof(_T) == 1 || sizeof(_T) == 2 || sizeof(_T) == 4);

    _T* _Array = NULL;

    myBitVector() {
        _Array = NULL;
    }
    myBitVector(int n) {
        _Array = (_T*)malloc(n * sizeof(_T));
        if (_Array == NULL) {
            printf("Failed to malloc an BitVector of %fG\n", n * sizeof(_T) / 1073741824.0);
            exit(-1);
        }
        memset(_Array, 0, n * sizeof(_T));
    }

    inline void setTrue(int pos) {
        _Array[pos >> _Posshr] |= _T{ 1 } << (pos & _Posmask);
    }

    inline bool checkIfTrue(int pos) {
        return _Array[pos >> _Posshr] & (_T{ 1 } << (pos & _Posmask));
    }

    void clear(int n) {
        memset(_Array, 0, n * sizeof(_T));
    }
};