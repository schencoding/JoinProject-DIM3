// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// DIM3 is licensed under Mulan PSL v2.

#pragma once


//=========↓↓↓============
// copy your server parameters here!!
const unsigned int __L1cache_size = 32768;
const unsigned int __L2cache_size = 262144;
const unsigned int __L3cache_size = 12582912; 
const unsigned int __L1cache_size_by_int = __L1cache_size / sizeof(int);
const unsigned int __L2cache_size_by_int = __L2cache_size / sizeof(int);
const unsigned int __L3cache_size_by_int = __L3cache_size / sizeof(int);
const unsigned int __cache_line_size = 64;
const unsigned int __partition_para = 10000;


const double Para_Radix_Flat_Hash_Table_cost = 38.006142; // 1e-9 seconds
const double Para_Radix_LSMapper_build_cost = 21.381317; // 1e-9 seconds
const double Para_Radix_LSMapper_probe_cost = 17.218456; // 1e-9 seconds

const double Para_Rx_Dense = 1.956135; // 1e-9 seconds
const double Para_Rx_Sparse = 4.879672; // 1e-9 seconds
const double Para_Sz_Bitmap = 3.484022; // 1e-9 seconds

const double Para_Sz_mem_seq_read = 0.121759; // 1e-9 seconds
const double Para_Sz_mem_rand_read_1k = 0.499693; // 1e-9 seconds
const double Para_Sz_mem_rand_read_write_1k = 0.443129; // 1e-9 seconds
const double Para_Sz_mem_rand_read_2L2 = 1.005135; // 1e-9 seconds
const double Para_Sz_mem_rand_read_write_2L2 = 1.324415; // 1e-9 seconds
const double Para_Sz_mem_rand_read_05L3 = 1.552501; // 1e-9 seconds
const double Para_Sz_mem_rand_read_write_05L3 = 1.914639; // 1e-9 seconds
const double Para_Sz_mem_rand_read_15L3 = 3.750952; // 1e-9 seconds
const double Para_Sz_mem_rand_read_write_15L3 = 5.359059; // 1e-9 seconds
const double Para_Sz_mem_rand_read_3L3 = 5.417254; // 1e-9 seconds
const double Para_Sz_mem_rand_read_write_3L3 = 7.911922; // 1e-9 seconds
const double Para_Sz_mem_rand_read_10L3 = 7.198486; // 1e-9 seconds
const double Para_Sz_mem_rand_read_write_10L3 = 11.245022; // 1e-9 seconds
//=========↑↑↑============

static double cal_mem_rand_read_cost(int array_size) {
    unsigned int x1 = 1000;
    unsigned int x2 = 2 * __L2cache_size_by_int;
    unsigned int x3 = 0.5 * __L3cache_size_by_int;
    unsigned int x4 = 1.5 * __L3cache_size_by_int;
    unsigned int x5 = 3 * __L3cache_size_by_int;
    unsigned int x6 = 10 * __L3cache_size_by_int;
    if (array_size <= x1) return Para_Sz_mem_rand_read_1k;
    else if (array_size <= x2) {
        return (double)(array_size - x1) / (double)(x2 - x1) * (Para_Sz_mem_rand_read_2L2 - Para_Sz_mem_rand_read_1k) + Para_Sz_mem_rand_read_1k;
    }
    else if (array_size <= x3) {
        return (double)(array_size - x2) / (double)(x3 - x2) * (Para_Sz_mem_rand_read_05L3 - Para_Sz_mem_rand_read_2L2) + Para_Sz_mem_rand_read_2L2;
    }
    else if (array_size <= x4) {
        return (double)(array_size - x3) / (double)(x4 - x3) * (Para_Sz_mem_rand_read_15L3 - Para_Sz_mem_rand_read_05L3) + Para_Sz_mem_rand_read_05L3;
    }
    else if (array_size <= x5) {
        return (double)(array_size - x4) / (double)(x5 - x4) * (Para_Sz_mem_rand_read_3L3 - Para_Sz_mem_rand_read_15L3) + Para_Sz_mem_rand_read_15L3;
    }
    else if (array_size <= x6) {
        return (double)(array_size - x5) / (double)(x6 - x5) * (Para_Sz_mem_rand_read_10L3 - Para_Sz_mem_rand_read_3L3) + Para_Sz_mem_rand_read_3L3;
    }
    else {
        return Para_Sz_mem_rand_read_10L3;
    }
}

static double cal_mem_rand_read_write_cost(int array_size) {
    unsigned int x1 = 1000;
    unsigned int x2 = 2 * __L2cache_size_by_int;
    unsigned int x3 = 0.5 * __L3cache_size_by_int;
    unsigned int x4 = 1.5 * __L3cache_size_by_int;
    unsigned int x5 = 3 * __L3cache_size_by_int;
    unsigned int x6 = 10 * __L3cache_size_by_int;
    if (array_size <= x1) return Para_Sz_mem_rand_read_write_1k;
    else if (array_size <= x2) {
        return (double)(array_size - x1) / (double)(x2 - x1) * (Para_Sz_mem_rand_read_write_2L2 - Para_Sz_mem_rand_read_write_1k) + Para_Sz_mem_rand_read_write_1k;
    }
    else if (array_size <= x3) {
        return (double)(array_size - x2) / (double)(x3 - x2) * (Para_Sz_mem_rand_read_write_05L3 - Para_Sz_mem_rand_read_write_2L2) + Para_Sz_mem_rand_read_write_2L2;
    }
    else if (array_size <= x4) {
        return (double)(array_size - x3) / (double)(x4 - x3) * (Para_Sz_mem_rand_read_write_15L3 - Para_Sz_mem_rand_read_write_05L3) + Para_Sz_mem_rand_read_write_05L3;
    }
    else if (array_size <= x5) {
        return (double)(array_size - x4) / (double)(x5 - x4) * (Para_Sz_mem_rand_read_write_3L3 - Para_Sz_mem_rand_read_write_15L3) + Para_Sz_mem_rand_read_write_15L3;
    }
    else if (array_size <= x6) {
        return (double)(array_size - x5) / (double)(x6 - x5) * (Para_Sz_mem_rand_read_write_10L3 - Para_Sz_mem_rand_read_write_3L3) + Para_Sz_mem_rand_read_write_3L3;
    }
    else {
        return Para_Sz_mem_rand_read_write_10L3;
    }
}