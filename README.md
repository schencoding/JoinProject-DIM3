# Density-optimized Intersection-free Mapping and Matrix Multiplication for Join-Project Operations

These are the source codes of the PVLDB'22 paper `'Density-optimized Intersection-free Mapping and Matrix Multiplication for Join-Project Operations'`.

This repo contains implementations of the Join-Project Algorithm `DIM3` and its parallel version in the paper.

If you use this work, please cite our paper as follows
```
@article{DBLP:journals/pvldb/HuangC22,
  author    = {Zichun Huang and
               Shimin Chen},
  title     = {Density-optimized Intersection-free Mapping and Matrix Multiplication
               for Join-Project Operations},
  journal   = {Proc. {VLDB} Endow.},
  volume    = {15},
  number    = {10},
  pages     = {2244--2256},
  year      = {2022},
  url       = {https://www.vldb.org/pvldb/vol15/p2244-chen.pdf},
  timestamp = {Tue, 30 Aug 2022 09:02:56 +0200},
  biburl    = {https://dblp.org/rec/journals/pvldb/HuangC22.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

## Building
`g++` is all you need. Besides, your server needs to support `AVX`. For paralle version, you need `openmp`.
```
make
```

## Running
Single-threaded `DIM3` for Join-Project:
```
./main
```

Note: The parameters in `server_parameters.hpp` is measured on our own server. To get the parameters for your server, please first check your L1,L2,L3 cache sizes and fill them in lines `16-18` of `server_parameter_estimator.cpp`. Then,
```
make
./server_parameter_estimator 2>/dev/null
```
Please copy the printed parameters to lines `6-36` of `server_parameters.hpp`. BTW, if you kwow your workloads and can get more accurate parameters, that would be better!

---

For parallel version,
```
./main_parallel
```
The `THREAD_NUM` is defined in `dim3_parallel.hpp`.

For Join-Aggregate:
```
./join_aggregate
```

## Baselines
We have three baselines in our paper: "Classical", "Gemm(MKL)", and "DHK". The first two baselines can be build by
```
make baseline_Classical
make baseline_GEMM # Intel MKL needed
```
For "DHK", we obtained the code from the original author of that paper. Therefore, we do not have the right to make it open source in our repo.
