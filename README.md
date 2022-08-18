# Density-optimized Intersection-free Mapping and Matrix Multiplication for Join-Project Operations

These are the source codes of the PVLDB'22 paper `'Density-optimized Intersection-free Mapping and Matrix Multiplication for Join-Project Operations'`.

This repo contains implementations of the Join-Project Algorithm `DIM3` and its parallel version in the paper.

<!-- If you use this work, please cite our paper as follows
```
@inproceedings{zhou2021spitfire,
  title={Spitfire: A Three-Tier Buffer Manager for Volatile and Non-Volatile Memory},
  author={Zhou, Xinjing and Arulraj, Joy and Pavlo, Andrew and Cohen, David},
  booktitle={Proceedings of the 2021 International Conference on Management of Data},
  pages={2195--2207},
  year={2021}
}
``` -->

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