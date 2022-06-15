# Density-optimized Intersection-free Mapping and Matrix Multiplication for Join-Project Operations

This is the source code of the upcoming paper `'Density-optimized Intersection-free Mapping and Matrix Multiplication for Join-Project Operations'`, which has been accpeted in PVLDB22.

This repo contains implementations of the ideas and experiments discussed in the paper: 
- Join-Project Algorithm `DIM3`.
- `DIM3` with partial result caching. [Coming soon]
- `DIM3` for Join-Aggregate operation. [Coming soon]
- `DIM3` with DP for line join projection. [Coming soon]
- ...

**This repo is under construction. We are now working on cleaning up and commenting the codes...**

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
`g++` is all you need. Besides, your machine needs to support `AVX`.
```
make
```

## Running
```
./main
```
