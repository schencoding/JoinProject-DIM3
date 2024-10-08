all:
	g++ --std=c++11 -mavx -O3 server_parameter_estimator.cpp -o server_parameter_estimator
	g++ --std=c++11 -mavx -O3 main.cpp -o main
	g++ --std=c++11 -mavx -O3 -fopenmp main_parallel.cpp -o main_parallel
	g++ --std=c++11 -mavx -O3 join_aggregate.cpp -o join_aggregate

baseline_Classical: baseline_Classical.cpp
	g++ --std=c++11 -mavx -O3 baseline_Classical.cpp -o baseline_Classical

baseline_GEMM: baseline_GEMM.cpp
	g++ --std=c++11 -mavx -O3 baseline_GEMM.cpp -o baseline_GEMM -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -DMKL_ILP64  -m64  -I"${MKLROOT}/include"
