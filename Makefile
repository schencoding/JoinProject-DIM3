all:
	g++ --std=c++11 -mavx -O3 server_parameter_estimator.cpp -o server_parameter_estimator
	g++ --std=c++11 -mavx -O3 main.cpp -o main
	g++ --std=c++11 -mavx -O3 -fopenmp main_parallel.cpp -o main_parallel
