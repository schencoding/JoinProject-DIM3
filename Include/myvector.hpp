// Copyright (c) 2021-2023 Institute of Computing Technology, Chinese Academy of Sciences
// DIM3 is licensed under Mulan PSL v2.

#pragma once

#include <bits/stdc++.h>
using namespace std;

template <class  T, typename _CountType = uint32_t>
struct myvector {
public:

	T* data;
	_CountType cur_num;
	_CountType max_num;

	myvector<T, _CountType>(_CountType ini_size = 8) {
		max_num = ini_size;
		data = (T*)malloc(sizeof(T) * max_num);
		cur_num = 0;
	}

	void init(_CountType ini_size = 8) {
		max_num = ini_size;
		data = (T*)malloc(sizeof(T) * max_num);
		cur_num = 0;
	}

	inline bool isEmpty() {
		return cur_num == 0;
	}

	void push_back(T x) {
		if (cur_num == max_num) {
			max_num <<= 1;
			auto tmp_p = (T*)realloc(data, sizeof(T) * max_num);
			if (tmp_p == NULL) {
				printf("Failed when resize myvector to %u\n", max_num);
				free(data);
				exit(1);
			}
			data = tmp_p;
		}
		data[cur_num] = x;
		cur_num++;
	}

	T pop_back() {
		if (cur_num == 0) return NULL;
		cur_num--;
		return data[cur_num];
	}

	void fixSize(_CountType size) {
		if (size == 0) return;
		auto tmp_p = (T*)realloc(data, sizeof(T) * size);
		if (tmp_p == NULL) {
			printf("Failed when fixsize myvector to %u\n", size);
			free(data);
			exit(1);
		}
		data = tmp_p;
		max_num = size;
		cur_num = size;
	}

	void reserve(_CountType size) {
		if (size < max_num) return;
		auto tmp_p = (T*)realloc(data, sizeof(T) * size);
		if (tmp_p == NULL) {
			printf("Failed when reserve myvector to %u\n", size);
			free(data);
			exit(1);
		}
		data = tmp_p;
		max_num = size;
	}

	inline _CountType size() const {
		return cur_num;
	}

	inline void clear() {
		cur_num = 0;
	}
};