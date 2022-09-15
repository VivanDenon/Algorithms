#include <iostream>
#include <vector>
#include <ctime>
#include <functional>
#include <random>
#include <math.h>
#include <fstream>

#define MAX_N  2000
#define MAX_ITER 100

using namespace std;

template <typename T>
T linear(T x) {
    return x;
}

template <typename T>
void linear_test(const vector<T> &v) {
    auto i = v.cbegin();
    while(i != v.cend()) {
        linear(*i);
        ++i;
    }
}

template <typename T>
double sum(const vector<T> &v) {
    double s = 0;
    auto i = v.cbegin();
    while(i != v.cend()) {
        s += *i;
        ++i;
    }
    return s;
}

template <typename T>
double product(const vector<T> &v) {
    double p = 1;
    auto i = v.cbegin();
    while(i != v.cend()) {
        p *= *i;
        ++i;
    }
    return p;
}

template <typename T>
double polynom(const vector<T> &v, double x) {
    double s = 0;
    for(size_t k = 0; k < v.size(); ++k) {
        s += v[k] * pow(x, k);
    }
    return s;
}

template <typename T>
double polynom_horner(const vector<T> &v, double x) {
    double r = 0;
    auto i = v.crbegin();
    while(i != v.crend()) {
        r *= x;
        r += *i;
        ++i;
    }
    return r;
}

template <typename T>
void bubble_sort(vector<T> &v) {
    for(size_t n = 0; n < v.size(); ++n) {
        for(size_t m = 0; m < v.size() - 1 - n; ++m) {
            if(v[m] > v[m + 1])
            {
                swap(v[m], v[m + 1]);
            }
        }
    }
}


template <typename T>
size_t partition (vector<T> &v, size_t left, size_t right) 
{ 
    T pivot = v[right]; 
    size_t i = left; 
    for (size_t j = left; j < right; j++) 
    { 
        if (v[j] <= pivot) 
        { 
            swap(v[i], v[j]);
            ++i;
        } 
    } 
    swap(v[i], v[right]); 
    return i; 
} 

template <typename T>
void quick_sort(vector<T> &v, size_t left, size_t right) {
    size_t pivot_i = partition<T>(v, left, right);
    if(pivot_i > left) 
        quick_sort<T>(v, left, pivot_i - 1);
    if(pivot_i < right) 
        quick_sort<T>(v, pivot_i + 1, right); 
}

template <typename T>
void quick_sort_test(vector<T> &v)
{
    quick_sort(v, 0, v.size());
}

const int block_size = 64;

template <typename T>
void insertion_sort(vector<T> &v, int left, int right)
{
	for (int i = left + 1; i <= right; ++i)
	{
		int temp = v[i];
		int j = i - 1;
		while (j >= left && v[j] > temp)
		{
			v[j + 1] = v[j];
			--j;
		}
		v[j + 1] = temp;
	}
}

template <typename T>
void merge_sort(vector<T> &v, int l, int m, int r)
{
	int len1 = m - l + 1, len2 = r - m;
    vector<T> left(len1), right(len2);
	for (int i = 0; i < len1; ++i)
		left[i] = v[l + i];
	for (int i = 0; i < len2; ++i)
		right[i] = v[m + 1 + i];

	int i = 0;
	int j = 0;
	int k = l;

	while (i < len1 && j < len2)
	{
		if (left[i] <= right[j])
		{
			v[k] = left[i];
			++i;
		}
		else
		{
			v[k] = right[j];
			++j;
		}
		++k;
	}

	while (i < len1)
	{
		v[k] = left[i];
		++k;
		++i;
	}

	while (j < len2)
	{
		v[k] = right[j];
		++k;
		++j;
	}
}

template <typename T>
void tim_sort(vector<T> &v)
{
	for (int i = 0; i < v.size(); i += block_size)
		insertion_sort<T>(v, i, min((i + block_size - 1), int(v.size() - 1)));

	for (int size = block_size; size < v.size(); size = 2 * size)
	{
		for (int left = 0; left < v.size(); left += 2 * size)
		{
			int mid = left + size - 1;
			int right = min((left + 2 * size - 1), int(v.size() - 1));
			if(mid < right)
				merge_sort<T>(v, left, mid, right);
		}
	}
}

template <typename T>
void product_matrix(const vector<vector<T>> &a, const vector<vector<T>> &b, vector<vector<double>> &r) {
    int n = a[0].size();
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            for(int k = 0; k < n; ++k) {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
} 


clock_t get_average_runtime(const function<void()> &f, const function<void()> &f_rand) {
    auto average_time = 0;
    for(int i = 0; i < MAX_ITER; ++i) {
        f_rand();
        auto t = clock();
        f();
        average_time += clock() - t;
    }
    return average_time / MAX_ITER;
} 

template <typename T>
void set_random(vector<T> &v, T min, T max){
    random_device rd;
    default_random_engine eng(rd());
    uniform_int_distribution<T> distr(min, max);

    auto i = v.begin();
    while(i != v.end()) {
        *i = distr(eng);
        ++i;
    }
}

template <typename T>
void set_random_matrix(vector<vector<T>> &m, T min, T max) {
    auto i = m.begin();
    while(i != m.end()) {
        set_random<T>(*i, min, max);
        ++i;
    }
} 


template <typename T>
void calculate_run_time(const string& file_name) {
    ofstream out(file_name);
    if(!out.is_open()) 
        return;

    out << "N,Linear,Sum,Product,Polynom,Polynom Horner,Bubble sort,Quick sort,Tim sort,Matrix product\n"; 

    vector<T> v;
    v.reserve(MAX_N);
    
    vector<vector<double>> r(MAX_N, vector<double>(MAX_N));
    vector<vector<T>> a(MAX_N, vector<T>()), b(MAX_N, vector<T>());
    for(int i = 0; i < MAX_N; ++i) {
        a[i].reserve(MAX_N);
        b[i].reserve(MAX_N);
    }
    
    vector<function<void()>> funcs = {bind(&linear_test<T>, ref(v)), bind(&sum<T>, ref(v)), bind(&product<T>, ref(v)), 
            bind(&polynom<T>, ref(v), 1.5), bind(&polynom_horner<T>, ref(v), 1.5), bind(&bubble_sort<T>, ref(v)), 
            bind(&quick_sort_test<T>, ref(v)), bind(&tim_sort<T>, ref(v))}; 

    auto f_set_random =  bind(&set_random<T>, ref(v), 1, numeric_limits<T>::max());
    auto f_set_random_m = bind([&] (vector<vector<T>> &a, vector<vector<T>> &b) 
                                { 
                                    set_random_matrix<T>(a, 1, numeric_limits<T>::max()); 
                                    set_random_matrix<T>(b, 1, numeric_limits<T>::max()); 
                                }, 
                                ref(a), ref(b));
    auto f_m = bind(&product_matrix<T>, ref(a), ref(b), ref(r));

    for(int n = 1; n <= MAX_N; ++n) 
    {
        v.push_back(0);
        out << n << ','; 
        for(auto f : funcs) {
            out << get_average_runtime(f, f_set_random) << ','; 
        }

        if(n <= 150) 
        {
            for(int i = 0; i < n; ++i) 
            {
                a[i].resize(n);
                b[i].resize(n);
            }
            out << get_average_runtime(f_m, f_set_random_m);
        }
        out << endl;
    }

    out.close();
}



int main() {
    calculate_run_time<u_int16_t>("uint16.csv");
    calculate_run_time<u_int32_t>("uint32.csv");
    calculate_run_time<u_int64_t>("uint64.csv");
    return 0;
}