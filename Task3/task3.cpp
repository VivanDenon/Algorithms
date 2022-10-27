#include <iostream>
#include <functional>
#include <math.h>
#include <vector>
#include <iomanip>
#include <random>
#include <fstream>
#include <math.h>
#include <algorithm>

using namespace std;

const size_t N = 100;
const double EPS  = 0.001;
size_t COUNT_ITERATION = 0;
size_t COUNT_CALL_FUNCTION = 0;

void print_vector(vector<double>::iterator iter_x, vector<double>::iterator end_x) {
    while(iter_x != end_x) {
        cout << *iter_x << ' ';
        ++iter_x;
    }
    cout << endl;
}

void print_vector(const vector<double> &v) {
    for(auto e : v) 
        cout << e << ' ';
    cout << endl;
}

inline double linear(double x, const vector<double> &c) {
    return c[0] * x + c[1];
}

inline vector<double> grad_linear(double x, const vector<double> &c) {
    return {x, 1.0};    
}

inline double real(double x, const vector<double> &c) {
    return c[0] / (1.0 + c[1] * x); 
}

inline vector<double> grad_real(double x, const vector<double> &c) {
    return { 1.0 / (1.0 + c[1] * x), -(x * c[0]) / pow(1 + c[1] * x, 2.0) };      
}

inline void vector_sub(vector<double>::iterator iter_v1, vector<double>::iterator end_v1, 
    vector<double>::const_iterator iter_v2) 
{
    while(iter_v1 != end_v1) {
        *iter_v1 -= *iter_v2;
        ++iter_v1; ++iter_v2;
    } 
}

inline vector<double> vector_sub(vector<double>::const_iterator iter_v1, vector<double>::const_iterator end_v1, 
    vector<double>::const_iterator iter_v2) 
{
    vector<double> result(end_v1 - iter_v1);
    for(auto &r : result) {
        r = *iter_v1 - *iter_v2;
        ++iter_v1; ++iter_v2;
    }
    return result; 
}

inline void vector_add(vector<double>::iterator iter_v1, vector<double>::iterator end_v1, 
    vector<double>::const_iterator iter_v2) 
{
    while(iter_v1 != end_v1) {
        *iter_v1 += *iter_v2;
        ++iter_v1; ++iter_v2;
    } 
}

inline void vector_mul(vector<double>::iterator iter_v, vector<double>::iterator end_v, 
    double c) 
{
    while(iter_v != end_v) {
        *iter_v *= c;
        ++iter_v;
    } 
}

inline vector<double> vector_mul(vector<double>::const_iterator iter_v1, vector<double>::const_iterator end_v1, 
    vector<double>::const_iterator iter_v2) 
{
    vector<double> result(end_v1 - iter_v1);
    for(auto &r : result) {
        r = *iter_v1 * (*iter_v2);
        ++iter_v1; ++iter_v2;
    } 
    return result;
}

inline void vector_mul(vector<double>::iterator iter_v1, vector<double>::iterator end_v1, 
    vector<double>::const_iterator iter_v2) 
{
    while(iter_v1 != end_v1) {
        *iter_v1 *= *iter_v2;
        ++iter_v1;
    } 
}

inline vector<double> vector_mul(vector<double>::const_iterator iter_v, vector<double>::const_iterator end_v, 
    const vector<vector<double>> &m) 
{
    auto n = end_v - iter_v;
    vector<double> result(n);
    for(auto i = 0; i < n; ++i) {
        auto iv = iter_v;
        for(auto j = 0; j < n; ++j) {
            result[i] += *iv * m[j][i];
            ++iv;
        }
    }
    return result;
}

inline void vector_div(vector<double>::iterator iter_v1, vector<double>::iterator end_v1, 
    vector<double>::const_iterator iter_v2) 
{
    while(iter_v1 != end_v1) {
        *iter_v1 /= *iter_v2;
        ++iter_v1;
    } 
}

inline void vector_div(vector<double>::iterator iter_v1, vector<double>::iterator end_v1, 
    double c) 
{
    while(iter_v1 != end_v1) {
        *iter_v1 /= c;
        ++iter_v1;
    } 
}

inline void matrix_div(double c, vector<vector<double>> &m) 
{
    for(auto &v : m) 
        for(auto &e : v)
            e = 1 / e;
}

inline void matrix_add(vector<vector<double>> &m1, const vector<vector<double>> &m2) 
{
    for(auto i = 0; i < m1.size(); ++i){
        vector_add(m1[i].begin(), m1[i].end(), m2[i].cbegin());
    }
}

inline void matrix_mul(vector<vector<double>> &m, double c) 
{
    for(auto &v : m) 
        vector_mul(v.begin(), v.end(), c);
}

inline double euclid_distance(vector<double>::const_iterator iter_v1, vector<double>::const_iterator end_v1, 
    vector<double>::const_iterator iter_v2)
{
    double sum = 0;
    while(iter_v1 != end_v1) {
        sum += pow(*iter_v1 - *iter_v2, 2.0);
        ++iter_v1; ++iter_v2;
    }
    return sqrt(sum);
}

inline double vector_distance(vector<double>::const_iterator iter_v1, vector<double>::const_iterator end_v1) {
    double sum = 0;
    while(iter_v1 != end_v1) {
        sum += pow(*iter_v1, 2.0);
        ++iter_v1;
    } 
    return sqrt(sum);
}

// Sum of squares error
inline double SSE(const vector<pair<double, double>> &data, const vector<double> &c, const function<double(double, vector<double>)> &f) {
    ++COUNT_CALL_FUNCTION;
    double sum = 0;
    for(auto &p : data) {
        sum += pow( f(p.first, c) - p.second, 2.0 );  
    }
    return sum;
}

// Get vector gradient of SSE
vector<double> grad_SSE(const vector<pair<double, double>> &data, const vector<double> &c, 
    const function<double(double, vector<double>)> &f, 
    const function< vector<double>(double, vector<double>) > &grad_f) 
{
    ++COUNT_CALL_FUNCTION;
    
    vector<double> result(c.size());
    for(auto &p : data) {
        auto grad = grad_f(p.first, c);
        vector_mul(grad.begin(), grad.end(), (f(p.first, c) - p.second));
        vector_add(result.begin(), result.end(), grad.begin());
    }
    vector_mul(result.begin(), result.end(), 2);
    return result;
}

vector<vector<double>> hessian_matrix(vector<double>::iterator iter_x, vector<double>::iterator end_x, 
    const function<vector<double>()> &grad_f) 
{ 
    auto n = end_x - iter_x;
    vector<vector<double>> result(n, vector<double>(n));
    auto grad = grad_f();
    double dx = EPS;
    
    for(auto i = 0; i < n; ++i) {
        *iter_x += dx;
        auto grad_dx = grad_f();
        *iter_x -= dx;
        
        vector_sub(grad_dx.begin(), grad_dx.end(), grad.cbegin());
        for(auto j = 0; j < n; ++j) {
            result[i][j] = grad_dx[j] / dx;
        }
        ++iter_x;
    }    
    return result;
}


// Expected value
inline double ED(const vector<pair<double, double>> &data) {
    double sum = 0;
    for(auto &p : data) {
        sum += p.second;  
    }
    return sum / data.size();
}

// Standart diviation
inline double SD(const vector<pair<double, double>> &data) {
    double mid = ED(data), sum = 0;

    for(auto &p : data) {
        sum += pow( p.second - mid, 2.0 );  
    }
    return sqrt(sum / data.size());
}

double dichotomy (double& value, double left, double right, double eps, const function<double()> &f) {    
    double mid = (right + left) / 2;
    if(right - left <= eps) {
        value = mid;
        return f();
    }

    value = (left + mid) / 2;
    double f1 = f();
    value = (right + mid) / 2;
    double f2 = f();

    if(f1 < f2) 
        return dichotomy(value, left, mid, eps, f);
    else
        return dichotomy(value, mid, right, eps, f);
}

enum type_inizial {
    begin,
    half,
    end,
    zero
};

inline void set_start_point(vector<double>::iterator iter_x, vector<double>::iterator end_x, 
    vector<pair<double, double>>::const_iterator iter_b, 
    type_inizial type) 
{
    while(iter_x != end_x) {
        switch(type) {
            case type_inizial::begin:
                *iter_x = iter_b->first;
                break;
            case type_inizial::half:
                *iter_x = (iter_b->second + iter_b->first ) / 2;
                break;
            case type_inizial::end:
                *iter_x = iter_b->second;
                break;
            case type_inizial::zero:
                *iter_x = 0;
                break;
            default: 
                *iter_x = 0;
                break;
        }
        ++iter_x; ++iter_b;
    }
}

// Calculate new point in gradient descent
inline void set_point_by_step(vector<double>::iterator iter_x, vector<double>::iterator end_x, 
    vector<pair<double, double>>::const_iterator iter_b,
    vector<double>::const_iterator iter_g, 
    double step_size) 
{
    while(iter_x != end_x) {
        *iter_x -= *iter_g * step_size;

        if(*iter_x < iter_b->first)
            *iter_x = iter_b->first;
        else if(*iter_x > iter_b->second)
            *iter_x = iter_b->second;

        ++iter_x; ++iter_g; ++iter_b;
    }
}

inline double func_by_step(vector<double>::iterator iter_x, vector<double>::iterator end_x, 
    vector<pair<double, double>>::const_iterator iter_b,
    vector<double>::const_iterator iter_g, 
    const vector<double> &x,
    const double &step_size, 
    const function<double()> &f) 
{
    copy(x.cbegin(), x.cend(), iter_x);
    set_point_by_step(iter_x, end_x, iter_b, iter_g, step_size);
    return f();
}

// Calculate Barzilai-Borwein coefficient, where a1 = An, a2 = An-1
double Barzilai_Borwein(vector<double>::const_iterator iter_a1, vector<double>::const_iterator end_a1,
    vector<double>::const_iterator iter_a2,
    vector<double>::const_iterator iter_grad1, vector<double>::const_iterator end_grad1,
    vector<double>::const_iterator iter_grad2) 
{
        auto a = vector_sub(iter_a1, end_a1, iter_a2);
        auto g = vector_sub(iter_grad1, end_grad1, iter_grad2);
        vector_mul(a.begin(), a.end(), g.cbegin());
        return vector_distance(a.cbegin(), a.cend()) / pow(vector_distance(g.cbegin(), g.cend()), 2.0);
}

double gradient_descent(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, 
    const double eps, 
    const function<double()> &f, const function<vector<double>()> &grad_f) 
{
    // Initialization
    set_start_point(iter_x, end_x, iter_b, type_inizial::begin);
    vector<double> x = { iter_x, end_x }, grad = grad_f();
    set_start_point(iter_x, end_x, iter_b, type_inizial::half);

    do {
        auto grad_last = grad;
        grad = grad_f();
        double step_size = Barzilai_Borwein(iter_x, end_x, x.cbegin(), grad.cbegin(),  grad.cend(), grad_last.cbegin());

        copy(iter_x, end_x, x.begin());
        set_point_by_step(iter_x, end_x, iter_b, grad.cbegin(), step_size);

        ++COUNT_ITERATION;
    } while(vector_distance(grad.cbegin(), grad.cend()) >= eps);
    return f();
}

// Calculate Fletcher-Reeves coefficient, where a1 = An, a2 = An-1
double Fletcher_Reeves(vector<double>::const_iterator iter_a1, vector<double>::const_iterator end_a1,
    vector<double>::const_iterator iter_a2,vector<double>::const_iterator end_a2)
{
        auto num = vector_mul(iter_a1, end_a1, iter_a1);
        auto denom = vector_mul(iter_a2, end_a2, iter_a2);
        return accumulate(num.cbegin(), num.cend(), 0.0) / accumulate(denom.cbegin(), denom.cend(), 0.0);
}

double conjugate_gradient(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, 
    const double eps, 
    const function<double()> &f, const function<vector<double>()> &grad_f) 
{
    // Initialization
    set_start_point(iter_x, end_x, iter_b, type_inizial::begin);
    vector<double> grad = grad_f(), s = grad;
    vector_mul(s.begin(), s.end(), -1.0);
    double result = f();

    while(true) {
        double step_size_a;
        vector<double> x = {iter_x, end_x};
        result = dichotomy(step_size_a, -1.0, 0.0, eps * eps, bind(func_by_step, iter_x, end_x, iter_b, s.cbegin(), cref(x), cref(step_size_a), f)); 

        auto grad_last = grad;
        grad = grad_f();
        if(vector_distance(grad.cbegin(), grad.cend()) < eps)
            break;

        double step_size_b = Fletcher_Reeves(grad.cbegin(), grad.cend(), grad_last.cbegin(), grad_last.cend());
        auto a = grad;
        vector_mul(a.begin(), a.end(), -1.0);
        set_point_by_step(a.begin(), a.end(), iter_b, s.cbegin(), step_size_b * -1.0);
        s = a;

        ++COUNT_ITERATION;
    } 
    return result;
}

inline bool is_positive(const vector<vector<double>> &m) {
    for(auto i = 0; i < m.size(); ++i)
        if(m[i][i] < 0) 
            return false;
    return true;
}   

double Newtons_method(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, 
    const double eps, 
    const function<double()> &f, const function<vector<double>()> &grad_f) 
{
    // Initialization
    vector<double> x(end_x - iter_x);
    set_start_point(iter_x, end_x, iter_b, type_inizial::begin);
    vector<double> grad = grad_f();
    double result = f();

    while(euclid_distance(iter_x, end_x, x.cbegin()) >= eps) {
        auto H = hessian_matrix(iter_x, end_x, grad_f);
        if(is_positive(H)) {
            matrix_div(1.0, H);
            grad = vector_mul(grad.cbegin(), grad.cend(), H);
        }
        double step_size;
        copy(iter_x, end_x, x.begin());
        result = dichotomy(step_size, 0.0, 1.0, eps * eps, bind(func_by_step, iter_x, end_x, iter_b, grad.cbegin(), cref(x), cref(step_size), f));
        grad = grad_f();

        ++COUNT_ITERATION;
    } 
    return result;
}

inline vector<vector<double>> matrix_diagonal(const vector<vector<double>> &m) {
    vector<vector<double>> result (m.size(), vector<double> (m.size()));
    for(auto i = 0; i < m.size(); ++i)
        result[i][i] = m[i][i];
    return result;
}

inline double vector_max(const vector<double> &v) {
    double max = v[0];
    for(auto &e : v) 
        if(e > max) 
            max = e;
    return max;
}

double Levenberg_Marquardt_method(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, 
    const double eps, 
    const function<double()> &f, const function<vector<double>()> &grad_f) 
{
    vector<double> x(end_x - iter_x);
    set_start_point(iter_x, end_x, iter_b, type_inizial::begin);
    vector<double> grad = grad_f();
    double result = f(), lamda = 1, v = 2; //vector_max(grad), v = 2;

    while(euclid_distance(iter_x, end_x, x.cbegin()) >= eps) {
        auto H = hessian_matrix(iter_x, end_x, grad_f);
        if(is_positive(H)) {
            auto I = matrix_diagonal(H);
            matrix_mul(I, lamda);
            matrix_add(H, I);
            matrix_div(1.0, H);
            grad = vector_mul(grad.cbegin(), grad.cend(), H);
        }
        double step_size;
        copy(iter_x, end_x, x.begin());
        double new_result = dichotomy(step_size, 0.0, 1.0, eps * eps, bind(func_by_step, iter_x, end_x, iter_b, grad.cbegin(), cref(x), cref(step_size), f));
        auto new_grad = grad_f();

        auto g = (result - new_result) / euclid_distance(grad.begin(), grad.end(), new_grad.begin());
        if( g > 0) {
            result = new_result;
            grad = new_grad;
            lamda *= max<double>(1 / 3, 1 - pow(2 * g - 1, 3.0));
            v = 2;
        }  
        else {
            copy(x.begin(), x.end(), iter_x);
            lamda *= v;
            v *= 2;
        }
        ++COUNT_ITERATION;
    } 
    return result;
}

void print(double x) {
    cout << setw(15) << x;
}

void print(const string &s, size_t w = 15) {
    cout << setw(w) << s;
}

int main() {
    vector<vector<double> > temp = { { 1, 2}, {3, 4}};
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(0, 1);

    vector<double> ab = { distr(eng), distr(eng)};
    vector<pair<double, double>> data(N);
    vector<vector<double>> results(11, vector<double>(N));

    for(auto i = 0; i < data.size(); ++i) {
        data[i].first = double(i) / N;
        data[i].second = linear(data[i].first, ab);
        results[0][i] = data[i].first; // vector x
        results[1][i] = data[i].second; // vector y
    }

    normal_distribution<double> distr_normal(0, SD(data));
    for(auto i = 0; i < data.size(); ++i){ 
        data[i].second += distr_normal(eng);
        results[2][i] = data[i].second; // vector y with noise
    }

    // Vector of parameters which we searched. In this solution - A and B 
    vector<double> c(2);
    vector<pair<double, double>> bounds = {make_pair(-1.0, 1.0), make_pair(-1.0, 1.0)}; // Bounds for A and B

    // Vector of functions and their gradients
    vector<string> appr_funcs_names = { "linear", "real" };
    vector< function<double(double, vector<double>)> > appr_funcs = {linear, real};
    vector< function<vector<double>(double, vector<double>)> > grad_appr_funcs = { grad_linear, grad_real };
   
    vector<string> opt_funcs_names = { "Gradient descent", "Conjugate grad desc", "Newton's method", "Levenberg Marquardt" };
    vector< function<double(vector<double>::iterator, vector<double>::iterator, vector<pair<double, double>>::const_iterator, 
        double, function<double()>, function<vector<double>()>)> > opt_funcs_nD = {gradient_descent, conjugate_gradient, Newtons_method, Levenberg_Marquardt_method};
    
    ofstream out("appr.csv");
    if(!out.is_open()) 
        return -1;

    print("Method\\Result", 20);
    for(auto &f_a_n : appr_funcs_names) {
        print("SSE of " + f_a_n); 
        print("a"); 
        print("b");
        print("iter");
        print("call");
    }
    cout << endl;
    print("Original func", 20); print(SSE(data, ab, linear));  print(ab[0]);   print(ab[1]);   cout << endl; // SSE orininal func after added noise
    out << "x,y,y with noise";

    for(auto i = 0; i < opt_funcs_nD.size(); ++i) {
        print(opt_funcs_names[i], 20);
        for(auto j = 0; j < appr_funcs.size(); ++j) {
            print( opt_funcs_nD[i](c.begin(), c.end(), bounds.cbegin(), EPS, bind(SSE, cref(data), cref(c), appr_funcs[j]), 
                   bind(grad_SSE, cref(data), cref(c), appr_funcs[j], grad_appr_funcs[j])) ); 
            for(auto &c_i : c) 
                print(c_i);

            print(COUNT_ITERATION);
            COUNT_ITERATION = 0;
            print(COUNT_CALL_FUNCTION);
            COUNT_CALL_FUNCTION = 0;

            out << "," << opt_funcs_names[i] << " " << appr_funcs_names[j];
            for(auto k = 0; k < N; ++k) 
                results[2 * i + j + 3][k] = appr_funcs[j](data[k].first, c); // vector y after approximation
        }
        cout << endl;
    }
    out << endl;

    for(auto i = 0; i < N; ++i) {
        for(auto j = 0; j < 9; ++j) {
            if(j > 0)
                out << ',';
            out << results[j][i];
        }
        out << endl;
    }
    out.close();
    return 0;
}
