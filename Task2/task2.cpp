#include <iostream>
#include <functional>
#include <math.h>
#include <vector>
#include <iomanip>
#include <random>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "task3.cpp"

using namespace std;

const size_t N = 100;
const double EPS  = 0.001;
const double Fi = (sqrt(5.0) + 1.0) / 2.0;
size_t COUNT_ITERATION = 0;
size_t COUNT_CALL_FUNCTION = 0;

double cube(double x) {
    ++COUNT_CALL_FUNCTION;
    return pow(x, 3.0);
}

double func1(double x) {
    ++COUNT_CALL_FUNCTION;
    return abs(x - 0.2);
}

double func2(double x) {
    ++COUNT_CALL_FUNCTION;
    return x * sin(1 / x);
}

double exhaustive_search_1D(double left, double right, double eps, const function<double(double)> &f) {
    size_t n = (right - left) / eps;
    double min = f(left);
    for(auto i = 1; i <= n; i++) {
        auto r =  f(left);
        if(r < min) 
            min = r;
        left += eps;

        ++COUNT_ITERATION;
    }
    return min; 
}

double dichotomy (double left, double right, double eps, const function<double(double)> &f) {
    ++COUNT_ITERATION;
    
    double mid = (right + left) / 2;
    if(right - left <= eps) 
        return f(mid);

    double x1 = (left + mid) / 2, 
        x2 = (right + mid) / 2;

    if(f(x1) < f(x2)) 
        return dichotomy(left, mid, eps, f);
    else 
        return dichotomy(mid, right, eps, f);
}

double golden_section(double left, double right, double eps, const function<double(double)> &f) {
    ++COUNT_ITERATION;
    
    double mid = (right + left) / 2;
    if(right - left <= eps) 
        return f(mid);
    
    double x1 = right - (right - left) / Fi,
        x2 = left + (right - left) / Fi;

    if(f(x1) < f(x2)) 
        return golden_section(left, x2, eps, f);
    else 
        return golden_section(x1, right, eps, f); 
}

double linear(double x, const vector<double> &c) {
    return c[0] * x + c[1];
}

double real(double x, const vector<double> &c) {
    return c[0] / (1 + c[1] * x); 
}

// Sum of squares error
double SSE(const vector<pair<double, double>> &data, const vector<double> &c, const function<double(double, vector<double>)> &f) {
    ++COUNT_CALL_FUNCTION;
    double sum = 0;
    for(auto &p : data) {
        sum += pow( f(p.first, c) - p.second, 2.0 );  
    }
    return sum;
}

// Expected value
double ED(const vector<pair<double, double>> &data) {
    double sum = 0;
    for(auto &p : data) {
        sum += p.second;  
    }
    return sum / data.size();
}

// Standart diviation
double SD(const vector<pair<double, double>> &data) {
    double mid = ED(data), sum = 0;

    for(auto &p : data) {
        sum += pow( p.second - mid, 2.0 );  
    }
    return sqrt(sum / data.size());
}

double exhaustive_search_nD(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, double eps, const function<double()> &f) {
    
    // Initialization
    *iter_x = iter_b->first;
    size_t n = (iter_b->second - iter_b->first) / eps;
    double min = f();
    vector<double> cmin(end_x - iter_x);

    // Find min(f) recursively for each vector of parameters where parameter in range (left bound, right bound) 
    for(size_t i = 0; i <= n; i++) {
        auto r = (iter_x < end_x - 1) ? exhaustive_search_nD(iter_x + 1, end_x, iter_b + 1, eps, f) : f();
        if(r < min) {
            for(size_t i = 0; i < end_x - iter_x; ++i) {
                cmin[i] = *(iter_x + i);
            }
            min = r;
        }
        *iter_x += eps;
        ++COUNT_ITERATION;
    }

    for(size_t i = 0; i < end_x - iter_x; ++i) {
        *(iter_x + i) = cmin[i];
    }
    return min;
}

double coordinate_descent(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, double eps, const function<double()> &f) {

    // Initialization
    size_t n_x = end_x - iter_x;
    for(size_t i = 0; i < n_x; ++i) {
       *(iter_x + i) = (iter_b + i)->first;
    }

    // Find min(f) with exhaustive_search for each c while (Xk + 1  - Xk) > EPS 
    double min = f(), max_div;
    do{ 
        max_div = 0;
        for(size_t i = 0; i < n_x; ++i) {
            double t = *(iter_x + i);
            min = exhaustive_search_nD(iter_x + i, iter_x + 1 + i, iter_b + i, eps, f);
            if( fabs(*(iter_x + i) - t) > max_div) 
                max_div = fabs(*(iter_x + i) - t);
            ++COUNT_ITERATION;    
        }
    }
    while(max_div >= EPS);
    return min; 
}

void set_start_point(vector<double>::iterator iter_x, vector<double>::iterator end_x, 
    vector<pair<double, double>>::const_iterator iter_b) 
{
    while(iter_x != end_x) {
        *iter_x = ( iter_b->second + iter_b->first ) / 2;
        ++iter_x; ++iter_b;
    }
}

double euclid_distance(vector<double>::const_iterator iter_v1, vector<double>::const_iterator end_v1, 
    vector<double>::const_iterator iter_v2)
{
    double sum = 0;
    auto iv1 = iter_v1;
    auto iv2 = iter_v2;
    while(iv1 != end_v1) {
        sum += pow(*iv1 - *iv2, 2.0);
        ++iv1; ++iv2;
    }
    return sqrt(sum);
}

double nelder_mead(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, double eps, const function<double()> &f) {
    
    // Reflection coefficient, effect coefficient, impact coefficient
    double alpha = 1, beta = 0.5, gamma = 2;
    size_t n_x = end_x - iter_x;
    
    // Vector pair of value of Fk(Xk) and k
    vector<pair<double, size_t>> Fk(n_x + 1);
    // Vector all points Xk
    vector<vector<double>> Xk(n_x + 1, vector<double>(n_x));
    // Xh is point of max Fk, Xg is point of second max Fk, Xl is point of smallest Fk
    vector<double>  Xc(n_x), Xr(n_x), Xs(n_x), Xe(n_x);

    // Initialization
    set_start_point(iter_x, end_x, iter_b);

    // Preparation
    auto ix = iter_x;
    auto ib = iter_b;
    
    copy(iter_x, end_x, Xk[0].begin());
    Fk[0].first = f();
    Fk[0].second = 0;

    for(size_t k = 1; k <= n_x; ++k) {
        auto temp = *ix;
        *ix += ( ib->second - ib->first ) * eps;

        copy(iter_x, end_x, Xk[k].begin());
        Fk[k].first = f();
        Fk[k].second = k;

        *ix = temp;
        ++ix; ++ib;
    }   

    double sd, min;
    do { 
        // Sorting 
        sort(Fk.begin(), Fk.end());  
        size_t h = Fk[n_x].second, g = Fk[n_x - 1].second, l = Fk[0].second;
        double &Fh = Fk[n_x].first, Fg = Fk[n_x - 1].first, Fl = Fk[0].first;
        vector<double> &Xh = Xk[h], &Xl = Xk[l];  
        min = Fl;

        // Find the gravity centre for all point except Xh
        for(auto k = 0; k <= n_x; ++k) {
            if(k == h) 
                continue;
            for(int j = 0; j < n_x; ++j) {
                Xc[j] += Xk[k][j];
            }
        }

        for(auto &x : Xc) {
            x /= n_x;
        }

        // Lamda func for generate new point Xk with any coefficient and return value of Fk
        auto fg = [&] (vector<double> &X1, const vector<double> &X2, const vector<double> &X3, double c) -> double { 
            for(int j = 0; j < X1.size(); ++j) 
                X1[j] = (1 - c) * X2[j] + c * X3[j];
            copy(X1.begin(), X1.end(), iter_x);
            return f();
        };

        // Reflect the point Xh with respect to Xc with the reflection coefficient (alpha)
        double Fr = fg(Xr, Xc, Xh, (alpha * -1));

        // Expansion
        if(Fr < Fl) {
            double Fe = fg(Xe, Xc, Xr, gamma); 
            if(Fe < Fr) {
                Xh = Xe;
                Fh = Fe;
            }
            else {
                Xh = Xr;
                Fh = Fr;
                copy(Xr.begin(), Xr.end(), iter_x);
            }
            l = h;
            Xl = Xh;
            min = Fh;
        }
        else if(Fr < Fg) {
            Xh = Xr;
            Fh = Fr;
            copy(Xl.begin(), Xl.end(), iter_x);
        }
        else {
            if(Fr < Fh) {
                Xh = Xr;
                Fh = Fr;
            }

            // Shrinking
            double Fs = fg(Xs, Xc, Xh, beta);
            if(Fs < Fh) {
                Xh = Xs;
                Fh = Fs;
            }
            else {
                // Global shrinking
                auto iter_Fk = Fk.begin() + 1;
                for(auto k = 0; k <= n_x; ++k) {
                    if(k == l) 
                        continue;
                    for(auto j = 0; j < n_x; ++j)
                        Xk[k][j] = Xl[j] + (Xk[k][j] - Xl[j]) / 2; 

                    copy(Xk[k].begin(), Xk[k].end(), iter_x);
                    iter_Fk->first= f();
                    iter_Fk->second = k;
                    iter_Fk++;
                }
            }
            copy(Xl.begin(), Xl.end(), iter_x);
        }

        // Calculate euclid distance all point to min point
        sd = 0;
        for(auto i = 0; i <= n_x; ++i) {
            if(i == l)
                continue;
            sd += euclid_distance(Xk[i].cbegin(), Xk[i].cend(), Xl.begin());
        }
        sd = sd / n_x;
        ++COUNT_ITERATION;

    } while(sd >= eps);
    return min;
}


void print(double x) {
    cout << setw(15) << x;
}

void print(const string &s, size_t w = 15) {
    cout << setw(w) << s;
}

int main() {
    // Test of approximation  
    vector<string> opt_funcs_names = { "Exhaustive search", "Dichotomy", "Golden Section"};
    vector< function<double (double, double, double, function<double(double)>) > > opt_funcs_1D = { exhaustive_search_1D, dichotomy, golden_section };   
    vector<string> test_funcs_names = { "x^3", "|x - 0.2|", "x * sin(1/x)"};
    vector< function<double(double)> > test_funcs = { cube, func1, func2 };
    vector<pair<double, double>> bounds = {make_pair(0.0, 1.0), make_pair(0.0, 1.0), make_pair(0.01, 1.0)};


    print("Type\\Func", 20);    
    for(auto &name : test_funcs_names) {
        print(name); print("iter"); print("call"); 
    }
    cout << endl;

    for(auto i = 0; i < opt_funcs_1D.size(); ++i) {
        print(opt_funcs_names[i], 20); 
        for(auto j = 0; j < test_funcs.size(); ++j) {
            print(opt_funcs_1D[i](bounds[j].first, bounds[j].second, EPS, test_funcs[j])); 
            
            print(COUNT_ITERATION);
            COUNT_ITERATION = 0;
            print(COUNT_CALL_FUNCTION);
            COUNT_CALL_FUNCTION = 0;
        }
        cout << endl;
    }
    cout << endl;

    // Test of approximation n-dimension function. In this ssolution - SSE
    vector<string> appr_funcs_names = { "linear", "real" };
    vector< function<double(double, vector<double>)> > appr_funcs = { linear, real };
    vector< function<vector<double>(double, vector<double>)> > grad_appr_funcs = { task3::grad_linear, task3::grad_real };

    opt_funcs_names = { "Exhaustive search", "Coordinate descent", "Nelder Mead"};
    vector< function<double(vector<double>::iterator, vector<double>::iterator, vector<pair<double, double>>::const_iterator, 
        double, function<double()>)> > opt_funcs_nD = {exhaustive_search_nD, coordinate_descent, nelder_mead };

    vector<string> opt_newtons_funcs_names = { "Gradient descent", "Conjugate grad desc", "Newton's method", "Levenberg Marquardt" };
    vector< function<double(vector<double>::iterator, vector<double>::iterator, vector<pair<double, double>>::const_iterator, 
        double, function<double()>, function<vector<double>()>)> > opt_newtons_funcs_nD = {task3::gradient_descent, task3::conjugate_gradient, 
        task3::Newtons_method, task3::Levenberg_Marquardt_method};
    
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(0, 1);

    vector<double> ab = { distr(eng), distr(eng)};
    vector<pair<double, double>> data(N);
    auto count_columns = 3 + (opt_funcs_nD.size() + opt_newtons_funcs_names.size()) * appr_funcs.size();
    vector<vector<double>> results(count_columns, vector<double>(N));

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
    // Bounds for A and B
    bounds = {make_pair(-1.0, 1.0), make_pair(-1.0, 1.0)};
    
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
            print(opt_funcs_nD[i](c.begin(), c.end(), bounds.cbegin(), EPS, bind(SSE, cref(data), cref(c), appr_funcs[j]))); 
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

      for(auto i = 0; i < opt_newtons_funcs_nD.size(); ++i) {
        print(opt_newtons_funcs_names[i], 20);
        for(auto j = 0; j < appr_funcs.size(); ++j) {
            print( opt_newtons_funcs_nD[i](c.begin(), c.end(), bounds.cbegin(), EPS, bind(SSE, cref(data), cref(c), appr_funcs[j]), 
                   bind(task3::grad_SSE, cref(data), cref(c), appr_funcs[j], grad_appr_funcs[j])) ); 
            for(auto &c_i : c) 
                print(c_i);

            print(task3::COUNT_ITERATION);
            task3::COUNT_ITERATION = 0;
            print(task3::COUNT_CALL_FUNCTION);
            task3::COUNT_CALL_FUNCTION = 0;

            out << "," << opt_newtons_funcs_names[i] << " " << appr_funcs_names[j];
            for(auto k = 0; k < N; ++k) 
                results[2 * i + j + 3 + opt_funcs_nD.size() * appr_funcs.size()][k] = appr_funcs[j](data[k].first, c); // vector y after approximation
        }
        cout << endl;
    }
    out << endl;

    for(auto i = 0; i < N; ++i) {
        for(auto j = 0; j < count_columns; ++j) {
            if(j > 0)
                out << ',';
            out << results[j][i];
        }
        out << endl;
    }
    out.close();

    return 0;
}
