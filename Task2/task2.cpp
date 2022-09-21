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

const double EPS  = 0.01;
const double Fi = (sqrt(5.0) + 1.0) / 2.0;

double cube(double x) {
    return pow(x, 3.0);
}

double func1(double x) {
    return abs(x - 0.2);
}

double func2(double x) {
    return x * sin(1 / x);
}

double exhaustive_search(double left, double right, double eps, function<double(double)> f) {
    int n = (right - left) / eps;
    double min = f(left);
    for(int i = 1; i <= n; i++) {
        auto r =  f(left);
        if(r < min) 
            min = r;
        left += eps;
    }
    return min; 
}


int count_exhaustive_search(double left, double right, double eps, function<double(double)> f) {
    return (right - left) / eps;
}

double dichotomy (double a, double b, double eps, function<double(double)> f) {
    double mid = (b + a) / 2;
    if(b - a <= eps) 
        return f(mid);

    double x1 = (a + mid) / 2, 
        x2 = (b + mid) / 2;

    if(f(x1) < f(x2)) 
        return dichotomy(a, mid, eps, f);
    else
        return dichotomy(mid, b, eps, f); 
}

int count_dichotomy (double a, double b, double eps, function<double(double)> f) {
    double mid = (b + a) / 2;
    if(b - a <= eps) 
        return 1;

    double x1 = (a + mid) / 2, 
        x2 = (b + mid) / 2;

    if(f(x1) < f(x2)) 
        return 1 + count_dichotomy(a, mid, eps, f);
    else
        return 1 + count_dichotomy(mid, b, eps, f); 
}

double golden_section(double a, double b, double eps, function<double(double)> f) {
    double mid = (b + a) / 2;
    if(b - a <= eps) 
        return f(mid);
    
    double x1 = b - (b - a) / Fi,
        x2 = a + (b - a) / Fi;

    if(f(x1) < f(x2)) 
        return golden_section(a, x2, eps, f);
    else
        return golden_section(x1, b, eps, f); 
}

int count_golden_section(double a, double b, double eps, function<double(double)> f) {
    double mid = (b + a) / 2;
    if(b - a <= eps) 
        return 1;
    
    double x1 = b - (b - a) / Fi,
        x2 = a + (b - a) / Fi;

    if(f(x1) < f(x2)) 
        return (1 + count_golden_section(a, x2, eps, f));
    else
        return (1 + count_golden_section(x1, b, eps, f)); 
}

double linear_appr(double x, const vector<double> &c) {
    return c[0] * x + c[1];
}

double real_appr(double x, const vector<double> &c) {
    return c[0] / (1 + c[1] * x); 
}

void print(double x) {
    cout << setw(20) << x;
}

void print(const string &s) {
    cout << setw(20) << s;
}

// Sum of squares diviations
double SSD(const vector<pair<double, double>> &data, const vector<double> &c, const function<double(double, vector<double>)> &f) {
    double sum = 0;
    for(auto &p : data) {
        sum += pow( f(p.first, c) - p.second, 2.0 );  
    }
    return sum;
}

// Standart diviation
double SD(const vector<pair<double, double>> &data, const vector<double> &c, const function<double(double, vector<double>)> &f) {
    return sqrt(SSD(data, c, f) / data.size());
}

// Expected value
double ED(const vector<pair<double, double>> &data, const vector<double> &c, const function<double(double, vector<double>)> &f) {
    double sum = 0;
    for(auto &p : data) {
        sum += f(p.first, c);  
    }
    return sum / data.size();
}

double exhaustive_search_nD(vector<double>::iterator iter_c, vector<double>::iterator end_c,
    vector<pair<double, double>>::const_iterator iter_b, double eps, const function<double()> &f) {
    
    // Initialization
    *iter_c = iter_b->first;
    int n = (iter_b->second - iter_b->first) / eps;
    double min = f();
    vector<double> cmin(end_c - iter_c);

    // Find min(f) recursively for each vector of parameters where parameter in range (left board, right board) 
    for(int i = 0; i <= n; i++) {
        auto r = (iter_c < end_c - 1) ? exhaustive_search_nD(iter_c + 1, end_c, iter_b + 1, eps, f) : f();
        if(r < min) {
            for(int i = 0; i < end_c - iter_c; ++i) {
                cmin[i] = *(iter_c + i);
            }
            min = r;
        }
        *iter_c += eps;
    }

    for(int i = 0; i < end_c - iter_c; ++i) {
        *(iter_c + i) = cmin[i];
    }
    return min;
}

double coordinate_descent(vector<double>::iterator iter_c, vector<double>::iterator end_c,
    vector<pair<double, double>>::const_iterator iter_b, double eps, const function<double()> &f) {

    // Initialization
    for(int i = 0; i < end_c - iter_c; ++i) {
       *(iter_c + i) = (iter_b + i)->first;
    }

    // Find min(f) with exhaustive_search for each c while f(Xk + 1) - f(Xk) > EPS 
    double min = f(), last_min;
    int i = 0;
    do{
        last_min = min;
        min = exhaustive_search_nD(iter_c + i, iter_c + 1 + i, iter_b + i, eps, f);
        i = (i + 1) % (end_c - iter_c);
    }
    while(abs(min - last_min) > EPS);
    return min; 
}

double nelder_mead(vector<double>::iterator iter_x, vector<double>::iterator end_x,
    vector<pair<double, double>>::const_iterator iter_b, double eps, const function<double()> &f) {
    
    // Reflection coefficient, effect coefficient, impact coefficient
    double alpha = 1, beta = 0.5, gamma = 2;
    int n = end_x - iter_x;
    
    // Vector pair of value of Fi(Xi) and i
    vector<pair<double, size_t>> Fi(n + 1);
    // Vector all points Xi
    vector<vector<double>> Xi(n + 1, vector<double>(n));
    vector<double> Xh(n), Xg(n), Xl(n), Xc(n), Xr(n), Xs(n);

    // Initialization
    auto ix = iter_x;
    auto ib = iter_b;
    for(int i = 0; i < n; ++i) {
       *ix = ib->first;
       ++ix; ++ib;
    }

    // Preparation
    for(int i = 0; i <= n; ++i) {
        ix = iter_x;
        ib = iter_b;
        auto iter_Xc = Xc.begin();
        auto iter_Xi = Xi[i].begin();
        while(ix != end_x) {
            // Uniform distribution of points on the interval (b.first, b.second) 
            *ix += (ib->second - ib->first) / (n + 2);
            *iter_Xi = *ix;
            *iter_Xc += *ix;
            ++ix; ++ib, ++iter_Xi;
        }
        Fi[i].first = f();
        Fi[i].second = i;
    }

    sort(Fi.begin(), Fi.end());  
    for(auto f : Fi) {
        cout << f.first << ' ' << f.second << endl;
    }

    size_t h = Fi[n].second, g = Fi[n - 1].second, l = Fi[0].second;
    Xh = Xi[h]; Xg = Xi[g], Xl = Xi[l];

    //Find the gravity centre for all point except xh
    auto iter_Xc = Xc.begin();
    auto iter_Xh = Xh.begin();
    while(iter_Xc != Xc.end()) {
        *iter_Xc = (*iter_Xc - *iter_Xh) / n;
        cout << *iter_Xc << ' ';
        ++iter_Xc, ++iter_Xh;
    }
    cout << endl;



}

int main() {
    vector< function<double (double, double, double, function<double(double)>) > > funcs = { exhaustive_search, dichotomy, golden_section, 
        count_exhaustive_search, count_dichotomy, count_golden_section };

    vector<string> f_names = { "exhaustive_search", "Dichotomy", "Golden Section", "Count i Num", "Count i Dich", "Count i Gold"};
    print("Type\\Func");
    print("Min(x^3)");
    print("Min(|x - 0.2|)");
    print("Min(x * sin(1/x))");   
    cout << endl;

    int i = 0;
    for(auto f : funcs) {
        print(f_names[i]);
        print(f(0.0, 1.0, EPS, &cube)); 
        print(f(0.0, 1.0, EPS, &func1)); 
        print(f(0.01, 1, EPS, &func2));
        cout << endl;
        ++i;
    }

    ofstream out("appr.csv");
    if(!out.is_open()) 
        return -1;


    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(0, 1);
    vector<double> ab = { distr(eng), distr(eng)};
    cout << ab[0] << ' ' << ab[1] << endl;
    vector<pair<double, double>> data(100);

    int k = 0;
    for(auto &p : data) {
        p.first = double(k) / 100;
        p.second = linear_appr(p.first, ab);
        ++k;
    }

    normal_distribution<double> distr_normal(ED(data, ab, linear_appr), SD(data, ab, linear_appr));
    for(auto &p : data) {
        p.second = linear_appr(p.first, ab) + distr_normal(eng);
    }

    vector<double> c(2);
    vector<pair<double, double>> boards = {make_pair(-1.0, 1.0), make_pair(-1.0, 1.0)};

    double r = exhaustive_search_nD(c.begin(), c.end(), boards.cbegin(), EPS, bind(SSD, cref(data), cref(c), linear_appr));
    cout << r << ' ' << c[0] << ' ' << c[1] << endl;

    //r = exhaustive_search_nD(c.begin(), c.end(), boards.cbegin(), EPS, bind(SSD, cref(data), cref(c), real_appr));
    cout << r << ' ' << c[0] << ' ' << c[1] << endl;

    r = coordinate_descent(c.begin(), c.end(), boards.cbegin(), EPS, bind(SSD, cref(data), cref(c), linear_appr));
    cout << r << ' ' << c[0] << ' ' << c[1] << endl;

    r = coordinate_descent(c.begin(), c.end(), boards.cbegin(), EPS, bind(SSD, cref(data), cref(c), real_appr));
    cout << r << ' ' << c[0] << ' ' << c[1] << endl;

    r = nelder_mead(c.begin(), c.end(), boards.cbegin(), EPS, bind(SSD, cref(data), cref(c), real_appr));
    cout << r << ' ' << c[0] << ' ' << c[1] << endl;

   /* auto ab_gauss = gauss_appr(data, -1, 1, EPS, linear_appr, dlinear_appr);
    cout <<  ab_gauss.first << ' ' << ab_gauss.second << endl;
    
    auto ab = exhaustive_search_appr(data, -1, 1, EPS, linear_appr);
    cout <<  ab.first << ' ' << ab.second << endl;
    auto ab_real = exhaustive_search_appr(data, -1, 1, EPS, real_appr);
    cout <<  ab.first << ' ' << ab.second;

    out << "x" << ',' << "f(x) = ax + b" << ',' << "y with noise" << ',' << "linear appr" 
            << ',' << "real appr" << endl;
    for(auto &p : data) {
        out << p.first << ',' << (a * p.first + b) << ',' << p.second << ',' << linear_appr(p.first, ab.first, ab.second) 
            << ',' << real_appr(p.first, ab_real.first, ab_real.second) << endl;
    }*/
    return 0;
}
