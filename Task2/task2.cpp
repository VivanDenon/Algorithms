#include <iostream>
#include <functional>
#include <math.h>
#include <vector>
#include <iomanip>
#include <random>

using namespace std;

const double EPS  = 0.001;
const double Fi = (sqrt(5) + 1) / 2;

double cube(double x) {
    return pow(x, 3.0);
}

double func1(double x) {
    return abs(x - 0.2);
}

double func2(double x) {
    return x * sin(1 / x);
}

double numeric(double a, double b, double eps, function<double(double)> f) {
    int n = (b - a) / eps;
    double min = f(a);
    for(int i = 1; i <= n; i++) {
        auto r =  f(a + eps * i);
        if(r < min) 
            min = r;
    }
    return min; 
}

int count_numeric(double a, double b, double eps, function<double(double)> f) {
    return (b - a) / eps;
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

double linear_appr(double x, double a, double b) {
    return a * x + b;
}

double real_appr(double x, double a, double b) {
    return a / (1 + b * x); 
}
 
double MSM(double a, double b, vector<pair<double, double>> data, function<double(double, double, double)> f) {
    double sum = 0;
    for(auto &p : data) {
        sum += pow( f(p.first, a, b) - p.second, 2.0 );  
    }
}

void print(double x) {
    cout << setw(20) << x;
}

void print(const string &s) {
    cout << setw(20) << s;
}



int main() {
    vector< function<double (double, double, double, function<double(double)>) > > funcs = { numeric, dichotomy, golden_section, 
        count_numeric, count_dichotomy, count_golden_section };

    vector<string> f_names = { "Numeric", "Dichotomy", "Golden Section", "Count i Num", "Count i Dich", "Count i Gold"};
    print("Type\\Func");
    print("Min(x^3)");
    print("Min(|x - 0.2|)");
    print("Min(x * sin(1/x))");   
    cout << endl;

    int i = 0;
    for(auto f : funcs) {
        print(f_names[i]);
        print(f(0, 1, EPS, &cube)); 
        print(f(0, 1, EPS, &func1)); 
        print(f(0.01, 1, EPS, &func2));
        cout << endl;
        ++i;
    }


    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(0, 1);
    double a = distr(eng), b = distr(eng);

    cout << a << ' ' << b << endl;

    normal_distribution<double> distr_normal(0, 1);
    vector<pair<double, double>> data(100);
    int k = 0;
    for(auto &p : data) {
        p.first = k / 100;
        p.second = a * p.first + b + distr_normal(eng);
    }
    
    return 0;
}
