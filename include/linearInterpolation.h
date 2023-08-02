#ifndef PROJECT_NAME_constants_hpp
#define PROJECT_NAME_constants_hpp
namespace constants {
  constexpr double ENERGYPARTON[20] = {10.7,11.4,13.4,15.7,22.9,29.7,40.8,56.1,61.8,89.6,
                           118.0,174.4,246.9,349.9,496.3,704.8,1001.8,1425.6,2030.6,2895.5};
  constexpr double NG[20] = {4.59,4.68,4.90,5.13,5.44,6.02,6.47,6.92,7.28,7.67,
                 8.25,8.60,9.04,9.49,9.93,10.37,10.81,11.26,11.70,12.14};

  constexpr double XS[20] = {4.922e9,3.652e9,1.671e9,728.9e6,85.94e6,17.25e6,2.121e6,229.0e3,72.97e3,2.733e3,
                 235.4,6.720,0.284,0.012,5.112e-4,21.65e-6,0.9017e-6,36.45e-9,1.419e-9,52.07e-12};
  constexpr int LEN = 20;

}
#endif

double Lerp(double, double, double);
int searchInterval(double, const double [], int);
double getLerp(double, const double[], const double [], int);
void getIntegral(double [], const double [], const double [], int);
void getCDF(double [], const double [], const double [], int);
void getPDF(double [], const double [], const double [], int);
double getInterpoCDF(double,const double [], const double [], const double [], int);
double invertCDF(double, const double [], const double [], const double [], int);
