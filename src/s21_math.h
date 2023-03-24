#ifndef s21_math_h
#define s21_math_h

#include <stdio.h>

#define MAX_INT 2147483647
#define MIN_INT -2147483646
#define S21_EXP 2.7182818284590452353602874713526624977572
#define S21_NAN 0.0 / 0.0
#define S21_INF 1.0 / 0.0
#define S21_EPS_T 1e-6
#define S21_EPS 1e-16
#define S21_PI 3.14159265358979323846264338327950288

int s21_abs(int x);
long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);
long double s21_ceil(double x);
long double s21_cos(double x);
long double s21_exp(double x);
long double s21_fabs(double x);
long double s21_floor(double x);
long double s21_fmod(double x, double y);
long double s21_log(double x);
long double s21_pow(double base, double exp);
long double s21_sin(double x);
long double s21_sqrt(double x);
long double s21_tan(double x);

long double factor(int i);
double binary_pow(double base, unsigned int exp);

#endif /* s21_math_h */
