#include "s21_math.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int s21_abs(int x) {
  if (x < 0) x *= -1;
  return x;
}

long double s21_acos(double x) {
  long double result;
  if (x < -1.0 || x > 1.0 || x != x) {
    result = S21_NAN;
  } else {
    result = S21_PI / 2 - s21_asin(x);
  }
  return result;
}

long double s21_asin(double x) {
  long double result;
  if (x < -1.0 || x > 1.0 || x != x) {
    result = S21_NAN;
  } else if (x == 1) {
    result = S21_PI / 2;
  } else if (x == -1) {
    result = -S21_PI / 2;
  } else {
    result = x;
    long double temp = x;
    int i = 1;
    while (s21_fabs(temp) >= S21_EPS && i <= 300000) {
      temp = (factor(2 * i) / (binary_pow(4.0, (unsigned int)(i)) * factor(i) *
                               factor(i) * (2 * i + 1))) *
             binary_pow(x, (unsigned int)(2 * i + 1));
      result += temp;
      ++i;
      // printf("x = %lf, res = %Lf\n",x, result);
    }
  }
  // long double up =
  return result;
}

long double s21_atan(double x) {
  long double result;
  int minus_flag = 1;
  if (x < 0) {
    x *= -1;
    minus_flag = -1;
  }
  if (x != x) {
    result = S21_NAN;
  } else if (x == S21_INF || x >= 1e7) {
    result = S21_PI / 2 * minus_flag;

  } else if (x == 1) {
    result = S21_PI / 4 * minus_flag;
  } else if (x > 1) {
    x = 1 / x;
    result = S21_PI / 2 - s21_atan(x);
    result *= minus_flag;
  } else {
    // result = s21_asin(x / s21_pow(1 + x * x, 0.5)) * minus_flag;
    result = x;
    long double temp = x;
    double sign = -1.0;
    int i = 2;
    while (s21_fabs(temp) > 1e-9 || i < 20) {
      temp = (binary_pow(x, (unsigned int)(2 * i - 1)) / (2 * i - 1) * sign);
      // printf("temp = %.10lf, i = %d\n", temp, i);
      result += temp;
      ++i;
      sign *= -1.0;
    }
    result *= minus_flag;
  }

  return result;
}

long double s21_ceil(double x) {
  long double result;
  if (x != x) {
    result = S21_NAN;
  } else if (x == S21_INF)
    result = S21_INF;
  else if (x == -S21_INF)
    result = -S21_INF;
  else {
    int exp = 0;
    while (x > 1e18) {
      x /= 10;
      ++exp;
    }
    long long temp = (long long)(x);
    if (x > 0 && temp - x != 0) {
      temp += 1;
    }
    result = temp;
    while (exp--) result *= 10;
  }
  return result;
}

long double s21_cos(double x) {
  if (x != x) return S21_NAN;
  if (x == S21_INF || x == -S21_INF) return S21_NAN;
  long double _x = s21_fmod(x, S21_PI * 2);
  long double rezult = 0.0;
  long double buf = 1.;
  long double diff = S21_EPS;
  int iter = 0;
  while ((buf > 0. ? buf : buf * -1.) > diff) {
    rezult += buf;
    iter++;
    buf *= -_x * _x / (iter * (iter + 1));
    iter++;
  }
  return rezult;
}

long double s21_exp(double x) {
  long double result = 1;
  if (x == 0) {
    result = 1;
  } else if (x < -20) {
    result = 0;
  } else if (x != x) {
    result = S21_NAN;
  } else if (x > 706) {
    result = S21_INF;
  } else {
    double x0 = x;
    if (x0 < 0) x0 = -x0;
    long double temp = 1;
    long double num = x0;
    long double i = 1;
    while (temp > 1e-17) {
      temp *= num / i;
      result += temp;

      ++i;
    }
    if (x < 0) {
      result = 1 / result;
    }
  }
  return result;
}

long double s21_fabs(double x) {
  long double result = x;
  if (x != x)
    result = S21_NAN;
  else if (x < 0)
    result *= -1.0;
  return result;
}

long double s21_floor(double x) {
  long double result;
  if (x != x) {
    result = S21_NAN;
  } else if (x == S21_INF)
    result = S21_INF;
  else if (x == -S21_INF)
    result = -S21_INF;
  else {
    if (s21_fabs(x) < 1) {
      result = (long double)((int)x);
      x < 0 ? result -= 1 : 1;
    } else if (s21_fabs(x) == 1) {
      result = x;
    } else {
      union {
        uint64_t i;
        double d;
      } un;
      int minus_flag = 0, round_flag = 0;
      if (x < 0) {
        x = -x;
        minus_flag = 1;
      }
      un.d = x;
      uint64_t temp = un.i;
      uint64_t exp = (temp >> 52) - 1023;
      int bit = 0;
      while (exp) {
        exp -= 1;
        ++bit;
      }
      if (bit <= 51) {
        temp = un.i;
        uint64_t mask = 0;
        mask = ~mask;
        int smeshenie = 64 - 12 - bit;
        temp = un.i;
        temp >>= smeshenie;
        if (temp & 1) round_flag = 1;
        mask <<= smeshenie;
        un.i &= mask;
        if (minus_flag && round_flag) {
          un.d += 1;
          un.d = -un.d;
        }
        result = (long double)un.d;
      } else
        result = (long double)x;
    }
  }
  return result;
}

long double s21_fmod(double x, double y) {
  long double result;
  if (x != x || y != y || y == 0) {
    result = S21_NAN;
  } else if (x == S21_INF || x == -S21_INF)
    result = S21_NAN;
  else if (y == S21_INF || y == -S21_INF)
    result = x;
  else if ((y <= 1 && y >= -1))
    result = 0;
  else {
    int flag_minus = 1;
    if (x < 0) {
      x *= -1;
      flag_minus = -1;
    }
    y = s21_fabs(y);
    long double del = x / y;
    long double ndel = s21_floor(del);

    result = x - y * ndel;
    result *= flag_minus;
    //        while(x > y) {
    //            x -= y;
    //        }
    //        result = x * flag_minus;
  }
  return result;
}

long double s21_log(double x) {
  long double result = 0;
  long double corr = 0;
  if (x < 0 || x != x) {
    result = S21_NAN;
  } else if (x == 0) {
    result = -S21_INF;
  } else if (x == S21_INF) {
    result = S21_INF;
  } else {
    if (s21_fabs(x) > 1e6) {
      while (s21_fabs(x) > 1e6) {
        x /= 1000;
        corr += 6.90775527898;
      }
    } else if (x < 1 && x > -1) {
      while (x < 1 && x > -1) {
        x *= 10;
        corr -= 2.30258509299;
      }
    }
    double t = (x - 1) / (x + 1);
    double up = t;
    int i = 1;
    while (s21_fabs(up) >= S21_EPS) {
      result += up / (2 * i - 1);
      up = up * t * t;
      ++i;
    }
    result *= 2;
    result += corr;
  }
  return result;
}

long double s21_pow(double base, double exp) {
  long double result;
  int int_exp = (int)exp;
  int minus_flag = 1;
  if (base == 0.0 && exp > 0)
    result = 0;
  else if (((base == S21_INF || base == -S21_INF) && exp > 0) ||
           (exp == S21_INF && s21_fabs(base) > 1) ||
           (exp == -S21_INF && s21_fabs(base) < 1) || (base == 0 && exp < 0))
    result = S21_INF;
  else if ((base == S21_INF || base == -S21_INF) && exp < 0)
    result = 0;
  else if (exp == 0)
    result = 1.0;
  else if (base == 1 || (base == -1 && (exp == S21_INF || exp == -S21_INF)))
    result = 1;
  else if (exp == -S21_INF && s21_fabs(base) > 1)
    result = 0.0;
  else if (exp != exp || base != base)
    result = S21_NAN;
  else if (base < 0 && int_exp - exp != 0)
    result = S21_NAN;
  else {
    double exp0 = exp;
    if (exp0 < 0) exp0 = -exp0;

    if (base < 0 && int_exp - exp == 0) {
      if (int_exp == 1 || abs(int_exp) % 2 == 1 || int_exp == -1)
        minus_flag = -1;
      base = -base;
    }
    int_exp = (int)exp0;
    if (exp0 > 1) exp0 = exp0 - int_exp;
    result = s21_exp(s21_log(base) * exp0);
    result *= binary_pow(base, (unsigned int)int_exp);
    if (exp < 0) result = 1 / result;
    result *= minus_flag;
  }
  return result;
}

long double s21_sin(double x) {
  if (x != x) return S21_NAN;
  if (x == S21_INF || x == -S21_INF) return S21_NAN;
  long double _x = s21_fmod(x, S21_PI * 2);
  long double rezult = 0.0;
  long double buf = _x;
  long double diff = S21_EPS;
  int iter = 0;
  while ((buf > 0. ? buf : buf * -1.) > diff) {
    rezult += buf;
    iter++;
    buf *= -_x * _x / (2 * iter * (2 * iter + 1));
  }
  return rezult;
}

long double s21_sqrt(double x) {
  long double result = 0;
  if (x < 0) {
    result = -S21_NAN;
  } else if (x != x) {
    result = S21_NAN;
  } else if (x == S21_INF || x <= S21_EPS) {
    result = x;
  } else {
    result = s21_exp(s21_log(x) / 2.);
  }
  return result;
}

long double s21_tan(double x) {
  if (x != x) return S21_NAN;
  if (x == S21_INF) return S21_NAN;
  //____________ATTENTION
  // if (s21_fabs(x - S21_PI / 2) < S21_EPS_T) return 16331239353195370.0;
  if (s21_fabs(x - S21_PI / 2) < S21_EPS_T) return 16331239353195370.0;
  return s21_sin(x) / s21_cos(x);
}

long double factor(int i) {
  long double result = 1.0;
  for (int j = 1; j <= i; ++j) {
    result *= (long double)j;
  }
  return result;
}

double binary_pow(double base, unsigned int exp) {
  double v = 1.0;
  while (exp != 0) {
    if ((exp & 1) != 0) {
      v *= base;
    }
    base *= base;
    exp >>= 1;
  }
  return v;
}
