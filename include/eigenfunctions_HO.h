/* * atus-pro testing - atus-pro testing playgroung
 * Copyright (C) 2017 Želimir Marojević <zelimir.marojevic@gmail.com>
 *
 * This file is part of atus-pro testing.
 *
 * atus-pro testing is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * atus-pro testing is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with atus-pro testing.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

inline double HO_0(const double a, const double x)
{
  double retval = pow(a, 0.1e1 / 0.4e1) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1);
  return retval;
}

inline double HO_1(const double a, const double x)
{
  double retval = sqrt(0.2e1) * pow(a, 0.3e1 / 0.4e1) * x * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1);
  return retval;
}

inline double HO_2(const double a, const double x)
{
  double retval = sqrt(0.2e1) * pow(a, 0.1e1 / 0.4e1) * (0.2e1 * a * x * x - 0.1e1) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.2e1;
  return retval;
}

inline double HO_3(const double a, const double x)
{
  double retval = sqrt(0.3e1) * pow(a, 0.3e1 / 0.4e1) * x * (0.2e1 * a * x * x - 0.3e1) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.3e1;
  return retval;
}

inline double HO_4(const double a, const double x)
{
  double retval = sqrt(0.6e1) * pow(a, 0.1e1 / 0.4e1) * (0.4e1 * a * a * pow(x, 0.4e1) - 0.12e2 * a * x * x + 0.3e1) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.12e2;
  return retval;
}

inline double HO_5(const double a, const double x)
{
  double retval = sqrt(0.15e2) * pow(a, 0.3e1 / 0.4e1) * x * (0.4e1 * a * a * pow(x, 0.4e1) - 0.20e2 * a * x * x + 0.15e2) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.30e2;
  return retval;
}

inline double HO_6(const double a, const double x)
{
  double retval = sqrt(0.5e1) * pow(a, 0.1e1 / 0.4e1) * (0.8e1 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.60e2 * a * a * pow(x, 0.4e1) + 0.90e2 * a * x * x - 0.15e2) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.60e2;
  return retval;
}

inline double HO_7(const double a, const double x)
{
  double retval = sqrt(0.70e2) * pow(a, 0.3e1 / 0.4e1) * x * (0.8e1 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.84e2 * a * a * pow(x, 0.4e1) + 0.210e3 * a * x * x - 0.105e3) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.420e3;
  return retval;
}

inline double HO_8(const double a, const double x)
{
  double retval = sqrt(0.70e2) * pow(a, 0.1e1 / 0.4e1) * (0.16e2 * pow(a, 0.4e1) * pow(x, 0.8e1) - 0.224e3 * pow(a, 0.3e1) * pow(x, 0.6e1) + 0.840e3 * a * a * pow(x, 0.4e1) - 0.840e3 * a * x * x + 0.105e3) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.1680e4;
  return retval;
}

inline double HO_9(const double a, const double x)
{
  double retval = sqrt(0.35e2) * pow(a, 0.3e1 / 0.4e1) * x * (0.16e2 * pow(a, 0.4e1) * pow(x, 0.8e1) - 0.288e3 * pow(a, 0.3e1) * pow(x, 0.6e1) + 0.1512e4 * a * a * pow(x, 0.4e1) - 0.2520e4 * a * x * x + 0.945e3) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.2520e4;
  return retval;
}

inline double HO_10(const double a, const double x)
{
  double retval = sqrt(0.7e1) * pow(a, 0.1e1 / 0.4e1) * (0.32e2 * pow(a, 0.5e1) * pow(x, 0.10e2) - 0.720e3 * pow(a, 0.4e1) * pow(x, 0.8e1) + 0.5040e4 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.12600e5 * a * a * pow(x, 0.4e1) + 0.9450e4 * a * x * x - 0.945e3) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.5040e4;
  return retval;
}

inline double HO_11(const double a, const double x)
{
  double retval = sqrt(0.154e3) * pow(a, 0.3e1 / 0.4e1) * x * (0.32e2 * pow(a, 0.5e1) * pow(x, 0.10e2) - 0.880e3 * pow(a, 0.4e1) * pow(x, 0.8e1) + 0.7920e4 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.27720e5 * a * a * pow(x, 0.4e1) + 0.34650e5 * a * x * x - 0.10395e5) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.55440e5;
  return retval;
}

inline double HO_12(const double a, const double x)
{
  double retval = sqrt(0.231e3) * pow(a, 0.1e1 / 0.4e1) * (0.64e2 * pow(a, 0.6e1) * pow(x, 0.12e2) - 0.2112e4 * pow(a, 0.5e1) * pow(x, 0.10e2) + 0.23760e5 * pow(a, 0.4e1) * pow(x, 0.8e1) - 0.110880e6 * pow(a, 0.3e1) * pow(x, 0.6e1) + 0.207900e6 * a * a * pow(x, 0.4e1) - 0.124740e6 * a * x * x + 0.10395e5) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.332640e6;
  return retval;
}

inline double HO_13(const double a, const double x)
{
  double retval = sqrt(0.6006e4) * pow(a, 0.3e1 / 0.4e1) * x * (0.64e2 * pow(a, 0.6e1) * pow(x, 0.12e2) - 0.2496e4 * pow(a, 0.5e1) * pow(x, 0.10e2) + 0.34320e5 * pow(a, 0.4e1) * pow(x, 0.8e1) - 0.205920e6 * pow(a, 0.3e1) * pow(x, 0.6e1) + 0.540540e6 * a * a * pow(x, 0.4e1) - 0.540540e6 * a * x * x + 0.135135e6) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.4324320e7;
  return retval;
}

inline double HO_14(const double a, const double x)
{
  double retval = sqrt(0.858e3) * pow(a, 0.1e1 / 0.4e1) * (0.128e3 * pow(a, 0.7e1) * pow(x, 0.14e2) - 0.5824e4 * pow(a, 0.6e1) * pow(x, 0.12e2) + 0.96096e5 * pow(a, 0.5e1) * pow(x, 0.10e2) - 0.720720e6 * pow(a, 0.4e1) * pow(x, 0.8e1) + 0.2522520e7 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.3783780e7 * a * a * pow(x, 0.4e1) + 0.1891890e7 * a * x * x - 0.135135e6) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.8648640e7;
  return retval;
}

inline double HO_15(const double a, const double x)
{
  double retval = sqrt(0.715e3) * pow(a, 0.3e1 / 0.4e1) * x * (0.128e3 * pow(a, 0.7e1) * pow(x, 0.14e2) - 0.6720e4 * pow(a, 0.6e1) * pow(x, 0.12e2) + 0.131040e6 * pow(a, 0.5e1) * pow(x, 0.10e2) - 0.1201200e7 * pow(a, 0.4e1) * pow(x, 0.8e1) + 0.5405400e7 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.11351340e8 * a * a * pow(x, 0.4e1) + 0.9459450e7 * a * x * x - 0.2027025e7) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.21621600e8;
  return retval;
}

inline double HO_16(const double a, const double x)
{
  double retval = sqrt(0.1430e4) * pow(a, 0.1e1 / 0.4e1) * (0.256e3 * pow(a, 0.8e1) * pow(x, 0.16e2) - 0.15360e5 * pow(a, 0.7e1) * pow(x, 0.14e2) + 0.349440e6 * pow(a, 0.6e1) * pow(x, 0.12e2) - 0.3843840e7 * pow(a, 0.5e1) * pow(x, 0.10e2) + 0.21621600e8 * pow(a, 0.4e1) * pow(x, 0.8e1) - 0.60540480e8 * pow(a, 0.3e1) * pow(x, 0.6e1) + 0.75675600e8 * a * a * pow(x, 0.4e1) - 0.32432400e8 * a * x * x + 0.2027025e7) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.172972800e9;
  return retval;
}

inline double HO_17(const double a, const double x)
{
  double retval = sqrt(0.12155e5) * pow(a, 0.3e1 / 0.4e1) * x * (0.256e3 * pow(a, 0.8e1) * pow(x, 0.16e2) - 0.17408e5 * pow(a, 0.7e1) * pow(x, 0.14e2) + 0.456960e6 * pow(a, 0.6e1) * pow(x, 0.12e2) - 0.5940480e7 * pow(a, 0.5e1) * pow(x, 0.10e2) + 0.40840800e8 * pow(a, 0.4e1) * pow(x, 0.8e1) - 0.147026880e9 * pow(a, 0.3e1) * pow(x, 0.6e1) + 0.257297040e9 * a * a * pow(x, 0.4e1) - 0.183783600e9 * a * x * x + 0.34459425e8) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.1470268800e10;
  return retval;
}

inline double HO_18(const double a, const double x)
{
  double retval = sqrt(0.12155e5) * pow(a, 0.1e1 / 0.4e1) * (0.512e3 * pow(a, 0.9e1) * pow(x, 0.18e2) - 0.39168e5 * pow(a, 0.8e1) * pow(x, 0.16e2) + 0.1175040e7 * pow(a, 0.7e1) * pow(x, 0.14e2) - 0.17821440e8 * pow(a, 0.6e1) * pow(x, 0.12e2) + 0.147026880e9 * pow(a, 0.5e1) * pow(x, 0.10e2) - 0.661620960e9 * pow(a, 0.4e1) * pow(x, 0.8e1) + 0.1543782240e10 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.1654052400e10 * a * a * pow(x, 0.4e1) + 0.620269650e9 * a * x * x - 0.34459425e8) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.8821612800e10;
  return retval;
}

inline double HO_19(const double a, const double x)
{
  double retval = sqrt(0.461890e6) * pow(a, 0.3e1 / 0.4e1) * x * (0.512e3 * pow(a, 0.9e1) * pow(x, 0.18e2) - 0.43776e5 * pow(a, 0.8e1) * pow(x, 0.16e2) + 0.1488384e7 * pow(a, 0.7e1) * pow(x, 0.14e2) - 0.26046720e8 * pow(a, 0.6e1) * pow(x, 0.12e2) + 0.253955520e9 * pow(a, 0.5e1) * pow(x, 0.10e2) - 0.1396755360e10 * pow(a, 0.4e1) * pow(x, 0.8e1) + 0.4190266080e10 * pow(a, 0.3e1) * pow(x, 0.6e1) - 0.6285399120e10 * a * a * pow(x, 0.4e1) + 0.3928374450e10 * a * x * x - 0.654729075e9) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.167610643200e12;
  return retval;
}

inline double HO_20(const double a, const double x)
{
  double retval = sqrt(0.46189e5) * pow(a, 0.1e1 / 0.4e1) * (0.1024e4 * pow(a, 0.10e2) * pow(x, 0.20e2) - 0.97280e5 * pow(a, 0.9e1) * pow(x, 0.18e2) + 0.3720960e7 * pow(a, 0.8e1) * pow(x, 0.16e2) - 0.74419200e8 * pow(a, 0.7e1) * pow(x, 0.14e2) + 0.846518400e9 * pow(a, 0.6e1) * pow(x, 0.12e2) - 0.5587021440e10 * pow(a, 0.5e1) * pow(x, 0.10e2) + 0.20951330400e11 * pow(a, 0.4e1) * pow(x, 0.8e1) - 0.41902660800e11 * pow(a, 0.3e1) * pow(x, 0.6e1) + 0.39283744500e11 * a * a * pow(x, 0.4e1) - 0.13094581500e11 * a * x * x + 0.654729075e9) * exp(-a * x * x / 0.2e1) * pow(0.3141592654e1, -0.1e1 / 0.4e1) / 0.335221286400e12;
  return retval;
}

static TEF EF_HO[] = {&HO_0, &HO_1, &HO_2, &HO_3, &HO_4, &HO_5, &HO_6, &HO_7, &HO_8, &HO_9, &HO_10, &HO_11, &HO_12, &HO_13, &HO_14, &HO_15, &HO_16, &HO_17, &HO_18, &HO_19, &HO_20};
