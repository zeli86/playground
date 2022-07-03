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

#include "gsl/gsl_sf.h"
#include "math.h"

typedef double (*TEF)(const double, const double);

inline double AiryEF_1(const double a, const double x)
{
  return 1.4261046287334950e+00 * sqrt(a) * gsl_sf_airy_Ai(a * x - 2.3381074104597670e+00, GSL_PREC_DOUBLE);
}

inline double AiryEF_2(const double a, const double x)
{
  return 1.2451573191271700e+00 * sqrt(a) * gsl_sf_airy_Ai(a * x - 4.0879494441309710e+00, GSL_PREC_DOUBLE);
}

inline double AiryEF_3(const double a, const double x)
{
  return 1.1557967485952720e+00 * sqrt(a) * gsl_sf_airy_Ai(a * x - 5.5205598280955510e+00, GSL_PREC_DOUBLE);
}

inline  double AiryEF_4(const double a, const double x)
{
  return 1.0978747223053990e+00 * sqrt(a) * gsl_sf_airy_Ai(a * x - 6.7867080900717590e+00, GSL_PREC_DOUBLE);
}

inline  double AiryEF_5(const double a, const double x)
{
  return 1.0555920040103590e+00 * sqrt(a) * gsl_sf_airy_Ai(a * x - 7.9441335871208530e+00, GSL_PREC_DOUBLE);
}

inline double AiryEF_6(const double a, const double x)
{
  return 1.0225755972118040e+00 * sqrt(a) * gsl_sf_airy_Ai(a * x - 9.0226508533409800e+00, GSL_PREC_DOUBLE);
}

inline double AiryEF_7(const double a, const double x)
{
  return 9.9564889221441920e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.0040174341558090e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_8(const double a, const double x)
{
  return 9.7300997897372820e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.1008524303733260e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_9(const double a, const double x)
{
  return 9.5354277742716440e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.1936015563236260e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_10(const double a, const double x)
{
  return 9.3651034928134670e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.2828776752865760e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_11(const double a, const double x)
{
  return 9.2140181626123890e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.3691489035210720e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_12(const double a, const double x)
{
  return 9.0784916103259980e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.4527829951775340e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_13(const double a, const double x)
{
  return 8.9557892079791080e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.5340755135978000e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_14(const double a, const double x)
{
  return 8.8438261710575260e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.6132685156945770e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_15(const double a, const double x)
{
  return 8.7409785312090940e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.6905633997429940e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_16(const double a, const double x)
{
  return 8.6459578462332320e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.7661300105697060e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_17(const double a, const double x)
{
  return 8.5577255660263640e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.8401132599207120e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_18(const double a, const double x)
{
  return 8.4754329372771240e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.9126380474246950e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_19(const double a, const double x)
{
  return 8.3983778518618640e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 1.9838129891721500e+01, GSL_PREC_DOUBLE);
}

inline double AiryEF_20(const double a, const double x)
{
  return 8.3259732318344970e-01 * sqrt(a) * gsl_sf_airy_Ai(a * x - 2.0537332907677570e+01, GSL_PREC_DOUBLE);
}

static TEF EF_AIRY[] = {&AiryEF_1, &AiryEF_2, &AiryEF_3, &AiryEF_4, &AiryEF_5, &AiryEF_6, &AiryEF_7, &AiryEF_8, &AiryEF_9, &AiryEF_10, &AiryEF_11, &AiryEF_12, &AiryEF_13, &AiryEF_14, &AiryEF_15, &AiryEF_16, &AiryEF_17, &AiryEF_18, &AiryEF_19, &AiryEF_20};

