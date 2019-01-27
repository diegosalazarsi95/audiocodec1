#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>

typedef struct {
	double re;
	double im;
} double_complex_t;

typedef struct {
	double_complex_t coeff_fft[32];
} double_complex_t_array;

typedef struct {
	int16_t re;
	int16_t im;
} int16_complex_t;

typedef struct {
	int16_complex_t coeff_fft[32];
} int16_complex_t_array;

typedef struct {
	int32_t re;
	int32_t im;
} int32_complex_t;

#endif /* TYPES_H */
