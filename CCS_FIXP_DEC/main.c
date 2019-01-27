#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "types.h"    // Floating-point header file
#include "kuo_fixp_fft.h"

#define N   64              // Number of FFT points
#define EXP 6               // EXP=log2(N)
#define PI  3.1415926535897 //PI
#define M 	90 				//Paquetes a procesar
#define SC 8
/* Prototypes */
extern void my_double_fft(int16_t *in, int16_complex_t *out, double_complex_t *tmp);
extern void my_fixp_fft(int16_t *in, int16_complex_t *out, int16_complex_t *tmp);
int16_complex_t output_i16cplx[N>>1];		/* Output of the 16-point real-valued FFT */

//Variables
FILE *f1;
int16_complex_t coeff_temp[N];			/* Temp coeficientes a aplicar inv FFT */
int16_complex_t twiddle_factor[EXP];	/* Temp storage for twiddle factors */
int16_complex_t coeff_array[N*M];			//Array de coeficientes
int16_t read_array[N*M];			//Array de coeficientes de escritura en bin
int16_t output[N*M];


/*
 * main.c
 */
int main(void){
	//definiciÃ³n de variables internas
	unsigned short i,k;
	unsigned short size,steps;

	//Se inicializa los coeficientes Twiddle
	//a partir de que se evaluara 64 muestras
	kuo_fixp_fft_init(twiddle_factor, EXP);

	// Se obtiene el conjugado de los coeficientes
	for(i=0;i<EXP;i++){
		twiddle_factor[i].im = -twiddle_factor[i].im;
	}

	//Leer archivo bin con datos
	f1 = fopen("data2.bin", "rb");
		if(f1==NULL){
			printf("No existe dichero");
		}
	//Se obtiene total de coeff *2
	size=fread(read_array,1,sizeof(read_array),f1);
	fclose(f1);
	printf("Finalizo lectura archivo bin!\n");
	steps = size/N;
	printf("Size:%d y Steps: %d\n",size,steps);

	//Se crea un array donde se guarden los 64 coeff por cada
	//paquete a crear
	for(k=0;k<steps;k++){
		coeff_array[k*N].re = read_array[k*(N/2)] << SC;
		coeff_array[k*N].im = read_array[k*(N/2)+1] << SC;
		for(i=1;i<=(N/2);i++){
			coeff_array[i+k*N].re = read_array[2*i+k*N] << SC;
			coeff_array[i+k*N].im = read_array[2*i+1+k*N] << SC;
			coeff_array[64+k*N-i].re = read_array[2*i+k*N] << SC;
			coeff_array[64+k*N-i].im = -read_array[2*i+1+k*N] << SC;
			//printf("%d --- %d\n",i+k*N,64+k*N-i);
		}
	}

	for(k=0;k<steps;k++){
		//Asigna valores coeff temp de paquete
		for(i=0;i<N;i++){
			coeff_temp[i].re =coeff_array[i+k*N].re;
			coeff_temp[i].im =coeff_array[i+k*N].im;
			//printf("%d +j %d\n",coeff_temp[i].re,coeff_temp[i].im);
		}
		//Funcion de FFT inversa
		kuo_fixp_fft(coeff_temp, EXP, twiddle_factor, 0);

		//Asignacion a matriz de salida
		for(i=0;i<N;i++){
			output[i+k*N] = coeff_temp[i].re;
			//printf("%d\n",coeff_temp[i].re);
		}
	}
	printf("FIN DECODER\n");

	return 0;
}
