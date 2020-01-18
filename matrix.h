#ifndef MAT_H_
#define MAT_H_

//インクルード関連
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

//関数群
void Cholesky_decomposition(int size1,int size2,double A[size1][size2],double out[size1][size2]);	//コレスキー分解
void inverse(int size1,int size2,double A[size1][size2],double out[size1][size2]);					//逆行列z
#endif