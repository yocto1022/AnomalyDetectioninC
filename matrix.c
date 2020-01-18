#include "matrix.h"

void Cholesky_decomposition(int size1,int size2,double A[size1][size2],double out[size1][size2]){
	int i,j,k;
	double temp;

		//行列を初期化(Lが0行列じゃないとまずいので)
		for(i=0;i<size1;i++){
			for(j=0;j<size2;j++){
				out[i][j] = 0.0;
			}
		}

	for(i=0;i<size1;i++){
		for(j=0;(j<i+1) && (j < size1);j++){
			if(i==j){
				if(j==0){
					out[i][j] = sqrt(A[i][j]);
				}else{
					temp = 0.0;
					for(k=0;k<=j-1;k++){
						temp += out[j][k] * out[j][k];
					}
					out[i][j] = sqrt(A[i][j]-temp);
				}
			}else{
				temp=0.0;
				for(k=0;k<=j-1;k++){
					temp += out[i][k] * out[j][k];
				}
				out[i][j] = (A[i][j] - temp)/out[j][j];
			}
		}
	}

	return;
}


void inverse(int size1,int size2,double A[size1][size2],double out[size1][size2]){
	int i,j,k,x_pv;
	double temp[size1][size2],temp2[size1],temp3[size1];
	double buf,pv_tmp;

	//まず掃き出し法のための単位行列をつくる&&tempに入力行列をコピー
	for(i=0;i<size1;i++){
		for(j=0;j<size2;j++){
			temp[i][j] = A[i][j];
			if(i==j){
				out[i][j] = 1.0;
			}else{
				out[i][j] = 0.0;
			}
		}
	}

	//掃き出しを開始
	for(i=0;i<size1;i++){

		//ピボット開始
		//対角最大成分の選択
		pv_tmp = temp[i][i];
		x_pv = 0;
		for(j=i;j<size1;j++){
			if(fabs(pv_tmp) < fabs(temp[j][i])){
				x_pv = j;
				pv_tmp = temp[i][j];
			}
		}

		//ピボット対象があったら
		if(x_pv != 0){
			if(fabs(pv_tmp) < 0.000001){
				fprintf(stderr,"逆行列:成分が小さすぎる(逆行列生成不可)");
				exit(1);
			}
		//ピボットの行入れ替え
			for(k=0;k<size1;k++){
				//tempのほう
				//退避
				temp2[k]=temp[i][k];
				//退避させたところに入れ替え対象を入れる
				temp[i][k] = temp[x_pv][k];
				//退避させたものを入れ替え先に入れる
				temp[x_pv][k] = temp2[k];
				
				//outのほう(同上)
				temp3[k]=out[i][k];
				out[i][k] = out[x_pv][k];
				out[x_pv][k] = temp3[k];
				
			}
		}

		buf = 1.0 /temp[i][i];
		for(j=0;j<size1;j++){
			temp[i][j] *= buf;
			out[i][j] *= buf;
		}

		//消してく
		for(j=0;j<size1;j++){
			if(i!=j){
				buf = temp[j][i];
				for(k=0;k<size1;k++){
					temp[j][k] -= temp[i][k] * buf;
					out[j][k] -= out[i][k] * buf; 
				}
			}
		}
	}
	return;
}

