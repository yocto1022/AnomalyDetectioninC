#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mt19973ar.h"
#include "matrix.h"

#define NumofSample 50
#define NumofTest 50

//学習のハイパーパラメータ
#define lambda 0.1//正則化
#define h 0.1//RBFカーネルバンド幅

//最適化変数
double theta[NumofSample];//各RBFのウェイト

void optimize(double test[NumofTest],double sample[NumofSample]);//ウェイト最適化
double RBF(double in,double sample);//RBFカーネルの1要素を計算(出力はスカラー))
double calc_abnormality(double in, double sample[NumofSample]);//異常度の計算
double gauss_distribution(double devitation);//いつもの


int main(int argc, char const * argv[]){
	double sample[NumofSample];
	double test[NumofTest];
	double abnormality[NumofTest];
	int i;

	//サンプルとテスト(最後が駄目なやつ)を生成
	for(i=0;i<NumofSample;i++){
		sample[i] = gauss_distribution(1.0);
	}
	for(i=0;i<NumofTest;i++){
		test[i] = gauss_distribution(1.0);
	}
	test[NumofTest-1] = 5;//外れ値が混ざる

	optimize(test,sample);//最適化

	//ここから推定
	for(i=0;i<NumofTest;i++){
		abnormality[i] = calc_abnormality(test[i],sample);
		printf("%f\t%f\n",test[i],abnormality[i]);
	}

	return 0;
}

void optimize(double sample[NumofSample],double test[NumofTest]){
	double G_hat[NumofSample][NumofSample]={};
	double h_hat[NumofSample]={};
	double temp[NumofSample]={};
	int testIndex,sampleIndex1,sampleIndex2;
	double tempofCalcTheta[NumofSample][NumofSample];
	double out[NumofSample];

	//まずG_hatを計算する
	for(testIndex=0;testIndex<NumofTest;testIndex++){
		//カーネルの出力ベクトルを生成
		for(sampleIndex1=0;sampleIndex1<NumofTest;sampleIndex1++){
			temp[sampleIndex1] = RBF(test[testIndex],sample[sampleIndex1]);
		}
		//内積を取ってG_hatへ足す
		for(sampleIndex1=0;sampleIndex1<NumofSample;sampleIndex1++){
			for(sampleIndex2=0;sampleIndex2<NumofSample;sampleIndex2++){
				G_hat[sampleIndex1][sampleIndex2] += temp[sampleIndex1]*temp[sampleIndex2];
			}
		}
	}
	//全部終わったらサンプル数で平均取り
	for(sampleIndex1=0;sampleIndex1<NumofSample;sampleIndex1++){
		for(sampleIndex2=0;sampleIndex2<NumofSample;sampleIndex2++){
			G_hat[sampleIndex1][sampleIndex2] = G_hat[sampleIndex1][sampleIndex2]/NumofTest;
		}
	}

	//次にh_hat
	for(sampleIndex1=0;sampleIndex1<NumofSample;sampleIndex1++){
		for(sampleIndex2=0;sampleIndex2<NumofSample;sampleIndex2++){
			h_hat[sampleIndex2] += RBF(sample[sampleIndex1],sample[sampleIndex2]);
		}
	}
	//平均化
	for(sampleIndex1=0;sampleIndex1<NumofSample;sampleIndex1++){
		h_hat[sampleIndex2] = h_hat[sampleIndex2]/NumofSample;
	}

	//解析解の計算
	//まずかけられる項(正則化項の加算)
	for(sampleIndex1=0;sampleIndex1<NumofSample;sampleIndex1++){
		G_hat[sampleIndex1][sampleIndex1] += lambda;
	}

	inverse(NumofSample,NumofSample,G_hat,tempofCalcTheta);//逆行列導出

	//２つをかけておわり
	for(sampleIndex1=0;sampleIndex1<NumofSample;sampleIndex1++){
		for(sampleIndex2=0;sampleIndex2<NumofSample;sampleIndex2++){
			out[sampleIndex1] += G_hat[sampleIndex1][sampleIndex2] * h_hat[sampleIndex2];
		}
	}

	//代入して終わり
	for(sampleIndex1=0;sampleIndex1<NumofSample;sampleIndex1++){
		theta[sampleIndex1] = out[sampleIndex1];
	}


}

double calc_abnormality(double in, double sample[NumofSample]){
	int i;
	double out=0.0;
	for(i=0;i<NumofSample;i++){
		out += (theta[i]*RBF(in,sample[i]));
	}
	out = -log(out);

	return out;
}

double RBF(double in,double sample){
	return exp(-pow(in-sample,2.0)/(2*h*h));
}


double gauss_distribution(double devitation){
	return sqrt(-2.0 * pow(devitation,2) * log(genrand_real3())) * cos(2.0 * M_PI * genrand_real3());
}
