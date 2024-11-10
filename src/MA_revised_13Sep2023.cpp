#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include<string.h>
#include<sstream>
#include<iomanip>
#include <time.h>
#include <ctime>
#include<math.h>
using namespace std;

double StartTime, FinishTime, Runtime;
int Runs;		//每个算例运行次数

int* BestS;
int* CluLen;	//每个cluster的size，swap前后size不变
double Objective;
double ObjBest;

//N,K,D在main函数里给出
int N;		//the number of vertices
int K;		//the number of clusters
int D;		//dimensions
int Type;		//算例类型
double** Point;	//顶点坐标，坐标维度为D
double** MatrixNK;
double* ObjClu;	//每个cluster的目标值，没有除以该cluster的size
double** DisSquare;	//欧几里得空间的任意两点的欧式距离的平方

//for the memetic algorithm
typedef struct Population {
	int* p;
	double cost;
}Population;

Population* Pop;
Population Child;					//子代
Population child_update;
int** arr1;
int** arr2;
int* len;
int* len2;
int* match;
int* flagC1;
int* flagC2;
int* flagV;
int* unassV;
int* addressUnaV;

typedef struct Neighborhood {
	int  type;
	int  v;
	int  g;
	int  x;
	int  y;
}Neighborhood;

//for responsive threshold search
double TR;					//threshold ratio
double TA, TB, TC;			//the parameters
double ThreshT;				//threshold T
double ObjP;
Neighborhood* Neighbors;
int* randN;
int* randK;
int* flagN;
int* flagK;

//the following are used parameters
int PopNum;										//parameter, population size
int LL;											//exploration strength of the TBE procedure
int RTS_ITER;										//search depth of RTS

#define DEBUG
#define MAXNUM 99999999999999
#define NUM1	50
#define DEVIATION	0.000001

void initialing(string file)
{
	double dd;
	int row = 0, col = 0;
	string line;
	stringstream ss;
	fstream fs(file);
	if (!fs)
	{
		cout << "error open, File_Name " << file << endl;
		getchar();
		exit(1);
	}
	Point = new double* [N];
	for (int index = 0; index < N; index++)
	{
		Point[index] = new double[D];
	}
	while (getline(fs, line))
	{
		col = 0;
		ss << line;
		while (ss >> dd)			//算例文件经过格式化处理，每列数据以制表符分割
		{
			Point[row][col++] = dd;
		}
		row++;
#ifdef DEBUG1			
		for (int i = 0; i < D; i++)
		{
			cout << Point[row - 1][i] << ",";
		}
		cout << endl;
#endif
		ss.clear();
	}

	BestS = new int[N];
	DisSquare = new double* [N];
	for (int i = 0; i < N; i++)
	{
		DisSquare[i] = new double[N];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			DisSquare[i][j] = 0;
	}
	Pop = new Population[PopNum];
	for (int i = 0; i < PopNum; i++)
		Pop[i].p = new int[N];
	Child.p = new int[N];
	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			int dimens = 0;
			while (dimens < D)
			{
				DisSquare[i][j] += (Point[i][dimens] - Point[j][dimens]) * (Point[i][dimens] - Point[j][dimens]);
				dimens++;
			}
			DisSquare[j][i] = DisSquare[i][j];
		}
	}
}

void allocateMemory()
{
	MatrixNK = new double* [N];
	for (int i = 0; i < N; i++)
	{
		MatrixNK[i] = new double[K];
	}
	CluLen = new int[K];
	Neighbors = new Neighborhood[N * (N - 1) / 2 + N * K];
	ObjClu = new double[K];
	randN = new int[N];
	randK = new int[K];
	flagN = new int[N];
	flagK = new int[K];
	arr1 = new int* [K];
	arr2 = new int* [K];
	len = new int[K];
	len2 = new int[K];
	match = new int[K];
	flagC1 = new int[K];
	flagC2 = new int[K];
	flagV = new int[N];
	unassV = new int[N];
	addressUnaV = new int[N];
	for (int i = 0; i < K; i++)
	{
		arr1[i] = new int[N];
		arr2[i] = new int[N];
	}
	child_update.p = new int[N];
}

void randomConstruction(int* ss)
{
	int i, len = 0, ver, k;
	int* realLen = new int[K];
	for (i = 0; i < N; i++)
		ss[i] = -1;
	for (i = 0; i < K; i++)
		realLen[i] = 0;
	int nk = N / K;
	int nk1 = N % K;

	if (nk1 == 0)
		for (i = 0; i < K; i++)
			CluLen[i] = nk;
	else
	{
		for (i = 0; i < nk1; i++)
			CluLen[i] = nk + 1;
		for (i = nk1; i < K; i++)
			CluLen[i] = nk;
	}
	while (len < N)
	{
		ver = rand() % N;
		k = rand() % K;
		while (ss[ver] == -1 && realLen[k] < CluLen[k])
		{
			ss[ver] = k;
			realLen[k]++;
			len++;
		}

	}
	delete[] realLen;
}

//计算目标函数值
double caculateObj(int* ss)
{
	double dd;
	double obj = 0;
	for (int i = 0; i < K; i++)
	{
		dd = 0;
		for (int j = 0; j < N; j++)
		{
			for (int k = j + 1; k < N; k++)
			{
				if (ss[j] == i && ss[k] == i)
				{
					dd += DisSquare[j][k];
				}
			}
		}
		dd /= CluLen[i];
		obj += dd;
		//cout << "i=" << i << ",dd=" << dd << ",len=" << len << endl;
	}
	//cout << "in caculateObj method,obj=" << obj << endl;
	return obj;
}

//与之前的xiangjing.lai文章中的gamma矩阵含义相同
void initialMatrixNKAndObjClu(int* ss)
{
	int i, j;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < K; j++)
			MatrixNK[i][j] = 0;
	}
	for (i = 0; i < K; i++)
		ObjClu[i] = 0;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			MatrixNK[i][ss[j]] += DisSquare[i][j];
	for (i = 0; i < K; i++)
	{
		for (j = 0; j < N; j++)
		{
			if (ss[j] == i)
			{
				ObjClu[i] += MatrixNK[j][i];
			}
		}
		ObjClu[i] /= 2;
	}
}

void updateMatrixNK(int i, int g0, int g1)
{
	for (int j = 0; j < N; j++)
	{
		if (j != i)
		{
			MatrixNK[j][g0] -= DisSquare[i][j];
			MatrixNK[j][g1] += DisSquare[i][j];
		}
	}
}


void checkMove(double obj, int* ss)
{
	double obj1 = caculateObj(ss);
	if (fabs(obj - obj1) > DEVIATION)
	{
		printf("heer,obj=%6f,obj1=%6f,", obj, obj1);
		cout << "obj!=obj1##############,error,please check the delta value" << endl;
		getchar();
	}
}

double tbe(int l_iter, int* ss, double maxTime)
{
	double delta;
	int i, len1 = 0, len2 = 0, ele, ele2, clu;
	Objective = caculateObj(ss);
	//one-move move
	int l = 0;
	while ((clock() - StartTime) / CLOCKS_PER_SEC <= maxTime && l < l_iter)
	{
		len1 = 0; len2 = 0;
		for (i = 0; i < N; i++)
			flagN[i] = 0;
		for (i = 0; i < K; i++)
			flagK[i] = 0;
		while (len1 < N)
		{
			int rn = rand() % N;
			if (flagN[rn] == 0)
			{
				randN[len1++] = rn;
				flagN[rn] = 1;
			}
		}
		while (len2 < K)
		{
			int rk = rand() % K;
			if (flagK[rk] == 0)
			{
				randK[len2++] = rk;
				flagK[rk] = 1;
			}
		}
		if (N % K != 0)
		{
			for (int ind = 0; ind < N; ind++)
			{
				ele = randN[ind];
				for (int jnd = 0; jnd < K; jnd++)
				{
					clu = randK[jnd];
					//cout << "jnd="<<jnd<<",clu=" << clu << endl;
					if (CluLen[ss[ele]] == CluLen[clu] + 1)
					{
						delta = (ObjClu[ss[ele]] - MatrixNK[ele][ss[ele]]) / (CluLen[ss[ele]] - 1) + (MatrixNK[ele][clu] + ObjClu[clu]) / (CluLen[clu] + 1)
							- (ObjClu[ss[ele]] / CluLen[ss[ele]] + ObjClu[clu] / CluLen[clu]);
						//cout << "delta=" << delta << endl;
						if (Objective + delta < ThreshT)
						{
							int temp = ss[ele];
							ss[ele] = clu;
							CluLen[temp]--;
							CluLen[clu]++;
							ObjClu[temp] -= MatrixNK[ele][temp];
							ObjClu[clu] += MatrixNK[ele][clu];
							Objective += delta;
							updateMatrixNK(ele, temp, clu);
							if (Objective < ObjBest)
							{
								ObjBest = Objective;
								for (int i = 0; i < N; i++)
									BestS[i] = ss[i];
								FinishTime = clock();
							}
							//checkMove(Objective, ss);
						}
					}
				}
			}

		}
		//swap move
		for (int ind = 0; ind < N; ind++)
		{
			ele = randN[ind];
			for (int ind2 = ind + 1; ind2 < N; ind2++)
			{
				ele2 = randN[ind2];
				if (ss[ele] != ss[ele2])
				{
					delta = (MatrixNK[ele2][ss[ele]] - DisSquare[ele][ele2] - MatrixNK[ele][ss[ele]]) / CluLen[ss[ele]] +
						(MatrixNK[ele][ss[ele2]] - DisSquare[ele][ele2] - MatrixNK[ele2][ss[ele2]]) / CluLen[ss[ele2]];
					if (Objective + delta < ThreshT)
					{
						int temp = ss[ele];
						int temp2 = ss[ele2];
						ss[ele] = temp2;
						ss[ele2] = temp;
						ObjClu[temp] += MatrixNK[ele2][temp] - DisSquare[ele][ele2] - MatrixNK[ele][temp];
						ObjClu[temp2] += MatrixNK[ele][temp2] - DisSquare[ele][ele2] - MatrixNK[ele2][temp2];
						updateMatrixNK(ele, temp, temp2);
						updateMatrixNK(ele2, temp2, temp);
						Objective += delta;
						if (Objective < ObjBest)
						{
							ObjBest = Objective;
							for (int i = 0; i < N; i++)
								BestS[i] = ss[i];
							FinishTime = clock();
						}
						//checkMove(Objective, ss);
					}
				}
			}
		}
		//}
		l++;
	}
	return Objective;
}

//descent based improvement
double dbi(int* ss, double maxTime)
{
	double delta;
	int i, len1 = 0, len2 = 0, ele, ele2, clu;
	int flag_move = 1;
	Objective = caculateObj(ss);
	//one-move move
	if (N % K != 0)
	{
		while (flag_move && (clock() - StartTime) / CLOCKS_PER_SEC <= maxTime)
		{
			flag_move = 0;
			len1 = 0;
			len2 = 0;
			for (i = 0; i < N; i++)
				flagN[i] = 0;
			for (i = 0; i < K; i++)
				flagK[i] = 0;
			while (len1 < N)
			{
				int rn = rand() % N;
				if (flagN[rn] == 0)
				{
					randN[len1++] = rn;
					flagN[rn] = 1;
				}
			}
			while (len2 < K)
			{
				int rk = rand() % K;
				if (flagK[rk] == 0)
				{
					randK[len2++] = rk;
					flagK[rk] = 1;
				}
			}
			for (int ind = 0; ind < N; ind++)
			{
				ele = randN[ind];
				for (int jnd = 0; jnd < K; jnd++)
				{
					clu = randK[jnd];
					if (CluLen[ss[ele]] == CluLen[clu] + 1)
					{
						delta = (ObjClu[ss[ele]] - MatrixNK[ele][ss[ele]]) / (CluLen[ss[ele]] - 1) + (MatrixNK[ele][clu] + ObjClu[clu]) / (CluLen[clu] + 1)
							- (ObjClu[ss[ele]] / CluLen[ss[ele]] + ObjClu[clu] / CluLen[clu]);
						if (delta < 0)
						{
							int temp = ss[ele];
							ss[ele] = clu;
							CluLen[temp]--;
							CluLen[clu]++;
							ObjClu[temp] -= MatrixNK[ele][temp];
							ObjClu[clu] += MatrixNK[ele][clu];
							Objective += delta;
							updateMatrixNK(ele, temp, clu);
							//checkMove(Objective, ss);
							flag_move = 1;
						}
					}
				}
				//printf("ele=%d,Objective=%6f,threshT=%6f\n",ele,Objective,ThreshT);
			}
		}
	}
	//swap move
	flag_move = 1;
	while (flag_move && (clock() - StartTime) / CLOCKS_PER_SEC <= maxTime)
	{
		flag_move = 0;
		len1 = 0;
		for (i = 0; i < N; i++)
			flagN[i] = 0;
		while (len1 < N)
		{
			int rn = rand() % N;
			if (flagN[rn] == 0)
			{
				randN[len1++] = rn;
				flagN[rn] = 1;
			}
		}
		for (int ind = 0; ind < N; ind++)
		{
			ele = randN[ind];
			for (int ind2 = ind + 1; ind2 < N; ind2++)
			{
				ele2 = randN[ind2];
				if (ss[ele] != ss[ele2])
				{
					delta = (MatrixNK[ele2][ss[ele]] - DisSquare[ele][ele2] - MatrixNK[ele][ss[ele]]) / CluLen[ss[ele]] +
						(MatrixNK[ele][ss[ele2]] - DisSquare[ele][ele2] - MatrixNK[ele2][ss[ele2]]) / CluLen[ss[ele2]];
					if (delta < 0)
					{
						int temp = ss[ele];
						int temp2 = ss[ele2];
						ss[ele] = temp2;
						ss[ele2] = temp;
						ObjClu[temp] += MatrixNK[ele2][temp] - DisSquare[ele][ele2] - MatrixNK[ele][temp];
						ObjClu[temp2] += MatrixNK[ele][temp2] - DisSquare[ele][ele2] - MatrixNK[ele2][temp2];
						updateMatrixNK(ele, temp, temp2);
						updateMatrixNK(ele2, temp2, temp);
						Objective += delta;
						//checkMove(Objective, ss);
						flag_move = 1;
					}
				}
			}
		}
	}
	if (Objective < ObjBest)
	{
		ObjBest = Objective;
		for (int i = 0; i < N; i++)
			BestS[i] = ss[i];
		FinishTime = clock();
	}
	return Objective;
}

//辅助randomshake函数：构造邻域
void buildNeighbors()
{
	int i, j, g;
	int count;
	count = 0;
	for (i = 0; i < N; i++)
	{
		for (g = 0; g < K; g++)
		{
			Neighbors[count].type = 1;
			Neighbors[count].v = i;
			Neighbors[count].g = g;
			count++;
		}
	}
	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			Neighbors[count].type = 2;
			Neighbors[count].x = i;
			Neighbors[count].y = j;
			count++;
		}
	}
}

//随机扰动，使用one-move和swap move邻域
void randomShake(int L, int p1[])
{
	int v, g, x, y;
	int NumberNeighbors, old_g1, old_g2;
	int cur_index, count = 0;
	NumberNeighbors = N * (N - 1) / 2 + N * K;
	do
	{
		cur_index = rand() % NumberNeighbors;
		//cout << "cur_index=" << cur_index <<",type="<< Neighbors[cur_index].type<<  endl;
		if (Neighbors[cur_index].type == 1 && N % K != 0)
		{
			v = Neighbors[cur_index].v;
			g = Neighbors[cur_index].g;
			if (p1[v] != g && CluLen[p1[v]] == CluLen[g] + 1)
			{
				old_g1 = p1[v];
				p1[v] = g;
				CluLen[g] ++;
				CluLen[old_g1]--;
				count++;
			}
		}
		else if (Neighbors[cur_index].type == 2)
		{
			x = Neighbors[cur_index].x;
			y = Neighbors[cur_index].y;
			if (p1[x] != p1[y])
			{
				old_g1 = p1[x];
				old_g2 = p1[y];
				p1[x] = old_g2;
				p1[y] = old_g1;
				count++;
			}
		}
		//cout << "count=" << count <<",L="<<L<< endl;
	} while (count < L);
}

int rts(int* ss, double maxTime)
{
	int iter = 0, w = 0;
	double l;
	Objective = caculateObj(ss);
	TA = 16.73;										
	TB = 76.56;
	TC = 0.0031;
	ObjP = Objective;
	//ObjBest = Objective;
	TR = 1 / (TA * (ObjP / 1000) + TB) + TC;
	ThreshT = (1 + TR) * ObjP;
	cout << "objective=" << Objective << ",ThreshT=" << ThreshT << endl;
	initialMatrixNKAndObjClu(ss);
	child_update.cost = MAXNUM;
	//while ((clock() - StartTime) / CLOCKS_PER_SEC <= maxTime)
	while (iter < RTS_ITER)
	{
		l = tbe(LL, ss, maxTime);								//threshold-based exploration		
		//cout << "tbe,l=" << l << endl;
		l = dbi(ss, maxTime);									//descent-based improvement			
		//cout << "dbi,l=" << l << endl;
		if (l < ObjP)
		{
			ObjP = l;
			FinishTime = clock();
			TR = 1 / (TA * (ObjP / 1000) + TB) + TC;
			ThreshT = (1 + TR) * ObjP;
			w = 0;
		}
		else
			w += 1;
		if (l < child_update.cost)
		{
			child_update.cost = l;
			memcpy(child_update.p, ss, sizeof(int) * N);
		}
		/*if (w == W)
		{
			randomShake(0.1*N, ss);
			initialMatrixNKAndObjClu(ss);
			Objective = caculateObj(ss);
			ObjP = Objective;
			TR = 1 / (TA*(ObjP / 1000) + TB) + TC;
			ThreshT = (1 + TR)*ObjP;
			w = 0;
		}*/
		iter++;
		//printf("i=%d,ObjBest=%6f\n", iter++, ObjBest);
	}
	return Objective;
}
void initialPopulation(double maxTime)
{
	for (int i = 0; i < PopNum; i++)
	{
		randomConstruction(Pop[i].p);
		initialMatrixNKAndObjClu(Pop[i].p);
		dbi(Pop[i].p, 0.01 * maxTime);
		Pop[i].cost = caculateObj(Pop[i].p);
	}
}

//backbone crossover:基于最大匹配
void crossover()
{
	int parent1, parent2;
	int i, j, m, n;
	int unassLen;
	int ver1, ver2, sharedV, sharedMax;
	int index;
	int can1 = -1 , can2 = -1;
	for (i = 0; i < N; i++)
	{
		Child.p[i] = -1;
		flagV[i] = 0;
	}
	parent1 = rand() % PopNum;
	parent2 = rand() % PopNum;
	while (parent1 == parent2)
		parent2 = rand() % PopNum;
	for (i = 0; i < K; i++)
	{
		len[i] = 0;
		len2[i] = 0;
		flagC1[i] = 0;
		flagC2[i] = 0;
	}
	for (i = 0; i < N; i++)
	{
		int clu = Pop[parent1].p[i];
		arr1[clu][len[clu]++] = i;
		int clu2 = Pop[parent2].p[i];
		arr2[clu2][len2[clu2]++] = i;
	}
	//最大匹配
	sharedMax = 0;
	for (int iter = 0; iter < K; iter++)
	{
		sharedMax = 0;
		for (i = 0; i < K; i++)
		{
			if (flagC1[i] == 0)
			{
				for (m = 0; m < K; m++)
				{
					if (flagC2[m] == 0)
					{
						sharedV = 0;
						index = 0;
						for (j = 0; j < len[i]; j++)
						{
							ver1 = arr1[i][j];
							if (flagV[ver1] == 0)
							{
								for (n = index; n < len2[m]; n++)
								{
									ver2 = arr2[m][n];
									if (flagV[ver2] == 0)
									{
										if (ver1 == ver2)
											sharedV++;
										if (ver2 > ver1)
										{
											index = n;
											break;
										}
									}
								}

							}
						}
						if (sharedV > sharedMax)
						{
							sharedMax = sharedV;
							can1 = i;
							can2 = m;
						}
					}
				}
			}
		}
		match[can1] = can2;
		flagC1[can1] = 1;
		flagC2[can2] = 1;
		/*if (sharedMax > CluLen[iter])
		{
			cout << "please move the cluster to the proper position" << endl;
			getchar();
		}*/
		index = 0;
		for (int x1 = 0; x1 < len[can1]; x1++)
		{
			int ver = arr1[can1][x1];
			for (int x2 = index; x2 < len2[can2]; x2++)
			{
				int ver2 = arr2[can2][x2];
				if (ver == ver2)					//标记该顶点已被分配到子代
				{
					flagV[ver] = 1;
					Child.p[ver] = iter;
					index = x2;
					break;
				}
				if (ver2 > ver)
				{
					break;
					index = x2;
				}
			}
		}
	}
	for (i = 0; i < K; i++)
		len[i] = 0;
	unassLen = 0;
	for (i = 0; i < N; i++)
	{
		if (Child.p[i] != -1)
		{
			int clu = Child.p[i];
			arr1[clu][len[clu]++] = i;
		}
		else
		{
			unassV[unassLen] = i;
			addressUnaV[i] = unassLen;
			unassLen++;
		}
	}
	/*int most_size_clu, least_size_clu;
	if (N % K) {
		most_size_clu = N / K + 1;
		least_size_clu = N / K;
	}

	double nk = N / K;
	int size_each_clu = N / K;*/

	//随机分配剩余顶点

	for (int i = 0; i < K; i++)
		cout << "i=" << i << ", CluLen_i=" << CluLen[i] << ", len_i=" << len[i] << endl;
	cout << "unassLen=" << unassLen << endl;

	int len111 = 0;
	int* rand_unass = new int[unassLen + 1];
	bool* flag_unass = new bool[unassLen + 1];
	for (int i = 0; i < unassLen + 1; i++)
		flag_unass[i] = false;

	//shuffle the unassigned vertices at random
	while (len111 < unassLen)
	{
		int rn = rand() % unassLen;
		if (!flag_unass[rn])
		{
			rand_unass[len111++] = unassV[rn];
			flag_unass[rn] = true;
		}
	}
	// assign the unassigned vertex to the cluster with minimum number of points
	for (int i = 0; i < len111; i++)
	{
		int ver = rand_unass[i];
		int sel_k = -1;
		int size_min = N + 1;
		/*select the cluster with the minimum number of points*/
		for (int k = 0; k < K; k++)
		{
			if (len[k] < size_min)
			{
				size_min = len[k];
				sel_k = k;
			}
		}
		// assign the vertex to the cluster
		Child.p[ver] = sel_k;
		len[sel_k]++;
	}
	// in the  former version, I forgot this assignment, leading to wrong results, So sorry.
	for (int i = 0; i < K; i++)
		CluLen[i] = len[i];
	delete[]rand_unass;
	delete[]flag_unass;
}

//种群更新：淘汰cost最差的个体,且新生的子代与种群中个体都不一样
void  updatePopulation(int* ch, double cost)
{
	int i, flag_diff;
	double maxCost = -MAXNUM;
	int select = -1;
	for (i = 0; i < PopNum; i++)
	{
		if (Pop[i].cost > maxCost)
		{
			maxCost = Pop[i].cost;
			select = i;
		}
	}
	flag_diff = 0;
	for (i = 0; i < PopNum; i++)
		if (fabs(Pop[i].cost - cost) < DEVIATION)
			flag_diff = 1;
	if (cost < maxCost && flag_diff == 0)
	{
		for (i = 0; i < N; i++)
			Pop[select].p[i] = ch[i];
		Pop[select].cost = cost;
	}
}

//计算两个个体之间的相似度
int Calculate_Sim_Between_Two(int* p1, int* p)
{

	int sum = 0, count, i, sum1, delta;
	int m1, m2, m3, m4, target1 = -1, target2 = -1;
	memset(flagC1, 0, sizeof(int) * K);
	memset(flagC2, 0, sizeof(int) * K);
	memset(len, 0, sizeof(int) * K);
	memset(len2, 0, sizeof(int) * K);
	for (i = 0; i < N; i++)
	{
		arr1[p1[i]][len[p1[i]]++] = i;
		arr2[p[i]][len2[p[i]]++] = i;
	}
	for (count = 0; count < K; count++)
	{
		delta = 0;
		for (m1 = 0; m1 < K; m1++)
		{

			if (flagC1[m1] != 1)
			{
				for (m2 = 0; m2 < K; m2++)
				{
					if (flagC2[m2] != 1)																		//表示染色体2 当前色族未被选中
					{
						sum1 = 0;
						for (m3 = 0; m3 < len[m1]; m3++)
						{
							for (m4 = 0; m4 < len2[m2]; m4++)
							{
								if (arr1[m1][m3] == arr2[m2][m4])
								{
									sum1++;
									break;
								}
							}

						}
						//cout << "sum1=" << sum1 <<" m1="<<m1<<" m2="<<m2<< endl;
						if (sum1 > delta)												//最大匹配
						{
							delta = sum1;
							target1 = m1;
							target2 = m2;
						}
					}

				}
			}
		}
		flagC1[target1] = 1;
		flagC2[target2] = 1;
		sum += delta;
		//diversity = MaxVtx - sum;
		//sumDiversity += diversity;
	}
	return sum;
}

//计算种群相似性
void compute_similarity()
{
	int i, j, sum_sim = 0;
	double similarity;
	for (i = 0; i < PopNum; i++)
		for (j = i + 1; j < PopNum; j++)
			sum_sim += Calculate_Sim_Between_Two(Pop[i].p, Pop[j].p);
	similarity = double(sum_sim) / (N * PopNum * (PopNum - 1) / 2);
	cout << "similarity=" << similarity << endl;
}

void memetic(double maxTime)
{
	ObjBest = MAXNUM;
	initialPopulation(maxTime);
	int iter = 0;
	while ((clock() - StartTime) / CLOCKS_PER_SEC <= maxTime)
	{
		crossover();
		rts(Child.p, maxTime);
		updatePopulation(child_update.p, child_update.cost);
		printf("generations=%d,objbest=%6f,spend time=%f\n", iter++, ObjBest, (clock() - StartTime) / CLOCKS_PER_SEC);
	}
	compute_similarity();
}

void freeMemory1()
{
	delete[] CluLen;
	delete[] ObjClu;
	for (int i = 0; i < N; i++)
		delete[] MatrixNK[i];
	delete[] MatrixNK;
	delete[] Neighbors;
	delete[] flagN;
	delete[] flagK;
	delete[] randN;
	delete[] randK;
	for (int i = 0; i < K; i++)
	{
		delete[] arr1[i];
		delete[] arr2[i];
	}
	delete[] arr1;
	delete[] arr2;
	delete[] len;
	delete[] len2;
	delete[] match;
	delete[] flagC1;
	delete[] flagV;
	delete[] unassV;
	delete[] addressUnaV;
}

void freeMemory()
{
	delete[] BestS;
	delete[] Child.p;
	for (int i = 0; i < N; i++)
	{
		delete[] Point[i];
		delete[] DisSquare[i];
	}
	delete[] Point;
	delete[] DisSquare;

}

/*
* param sol: to be verified solution, n: the number of points, k: the number of clusters
* verify the obtained solution is feasible or not
*/
bool verify_sol(int* sol, double sol_f, int n, int k)
{
	int* len_arr = new int[k];
	double obj = 0;
	bool feasible = true;
	memset(len_arr, 0, sizeof(int) * k);
	for (int i = 0; i < n; i++)
		len_arr[sol[i]]++;
	if (n % k) {
		for (int i = 0; i < k; i++) {
			if (feasible) {
				for (int j = i + 1; j < k; j++) {
					if (len_arr[i] - len_arr[j] > 1 || len_arr[j] - len_arr[i] > 1) {
						feasible = false;
						break;
					}
				}
			}
			else
				break;
		}
	}
	else {
		for (int i = 0; i < k; i++) {
			if (feasible) {
				for (int j = i + 1; j < k; j++) {
					if (len_arr[i] != len_arr[j]) {
						feasible = false;
						break;
					}
				}
			}
			else
				break;
		}
	}
	if (!feasible) {
		cout << "in verify_sol method, an error is checked!!!!!!!!!!!!!!!!!!!!" << endl;
		for (int i = 0; i < k; i++)
			cout << "i=" << i << ", len_i=" << len_arr[i] << endl;
		getchar();
	}
	for (int i = 0; i < K; i++)
		CluLen[i] = len_arr[i];
	obj = caculateObj(sol);
	if (fabs(obj - sol_f) > 1.0e-3) {
		cout << "in verify_sol method, an error is checked!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "obj=" << obj << ", sol_f=" << sol_f << endl;
		getchar();
	}

	delete[]len_arr;
	return feasible;
}

int main(int argc, char* argv[])
{	
	if(argc < 8){
		cout << "MA_BMSSC usage: input_file n d k time seed output_res_file" << endl;
		cout << "(where input_file is the instance name, n is the numbe of points of the instance, d is the dimension of the point, k is the number of clusterrs, time is the cutoff of an execution, seed is the random seed, such as 1, 2, ..., 10., output_res_file is a file used to store the running information)" << endl;
		exit(-1);
	}
	//the following are used parameters
	PopNum = 15;
	LL = 5;
	RTS_ITER = 50;
		
	ofstream resultsFile;	
	ObjBest = MAXNUM;

	char *input_file = argv[1];
	int n = atoi(argv[2]);
	int d = atoi(argv[3]);
	int k = atoi(argv[4]);
	double time = atof(argv[5]);
	int seed = atoi(argv[6]);
	char *output_res_file = argv[7];

	srand(seed);	
	resultsFile.open(output_res_file, ofstream::app);
		
	N = n;
	D = d;
	initialing(input_file);	
	K = k;
	allocateMemory();
	buildNeighbors();
	
	StartTime = clock();
	memetic(time);
	verify_sol(BestS, ObjBest, N, K);
	Runtime = (FinishTime - StartTime) / CLOCKS_PER_SEC;
	cout << endl << setprecision(6) << scientific << "Objective Function value: " << ObjBest << " in " << setprecision(2) << fixed << Runtime << " seconds" << endl;
			
	resultsFile << input_file << " ,K=" << K << " ,bestV=" << setprecision(6) << scientific << ObjBest;
	resultsFile << " ,run_time=" << setprecision(2) << fixed << Runtime << endl;
	//cout << endl;
	freeMemory1();	
	freeMemory();
	
}
