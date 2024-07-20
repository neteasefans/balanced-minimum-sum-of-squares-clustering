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
using namespace std;

double StartTime, FinishTime, Runtime;
int Runs;		//ÿ���������д���

int* BestS;
int* CluLen;	//ÿ��cluster��size��swapǰ��size����
double Objective;
double ObjBest;

//N,K,D��main���������
int N;		//the number of vertices
int K;		//the number of clusters
int D;		//dimensions
int Type;		//��������
double** Point;	//�������꣬����ά��ΪD
double** MatrixNK;
double* ObjClu;	//ÿ��cluster��Ŀ��ֵ��û�г��Ը�cluster��size
double** DisSquare;	//ŷ����ÿռ�����������ŷʽ�����ƽ��

//for the memetic algorithm
typedef struct Population {
	int* p;
	double cost;
}Population;
Population* Pop;
Population Child;		//�Ӵ�
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
int RTS_ITER;				//ÿ����RTS��������
Neighborhood* Neighbors;
int* randN;
int* randK;
int* flagN;
int* flagK;

#define PopNum 5										//��Ⱥ��С
#define DEBUG
#define MAXNUM 99999999999999
#define NUM1	50
#define DEVIATION	0.000001
#define DATASETNUM 16
#define CLUTERNUM 10
#define RUNS 10
double MaxTimes[DATASETNUM][CLUTERNUM][RUNS];
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
		while (ss >> dd)			//�����ļ�������ʽ��������ÿ���������Ʊ����ָ�
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

void readTimeFile(string file)
{
	string line;
	stringstream ss;
	fstream fs(file);
	int col1, col2, col3;
	double dd;
	col1 = 0;
	//col1 = 14;
	col2 = 0;

	cout << "file=" << file << endl;
	if (!fs)
	{
		cout << "error open, File_Name: " << file << endl;
		getchar();
		exit(1);
	}
	while (getline(fs, line))
	{
		col3 = 0;
		ss << line;
		while (ss >> dd)			//�����ļ�������ʽ��������ÿ���������Ʊ����ָ�
		{
			MaxTimes[col1][col2][col3++] = dd;
		}
		col2++;
		if (col2 % CLUTERNUM == 0)
		{
			col1++;
			col2 = 0;
		}
		ss.clear();
	}
#ifdef DEBUG1			
	for (int i = 0; i < col1; i++)
	{
		for (int k = 0; k < CLUTERNUM; k++)
		{
			for (int j = 0; j < RUNS; j++)
				cout << "i=" << i << ",k=" << k << ",j=" << j << "," << MaxTimes[i][k][j] << " " << endl;
			cout << endl;
		}
	}
	cout << endl;
#endif
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

//����Ŀ�꺯��ֵ
double caculateObj(int* ss)
{
	double dd;
	double obj = 0;
	int dimens = 0;
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

//��֮ǰ��xiangjing.lai�����е�gamma��������ͬ
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

//����randomshake��������������
void buildNeighbors()
{
	int i, j, g;
	int count;
	int SN = N * (N - 1) / 2 + N * K;
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

//����Ŷ���ʹ��one-move��swap move����
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
	int LL = 10;
	int W = 20;
	Objective = caculateObj(ss);
	TA = 16.73;										//ta,tb,tc ������������ôȷ����
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

//backbone crossover:�������ƥ��
void crossover()
{
	int parent1, parent2;
	int i, j, m, n;
	int unassLen;
	int ver1, ver2, sharedV, sharedMax;
	int index;
	int can1, can2;
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
	//���ƥ��
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
				if (ver == ver2)					//��Ǹö����ѱ����䵽�Ӵ�
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

	//�������ʣ�ඥ��

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

//��Ⱥ���£���̭cost���ĸ���,���������Ӵ�����Ⱥ�и��嶼��һ��
void  updatePopulation(int* ch, double cost)
{
	int i, flag_diff;
	double maxCost = -MAXNUM;
	int select;
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

//������������֮������ƶ�
int Calculate_Sim_Between_Two(int* p1, int* p)
{

	int sum = 0, count, i, sum1, delta;
	int m1, m2, m3, m4, target1, target2;
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
					if (flagC2[m2] != 1)																		//��ʾȾɫ��2 ��ǰɫ��δ��ѡ��
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
						if (sum1 > delta)												//���ƥ��
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

//������Ⱥ������
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
	if (fabs(obj - sol_f) > 1.0e-6) {
		cout << "in verify_sol method, an error is checked!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "obj=" << obj << ", sol_f=" << sol_f << endl;
		getchar();
	}

	delete[]len_arr;
	return feasible;
}

int main(int argc, char* argv[])
{
	string instances[DATASETNUM];
	int nClusters[DATASETNUM];
	int nPoints[DATASETNUM];
	int nDimensions[DATASETNUM];
	string instanceName[DATASETNUM];
	int nCluterK[DATASETNUM][CLUTERNUM];
	int cluK1[CLUTERNUM] = { 2, 3, 4, 6, 7, 10, 11, 13, 15, 20 };
	string PATH = "F:\\clustering problem\\BMSSC(balanced minimum sum-of-squares of clustering)\\data2\\";
	string timeFile = PATH + "time1.txt";


	instances[0] = PATH + "iris.txt";	//150
	nPoints[0] = 150; nDimensions[0] = 4; nClusters[0] = 6;
	instanceName[0] = "iris";

	instances[1] = PATH + "wine.txt";	//178
	nPoints[1] = 178; nDimensions[1] = 13; nClusters[1] = 3;
	instanceName[1] = "wine";

	instances[2] = PATH + "glass.txt";	//214
	nPoints[2] = 214; nDimensions[2] = 9; nClusters[2] = 7;
	instanceName[2] = "glass";

	instances[3] = PATH + "thyroid.txt";//215
	nPoints[3] = 215; nDimensions[3] = 5; nClusters[3] = 3;
	instanceName[3] = "thyroid";

	instances[4] = PATH + "ionosphere.txt"; //351
	nPoints[4] = 351; nDimensions[4] = 34; nClusters[4] = 2;
	instanceName[4] = "ionosphere";

	instances[5] = PATH + "libra.txt";//360
	nPoints[5] = 360; nDimensions[5] = 90; nClusters[5] = 15;
	instanceName[5] = "libra";

	instances[6] = PATH + "user_knowledge.txt";//403
	nPoints[6] = 403; nDimensions[6] = 5; nClusters[6] = 4;
	instanceName[6] = "user_knowledge";

	instances[7] = PATH + "body.txt";//507
	nPoints[7] = 507; nDimensions[7] = 5; nClusters[7] = 2;
	instanceName[7] = "body";

	instances[8] = PATH + "water.txt";//527
	nPoints[8] = 527; nDimensions[8] = 38; nClusters[8] = 13;
	instanceName[8] = "water";

	instances[9] = PATH + "breast.txt";//569
	nPoints[9] = 569; nDimensions[9] = 30; nClusters[9] = 2;
	instanceName[9] = "breast";

	instances[10] = PATH + "synthetic.txt";//600
	nPoints[10] = 600; nDimensions[10] = 60; nClusters[10] = 6;
	instanceName[10] = "synthetic";

	instances[11] = PATH + "vehicle.txt";//846
	nPoints[11] = 846; nDimensions[11] = 18; nClusters[11] = 6;
	instanceName[11] = "vehicle";

	instances[12] = PATH + "vowel.txt";//990
	nPoints[12] = 990; nDimensions[12] = 10; nClusters[12] = 11;
	instanceName[12] = "vowel";

	instances[13] = PATH + "yeast.txt"; //1484
	nPoints[13] = 1484; nDimensions[13] = 8; nClusters[13] = 10;
	instanceName[13] = "yeast";

	instances[14] = PATH + "multiple.txt";//2000
	nPoints[14] = 2000; nDimensions[14] = 240; nClusters[14] = 7;
	instanceName[14] = "multiple";

	instances[15] = PATH + "image.txt";//2310
	nPoints[15] = 2310; nDimensions[15] = 19; nClusters[15] = 7;
	instanceName[15] = "image";
	srand(time(NULL));
	ofstream resultsFile;
	ofstream valuesFile;
	resultsFile.open("E:\\resultados_MA_RTS.txt", ofstream::app);
	valuesFile.open("E:\\solucoes_MA_RTS.txt", ofstream::app);
	Runs = 10;
	RTS_ITER = 100;
	double bestValue = MAXNUM;
	double avgValue = 0;
	double avgTime = 0;
	double bestTime;
	readTimeFile(timeFile);
	ObjBest = MAXNUM;
	for (int i = 0; i < DATASETNUM; i++)
		for (int j = 0; j < CLUTERNUM; j++)
			nCluterK[i][j] = cluK1[j];
	for (int i = 0; i < 16; i++)
	{
		N = nPoints[i];
		D = nDimensions[i];
		//K = nClusters[i];
		initialing(instances[i]);
		for (int k = 0; k < CLUTERNUM; k++)
		{
			K = nCluterK[i][k];
			allocateMemory();
			buildNeighbors();
			bestValue = MAXNUM;
			bestTime = MAXNUM;
			avgValue = 0;
			avgTime = 0;
			cout << "=================================================INSTANCE " << i + 1 << "=================================================" << endl;
			cout << "Location: " << instances[i] << endl;
			cout << "Clusters: " << K << endl;
			valuesFile << instanceName[i] << ":";
			for (int j = 0; j < Runs; j++)
			{
				cout << "------------------------------------- Execution " << j + 1 << "-----------------------------------------" << endl;
				cout << "maxTime = " << MaxTimes[i][k][j] << endl;
				StartTime = clock();
				memetic(MaxTimes[i][k][j]);
				Runtime = (FinishTime - StartTime) / CLOCKS_PER_SEC;
				cout << endl << setprecision(6) << scientific << "Objective Function value: " << ObjBest << " in " << setprecision(2) << fixed << Runtime << " seconds" << endl;
				if (ObjBest < bestValue)
					bestValue = ObjBest;
				if (Runtime < bestTime)
					bestTime = Runtime;
				avgValue += ObjBest;
				avgTime += Runtime;
				valuesFile << setprecision(6) << scientific << ObjBest << ";";
				verify_sol(BestS, ObjBest, N, K);
			}
			resultsFile << instanceName[i] << ":bestV=" << setprecision(6) << scientific << bestValue;
			resultsFile << ",avgV=" << setprecision(6) << scientific << avgValue / Runs;
			resultsFile << ",bestTime=" << setprecision(2) << fixed << bestTime;
			resultsFile << ",avgTime=" << setprecision(2) << fixed << avgTime / Runs << ",K=" << K << endl;
			valuesFile << ",K=" << K << endl;
			cout << endl;
			freeMemory1();
		}
		freeMemory();
	}
}