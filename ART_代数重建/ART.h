# pragma once
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# define pi 3.1415926 //����Բ����

typedef struct node//�������ߺ;��η���Ľ���Ľṹ��
{
	double result;//�������̵�t
	int XI;//��ڼ���X(i)�ཻ
	int YJ;//��ڼ���Y(j)�ཻ
}NodeT;

//����Main�����е�ȫ�ֱ���
extern int m;
extern int n;
extern double rmax;

//������η���ÿ�����ӵķ���
void ComputeXOY(double *X, double *Y)
{
	int i;//ѭ������
	double dr;//̽������Ԫ�ļ��
	//����dr
	dr = 2.0 * rmax / m;
	for (i = 0; i <= m; i++)
	{
		X[i] = -rmax + i * dr;//����
		Y[i] = rmax - i * dr;//����
	}
}

//�������������X ��Y��Ľ���
void LineMethod(double r_i, double angle, double *X, double *Y, struct node *X_T, struct node *Y_T, int &lenX, int &lenY)
{
	double Y_M = 5.0 * rmax; //��ʾ�������ߵ�������
	double x_m; //�洢�ڽǶ�angle�µ�i�����ߵĺ�����
	double y_m; //�洢�ڽǶ�angle�µ�i�����ߵ�������
	double temp, temp2;//�м����
	int i, j, k = 0;//ѭ������

	//��������Ƕ�angle�����ߵ�x_m��y_m
	x_m = r_i * cos(angle) - Y_M * sin(angle);
	y_m = r_i * sin(angle) + Y_M * cos(angle);

	//printf("y_m=%lf,cos(angle)=%lf\n",y_m, cos(angle));
	//������X�ύ��
	for (i = 0; i <= m; i++)
	{
		temp = (X[i] - x_m) / sin(angle);
		//printf("temp=%lf\t",temp);
		temp2 = y_m - temp * cos(angle);
		//printf("temp2 = %.10lf\n",temp2);
		if (temp2 >= Y[m] && temp2 <= Y[0])
		{
			X_T[k].result = temp;
			X_T[k].XI = i;
			//printf("X_T[%d].result=%lf,X_T[%d].XI=%d\n", k,X_T[k].result, k, X_T[k].XI);
			k = k + 1;
		}
	}
	lenX = k;

	k = 0;
	//������Y�ύ��
	for (j = 0; j <= m; j++)
	{
		temp = (y_m - Y[j]) / cos(angle);
		temp2 = x_m + temp * sin(angle);
		if (temp2 >= X[0] && temp2 <= X[m])
		{
			Y_T[k].result = temp;
			Y_T[k].YJ = j;
			//printf("Y_T[%d].result=%lf,Y_T[%d].XI=%d\n", k, Y_T[k].result, k, Y_T[k].XI);
			k = k + 1;
		}
	}
	lenY = k;
}

// �������򷨶���õ�t���д�С��������
void ChooseSort(struct node *X_T, struct node *Y_T, struct node *sort, int lenX, int lenY) 
{
	int i, j;//ѭ������
	struct node temp;//�м����
	//�ϲ�
	for (i = 0; i < lenX; i++) 
	{
		sort[i] = X_T[i];
	}
	for (i = 0; i < lenY; i++) 
	{
		sort[lenX + i] = Y_T[i];
	}
	//����
	for (i = 0; i < lenX + lenY - 1; i++)
	{
		temp = sort[i + 1];
		for (j = i; j >= 0; j--)
		{
			if (sort[j].result > temp.result)
			{
				sort[j + 1] = sort[j];
				sort[j] = temp;
			}
		}
	}
}
//��������������С��һ��
int Min(int num1, int num2) 
{
	if (num1 > num2) 
	{
		return num2;
	}
	else 
	{
		return num1;
	}
}
//���㽻�߳�
void ComputeLine(double **R, struct node *sort, double r_i, double angle, int lenX, int lenY)
{
	double Y_M = 5.0 * rmax; //��ʾ�������ߵ�������
	double x_m; //�洢�ڽǶ�angle�µ�i�����ߵĺ�����
	double y_m; //�洢�ڽǶ�angle�µ�i�����ߵ�������
	int i;//ѭ������
	double temp;//�м����
	int row, col;//�м����
	double dr;//̽������Ԫ�ļ��

	//����dr
	dr = 2.0 * rmax / m;

	//��������Ƕ�angle�����ߵ�x_m��y_m
	x_m = r_i * cos(angle) - Y_M * sin(angle);
	y_m = r_i * sin(angle) + Y_M * cos(angle);

	//printf("���е��˴�1��\n");
	//���㽻�߳�
	for (i = 0; i < lenX + lenY - 1; i++)
	{
		//printf("ִ�е��˴���\n");
		if (sort[i].XI == -1 && sort[i+1].XI == -1) 
		{
			temp = x_m + sort[i].result * sin(angle);
			//r[i]=-rmax + i * dr;
			temp = (temp + rmax) / dr;
			//printf("temp=%lf\n",temp);
			col = (int)temp;
			row = Min(sort[i].YJ, sort[i+1].YJ);
			//printf("row=%d,col=%d\n", row, col);
			//printf("ִ�е�1��\n");
			//printf("x_m=%lf,sort[i].result=%lf,sort[i].YJ=%d,col=%d\n", x_m, sort[i].result,sort[i].YJ,col);
			R[row][col] = sort[i+1].result - sort[i].result;
		}
		else if (sort[i].YJ ==-1 && sort[i+1].YJ == -1) 
		{
			temp = y_m - sort[i].result * cos(angle);
			//Y[i] = rmax - i * dr;
			temp = (rmax - temp) / dr;
			row = (int)temp;
			col = Min(sort[i].XI, sort[i + 1].XI);
			//printf("row=%d,col=%d\n", row, col);
			//printf("ִ�е�2��\n");
			R[row][col] = sort[i + 1].result - sort[i].result;
		}
		else if (sort[i].XI == -1 && sort[i+1].YJ == -1) 
		{
			if (fabs(sort[i + 1].result - sort[i].result) < 1e-8)
			{
				continue;
			}
			if ( angle < pi/2) 
			{
				row = sort[i].YJ;
				col = sort[i + 1].XI - 1;
				//printf("row=%d,col=%d\n",row,col);
				//printf("ִ�е�3��\n");
				R[row][col] = sort[i + 1].result - sort[i].result;
			}
			else 
			{
				row = sort[i].YJ - 1;
				col = sort[i + 1].XI - 1;
				//printf("row=%d,col=%d\n", row, col);
				//printf("ִ�е�4��\n");
				R[row][col] = sort[i + 1].result - sort[i].result;
			}
		}
		else if (sort[i].YJ == -1 && sort[i+1].XI == -1) 
		{
			if (fabs(sort[i + 1].result - sort[i].result) < 1e-8)
			{
				continue;
			}
			if (angle < pi/2) 
			{
				row = sort[i + 1].YJ - 1;
				col = sort[i].XI;
				//printf("row=%d,col=%d\n", row, col);
				//printf("ִ�е�5��\n");
				R[row][col] = sort[i + 1].result - sort[i].result;
			}
			else 
			{
				row = sort[i + 1].YJ;
				col = sort[i].XI;
				//printf("row=%d,col=%d\n", row, col);
				//printf("ִ�е�6��\n");
				R[row][col] = sort[i + 1].result - sort[i].result;
			}
		}
	}
	//printf("ִ��\n");
}

//�ڽǶ�Ϊ0���90���ʱ��Ľ��߳�
void ComputeZero(double **R,double angle, int number)
{
	int i;//ѭ������
	double dr;//̽������Ԫ�ļ��

	//����dr
	dr = 2.0 * rmax / m;
	if (angle == 0) 
	{
		for (i = 0; i < m; i++)
		{
			R[i][number] = dr;
		}
	}
	else 
	{
		for (i = 0; i < m; i++)
		{
			R[number][i] = dr;
		}
	}
}

//һ�ε�������
void OneIterator(double **R, double **P, double **X0, double **XX, int r_i, int theta_i) 
{
	int i, j;//ѭ������
	double temp1, temp2;//�м����
	temp1 = 0; temp2 = 0;
	for (i = 0; i < m; i++) 
	{
		for (j = 0; j < m; j++) 
		{
			temp1 = temp1 + R[i][j] * R[i][j];
			temp2 = temp2 + R[i][j] * X0[i][j];
		}
	}
	temp2 = P[theta_i][r_i] - temp2;
	//����X
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			XX[i][j] = X0[i][j] + temp2 * R[i][j] / temp1;
		}
	}
}

//��ʼ��R
void initR(double **R) 
{
	int i, j;//ѭ������
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			R[i][j] = 0;
		}
	}
}
void Iterator(double *r, double *theta, double **P, double **X0, double **XX)
{
	FILE *fp;//�����ļ�ָ��
	double *X, *Y;//����ָ������X��Y�洢���η��������
	double **R;//�������R�洢���߳�
	struct node *X_T, *Y_T;//�洢��X��Y�Ľ���
	struct node *sort;//�洢�����Ľ���
	int lenX, lenY;
	int i, j, k; //ѭ������

	/*****************���ٿռ�***********************/
	X = (double *)malloc((m + 1) * sizeof(double));
	Y = (double *)malloc((m + 1) * sizeof(double));
	R = (double **)malloc(m * sizeof(double *));
	for (i = 0; i < m; i++)
	{
		R[i] = (double *)malloc(m * sizeof(double));
	}
	X_T = (struct node *)malloc((m + 1) * sizeof(struct node));
	Y_T = (struct node *)malloc((m + 1) * sizeof(struct node));
	sort = (struct node *)malloc(m*m * sizeof(struct node));

	/********************��ʼ��**********************/
	for (i = 0; i <= m; i++)
	{
		X_T[i].result = -1;
		X_T[i].XI = -1;
		X_T[i].YJ = -1;

		Y_T[i].result = -1;
		Y_T[i].XI = -1;
		Y_T[i].YJ = -1;
	}
	//��ʼ��R
	initR(R);

	/******************************************************************************/
	//������η��������
	ComputeXOY(X, Y);
	k = 0;
	while (k < 1) 
	{
		//printf("-----------��%d�θ���--------\n",k+1);
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < m; j++)
			{
				//�ڽǶ�Ϊ0���90���ʱ��Ľ��߳�
				if (theta[i] == 0 || theta[i] == pi / 2.0)
				{
					ComputeZero(R, theta[i], j);
				}
				else
				{
					//���㽻��
					LineMethod(r[j], theta[i], X, Y, X_T, Y_T, lenX, lenY);
					//����
					ChooseSort(X_T, Y_T, sort, lenX, lenY);
					//printf("i=%d,j=%d\n", i, j);
					//���㽻�߳�
					ComputeLine(R, sort, r[j], theta[i], lenX, lenY);
				}
				//����һ��
				OneIterator(R, P, X0, XX, j, i);
				//ÿ�θ���֮��Ѹ���֮���ֵ������ֵ
				X0 = XX;
				//ÿ�θ���֮���ʼ��R
				initR(R);
			}
		}
		k = k + 1;
	}
	//�����д���ļ�
	//��Radon�任������д���ļ�����
	fp = fopen("D://Matlab/workspace/CT-2019-7/chongjian.txt", "w");
	if (!fp)
	{
		printf("can not open the file!\n");
		exit(-1);
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			fprintf(fp, "%15.9lf", XX[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("Success!\n");
	/*****************************************************************************/

	/****************�ͷſռ�**********************/
	free(X);
	free(Y);
	for (i = 0; i < m; i++)
	{
		free(R[i]);
	}
	free(R);
	free(X_T);
	free(Y_T);
	free(sort);
}