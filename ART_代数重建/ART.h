# pragma once
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# define pi 3.1415926 //定义圆周率

typedef struct node//定义射线和矩形方格的交点的结构体
{
	double result;//参数方程的t
	int XI;//与第几个X(i)相交
	int YJ;//与第几个Y(j)相交
}NodeT;

//引用Main函数中的全局变量
extern int m;
extern int n;
extern double rmax;

//计算矩形方格每条格子的方程
void ComputeXOY(double *X, double *Y)
{
	int i;//循环变量
	double dr;//探测器单元的间隔
	//计算dr
	dr = 2.0 * rmax / m;
	for (i = 0; i <= m; i++)
	{
		X[i] = -rmax + i * dr;//横轴
		Y[i] = rmax - i * dr;//纵轴
	}
}

//计算给定射线与X 、Y轴的交点
void LineMethod(double r_i, double angle, double *X, double *Y, struct node *X_T, struct node *Y_T, int &lenX, int &lenY)
{
	double Y_M = 5.0 * rmax; //表示所有射线的纵坐标
	double x_m; //存储在角度angle下第i条射线的横坐标
	double y_m; //存储在角度angle下第i条射线的纵坐标
	double temp, temp2;//中间变量
	int i, j, k = 0;//循环变量

	//计算给定角度angle下射线的x_m和y_m
	x_m = r_i * cos(angle) - Y_M * sin(angle);
	y_m = r_i * sin(angle) + Y_M * cos(angle);

	//printf("y_m=%lf,cos(angle)=%lf\n",y_m, cos(angle));
	//计算与X轴交点
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
	//计算与Y轴交点
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

// 插入排序法对求得的t进行从小到大排序
void ChooseSort(struct node *X_T, struct node *Y_T, struct node *sort, int lenX, int lenY) 
{
	int i, j;//循环变量
	struct node temp;//中间变量
	//合并
	for (i = 0; i < lenX; i++) 
	{
		sort[i] = X_T[i];
	}
	for (i = 0; i < lenY; i++) 
	{
		sort[lenX + i] = Y_T[i];
	}
	//排序
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
//返回两个数中最小的一个
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
//计算交线长
void ComputeLine(double **R, struct node *sort, double r_i, double angle, int lenX, int lenY)
{
	double Y_M = 5.0 * rmax; //表示所有射线的纵坐标
	double x_m; //存储在角度angle下第i条射线的横坐标
	double y_m; //存储在角度angle下第i条射线的纵坐标
	int i;//循环变量
	double temp;//中间变量
	int row, col;//中间变量
	double dr;//探测器单元的间隔

	//计算dr
	dr = 2.0 * rmax / m;

	//计算给定角度angle下射线的x_m和y_m
	x_m = r_i * cos(angle) - Y_M * sin(angle);
	y_m = r_i * sin(angle) + Y_M * cos(angle);

	//printf("运行到此处1！\n");
	//计算交线长
	for (i = 0; i < lenX + lenY - 1; i++)
	{
		//printf("执行到此处！\n");
		if (sort[i].XI == -1 && sort[i+1].XI == -1) 
		{
			temp = x_m + sort[i].result * sin(angle);
			//r[i]=-rmax + i * dr;
			temp = (temp + rmax) / dr;
			//printf("temp=%lf\n",temp);
			col = (int)temp;
			row = Min(sort[i].YJ, sort[i+1].YJ);
			//printf("row=%d,col=%d\n", row, col);
			//printf("执行到1处\n");
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
			//printf("执行到2处\n");
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
				//printf("执行到3处\n");
				R[row][col] = sort[i + 1].result - sort[i].result;
			}
			else 
			{
				row = sort[i].YJ - 1;
				col = sort[i + 1].XI - 1;
				//printf("row=%d,col=%d\n", row, col);
				//printf("执行到4处\n");
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
				//printf("执行到5处\n");
				R[row][col] = sort[i + 1].result - sort[i].result;
			}
			else 
			{
				row = sort[i + 1].YJ;
				col = sort[i].XI;
				//printf("row=%d,col=%d\n", row, col);
				//printf("执行到6处\n");
				R[row][col] = sort[i + 1].result - sort[i].result;
			}
		}
	}
	//printf("执行\n");
}

//在角度为0°和90°的时候的交线长
void ComputeZero(double **R,double angle, int number)
{
	int i;//循环变量
	double dr;//探测器单元的间隔

	//计算dr
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

//一次迭代过程
void OneIterator(double **R, double **P, double **X0, double **XX, int r_i, int theta_i) 
{
	int i, j;//循环变量
	double temp1, temp2;//中间变量
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
	//更新X
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			XX[i][j] = X0[i][j] + temp2 * R[i][j] / temp1;
		}
	}
}

//初始化R
void initR(double **R) 
{
	int i, j;//循环变量
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
	FILE *fp;//定义文件指针
	double *X, *Y;//定义指针数组X和Y存储矩形方格的坐标
	double **R;//定义矩阵R存储交线长
	struct node *X_T, *Y_T;//存储与X、Y的交点
	struct node *sort;//存储排序后的交点
	int lenX, lenY;
	int i, j, k; //循环变量

	/*****************开辟空间***********************/
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

	/********************初始化**********************/
	for (i = 0; i <= m; i++)
	{
		X_T[i].result = -1;
		X_T[i].XI = -1;
		X_T[i].YJ = -1;

		Y_T[i].result = -1;
		Y_T[i].XI = -1;
		Y_T[i].YJ = -1;
	}
	//初始化R
	initR(R);

	/******************************************************************************/
	//计算矩形方格的坐标
	ComputeXOY(X, Y);
	k = 0;
	while (k < 1) 
	{
		//printf("-----------第%d次更新--------\n",k+1);
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < m; j++)
			{
				//在角度为0°和90°的时候的交线长
				if (theta[i] == 0 || theta[i] == pi / 2.0)
				{
					ComputeZero(R, theta[i], j);
				}
				else
				{
					//计算交点
					LineMethod(r[j], theta[i], X, Y, X_T, Y_T, lenX, lenY);
					//排序
					ChooseSort(X_T, Y_T, sort, lenX, lenY);
					//printf("i=%d,j=%d\n", i, j);
					//计算交线长
					ComputeLine(R, sort, r[j], theta[i], lenX, lenY);
				}
				//迭代一次
				OneIterator(R, P, X0, XX, j, i);
				//每次更新之后把更新之后的值赋给初值
				X0 = XX;
				//每次更新之后初始化R
				initR(R);
			}
		}
		k = k + 1;
	}
	//将结果写入文件
	//将Radon变换的数据写入文件当中
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

	/****************释放空间**********************/
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