# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# include "Radon.h"
# include "ART.h"

int m, n; //距离r和角度theta的离散化数量
double rmax; //扫描范围

void main() 
{

	time_t start, end;
	double cost;
	time(&start);
	
	double *r, *theta;//定义指针r和theta存储离散化后的r和theta
	double **P;//定义指针数组P存储Radon变换的结果,矩阵P中每个元素表示(theta,r)
	double **X0, **XX;//存储初始值和更新值
	double x0, y0, r0, miou;/*x0、y0是小圆的圆心坐标，r0是小圆半径，miou是小圆的吸收率*/
	int i, j;//循环变量

	/**************输入已知量*****************************/
	printf("请输入小圆的圆心坐标和半径，元素之间用空格分开：\n");
	scanf_s("%lf%lf%lf", &x0, &y0, &r0);
	printf("请输入小圆的吸收率：\n");
	scanf_s("%lf", &miou);
	printf("请分别输入距离r和角度theta的离散化数量，元素之间用空格分开：\n");
	scanf_s("%d%d",&m,&n);

	/************开辟存储空间*****************************/
	r = (double *)malloc(m * sizeof(double));
	theta = (double *)malloc(n * sizeof(double));
	P = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++) 
	{
		P[i] = (double *)malloc(m * sizeof(double));
	}
	X0 = (double **)malloc(m * sizeof(double *));
	for (i = 0; i < m; i++)
	{
		X0[i] = (double *)malloc(m * sizeof(double));
	}
	XX= (double **)malloc(m * sizeof(double *));
	for (i = 0; i < m; i++)
	{
		XX[i] = (double *)malloc(m * sizeof(double));
	}


	//计算扫描范围
	rmax = sqrt(x0 * x0 + y0 * y0) + 3 * r0;
	//初始化X0
	for (i = 0; i < m; i++) 
	{
		for (j = 0; j < m; j++) 
		{
			X0[i][j] = 0;
			XX[i][j] = 0;
		}
	}
	/*****************************************************/
	/*离散化r和theta*/
	ComputeR_Theta(r, theta);
	/*将theta随机化*/
	randperm(theta, n);
	/*Radon变换*/
	RadonMethod(r, theta, P, x0, y0, r0, miou);
	/*ART迭代*/
	Iterator(r, theta, P, X0, XX);

	/************释放空间*****************************/
	free(r);
	free(theta);
	for (i = 0; i < n; i++)
	{
		free(P[i]);
	}
	free(P);
	for (i = 0; i < m; i++)
	{
		free(X0[i]);
		free(XX[i]);
	}
	free(X0);
	free(XX);

	time(&end);
	cost = difftime(end, start);
	printf("%f\n", cost);
}