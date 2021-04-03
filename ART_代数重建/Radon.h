#pragma once
# include <stdio.h>
# include <math.h>
# define pi 3.1415926 //定义圆周率

extern int m;
extern int n;
extern double rmax;

/*将theta随机化存储*/
int randperm(double * theta,int n) {
	
	for (int i = 0; i < n; i++)
	{
		int a = int(rand()%n);
		double temp = theta[i];
		theta[i] = theta[a];
		theta[a] = temp;
	}
	return 0;
}


/*离散化r和theta*/
void ComputeR_Theta(double *r, double *theta) 
{
	int i;//循环变量
	double dr;//定义变量dr存储探测器单元格的长度
	
	dr = 2 * rmax / m;
	//计算向量theta
	for (i = 0; i < n; i++) 
	{
		theta[i] = 0 + i * pi / n;
	}
	//计算向量r
	for (i = 0; i < m; i++) 
	{
		r[i] = -rmax + 0.5 * dr + i * dr;
	}
}
/*Radon变换*/
void RadonMethod(double *r, double *theta, double **P, double x0, double y0, double r0, double miou)
{
	FILE *fp;//定义文件指针
	int i, j;//循环变量
	double distance; //存储圆心到射线的距离
	double RL = 0; //定义变量RL表示在给定的(theta,r)下的Radon变换结果
	for (i = 0; i < n; i++) 
	{
		for (j = 0; j < m; j++) 
		{
			distance = fabs(x0 * cos(theta[i]) + y0 * sin(theta[i]) - r[j]);
			//判断射线是否穿过小圆
			if (distance < r0) //穿过小圆
			{
				RL = miou * 2*sqrt(r0 * r0 - distance * distance);
			}
			else
			{
				RL = 0;
			}
			P[i][j] = RL;
		}
	}

	//将Radon变换的数据写入文件当中
	fp=fopen("D://Matlab/workspace/CT-2019-7/sector.txt","w");
	if (!fp)
	{
   	    printf("can not open the file!\n");
	    exit(-1);
	}
	for (i=0;i<n;i++)
	{
	    for (j=0;j<m;j++)
	   {
	        fprintf(fp,"%9.6lf",P[i][j]);
	   }
	   fprintf(fp,"\n");
	}
	fclose(fp);
	printf("Success!\n");
}