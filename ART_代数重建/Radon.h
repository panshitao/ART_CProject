#pragma once
# include <stdio.h>
# include <math.h>
# define pi 3.1415926 //����Բ����

extern int m;
extern int n;
extern double rmax;

/*��theta������洢*/
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


/*��ɢ��r��theta*/
void ComputeR_Theta(double *r, double *theta) 
{
	int i;//ѭ������
	double dr;//�������dr�洢̽������Ԫ��ĳ���
	
	dr = 2 * rmax / m;
	//��������theta
	for (i = 0; i < n; i++) 
	{
		theta[i] = 0 + i * pi / n;
	}
	//��������r
	for (i = 0; i < m; i++) 
	{
		r[i] = -rmax + 0.5 * dr + i * dr;
	}
}
/*Radon�任*/
void RadonMethod(double *r, double *theta, double **P, double x0, double y0, double r0, double miou)
{
	FILE *fp;//�����ļ�ָ��
	int i, j;//ѭ������
	double distance; //�洢Բ�ĵ����ߵľ���
	double RL = 0; //�������RL��ʾ�ڸ�����(theta,r)�µ�Radon�任���
	for (i = 0; i < n; i++) 
	{
		for (j = 0; j < m; j++) 
		{
			distance = fabs(x0 * cos(theta[i]) + y0 * sin(theta[i]) - r[j]);
			//�ж������Ƿ񴩹�СԲ
			if (distance < r0) //����СԲ
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

	//��Radon�任������д���ļ�����
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