# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# include "Radon.h"
# include "ART.h"

int m, n; //����r�ͽǶ�theta����ɢ������
double rmax; //ɨ�跶Χ

void main() 
{

	time_t start, end;
	double cost;
	time(&start);
	
	double *r, *theta;//����ָ��r��theta�洢��ɢ�����r��theta
	double **P;//����ָ������P�洢Radon�任�Ľ��,����P��ÿ��Ԫ�ر�ʾ(theta,r)
	double **X0, **XX;//�洢��ʼֵ�͸���ֵ
	double x0, y0, r0, miou;/*x0��y0��СԲ��Բ�����꣬r0��СԲ�뾶��miou��СԲ��������*/
	int i, j;//ѭ������

	/**************������֪��*****************************/
	printf("������СԲ��Բ������Ͱ뾶��Ԫ��֮���ÿո�ֿ���\n");
	scanf_s("%lf%lf%lf", &x0, &y0, &r0);
	printf("������СԲ�������ʣ�\n");
	scanf_s("%lf", &miou);
	printf("��ֱ��������r�ͽǶ�theta����ɢ��������Ԫ��֮���ÿո�ֿ���\n");
	scanf_s("%d%d",&m,&n);

	/************���ٴ洢�ռ�*****************************/
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


	//����ɨ�跶Χ
	rmax = sqrt(x0 * x0 + y0 * y0) + 3 * r0;
	//��ʼ��X0
	for (i = 0; i < m; i++) 
	{
		for (j = 0; j < m; j++) 
		{
			X0[i][j] = 0;
			XX[i][j] = 0;
		}
	}
	/*****************************************************/
	/*��ɢ��r��theta*/
	ComputeR_Theta(r, theta);
	/*��theta�����*/
	randperm(theta, n);
	/*Radon�任*/
	RadonMethod(r, theta, P, x0, y0, r0, miou);
	/*ART����*/
	Iterator(r, theta, P, X0, XX);

	/************�ͷſռ�*****************************/
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