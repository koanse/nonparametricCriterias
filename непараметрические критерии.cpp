#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define REPEATS 100.0

// Варианты критической области
#define DOUBLE_SIDED 1
#define ONE_SIDED_L 2
#define ONE_SIDED_R 3

float C(int n, int m);									// Число сочетаний из n по m
float Criterion1(float *x, int nx, float *y, int ny);	// Критерий рандомизации компонент
														// Фишера для независимых выборок

float Criterion2(float *x, int nx, float *y, int ny);	// Критерий ранговой рандомизации
														// Ансари-Брэдли
int main()
{
	float x[10] = { 1.2, 1.5, 2.4, 1.5, 1.6, 1.5, 2.1, 1.7, 1.4, 1.2 };
	float y[12] = { 2.1, 2.4, 2.1, 1.5, 1.6, 1.4, 1.4, 1.5, 1.9, 1.8, 1.8, 1.3 };
	
	printf("%.8f", Criterion1(x, 10, y, 12));
	getch();
	
	return 0;
}

float Criterion1(float *x, int nx, float *y, int ny)
{
	int i, j, k, n = nx + ny, m;
	float comb, mins, s1 = 0, s2 = 0, sum, N = 0;
	bool *z;

	// Определение выборки с минимальной суммой
	for(i = 0; i < nx; i++)
		s1 += x[i];

	for(j = 0; j < ny; j++)
		s2 += y[j];

	if(s1 < s2)
	{
		mins = s1;
		m = nx;
	}
	else
	{
		mins = s2;
		m = ny;
	}

	z = new bool [n];
	comb = C(n, m);
		
	for(i = 0; i < comb; i++)
	{
		// Случайный выбор элементов из объединенной выборки
		for(j = 0; j < n; j++)
			z[j] = false;

		sum = 0;
		for(j = 0; j < m; j++)
		{
			do k = rand() % n;
			while(z[k] == true);
			
			z[k] = true;

			if(k < nx) sum += x[k];
			else sum += y[k - nx];
		}

		// Подсчет числа благоприятных исходов
		if(sum > mins) N++;
		if(sum == mins) N += 0.5;
	}

	return (float) N / comb;
}

float C(int n, int m)
{
	int i;
	double res = 1;
	
	for(i = 0; i < m; i++)
	{
		res *= n - i;
		res /= m - i;
	}

    return res;
}
	
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////



void W(float *x, int n, float *y, int m, float alpha, int critical_area);	// Критерий Вилкоксона
void Q(float *x, int m, float *y, int n, float alpha, float critical_value);// Критерий Крамера-Уэлча


void W(float *x, int m, float *y, int n, float alpha, int critical_area)
{
	int *r;			// Массив с текущими рангами xi
	float *z;		// Объединенный массив
	int sum, min_sum, max_sum;
	int i, j;
	float p;		// Вероятность отклонения нулевой гипотезы
	int *ns;		// Массив с частотами появления сумм рангов
	
	for(min_sum = 0, i = 1; i <= m; i++)
		min_sum += i;

	for(max_sum = 0, i = n + 1; i <= m + n; i++)
		max_sum += i;

	r = new int[m];
	ns = new int[max_sum - min_sum + 1];
	z = new float[m + n];
	
	for(i = 0; i < m; i++)
		r[i] = i + 1;

	for(i = 0; i <= max_sum - min_sum; i++)
		ns[i] = 0;
	
	// Построение закона распределения для суммы рангов
	i = m - 1;
	while(1)
	{
		for(sum = 0, j = 0; j < m; j++)
			sum += r[j];

		ns[sum - min_sum] += 1;

		if(r[m - 1] == m + n)
		{
			i = m - 2;
			while(r[i + 1] - r[i] < 2 && i != -1)
				i--;
			if(i == -1) break;

			r[i]++;

			for(j = i + 1; j < m; j++)
				r[j] = r[j - 1] + 1;
		}
		else
			r[m - 1]++;
	}

	// Общее количество комбинаций
	for(i = 0; i <= max_sum - min_sum; i++)
		j += ns[i];
		
	// Вычисление доверительного интервала
	p = 0;
	
	if(critical_area == DOUBLE_SIDED)
	{
		for(i = 0; p < alpha; i++)
			p += (float) (ns[i] + ns[max_sum - min_sum - i]) / j;
		i--;
		printf("Interval: [%d; %d]\n", min_sum + i, max_sum - i);
	}

	if(critical_area == ONE_SIDED_L)
	{
		for(i = 0; p < alpha; i++)
			p += (float) ns[max_sum - min_sum - i] / j;
		i--;
		printf("Interval: [-infinity; %d]\n", max_sum - i);
	}
	
	if(critical_area == ONE_SIDED_R)
	{
		for(i = 0; p < alpha; i++)
			p += (float) ns[i] / j;
		i--;
		printf("Interval: [%d; +infinity]\n", min_sum + i);
	}

	printf("M(S) = %.1f\n", (float) m * (m + n + 1) / 2);


	// Нахождение наблюдаемого значения критерия
	srand((unsigned)time(NULL));

	
	sum = 0;
	for(int k = 0; k < REPEATS; k++)
	{
        for(i = 0; i < m; i++)
			z[i] = x[i];
	
		for(i = 0; i < n; i++)
			z[i + m] = y[i];
	
		for(i = 0; i < m; i++)
			r[i] = i + 1;

		for(i = 0; i < m + n; i++)
			for(j = 0; j < m + n - 1; j++)
				if(z[j] > z[j + 1] || (z[j] == z[j + 1] && rand()%2 == 0))
				{
					int k1, k2;
				
					for(k1 = 0; k1 < m; k1++)
						if(r[k1] == j + 1) break;
				
					for(k2 = 0; k2 < m; k2++)
						if(r[k2] == j + 2) break;
				
					if(k1 < m) r[k1]++;
					if(k2 < m) r[k2]--;
					
					p = z[j];
					z[j] = z[j + 1];
					z[j + 1] = z[j];
				}
					
		for(i = 0; i < m; i++)
			sum += r[i];
	}
	
	printf("S = %f", sum / REPEATS);
	
	getch();

	return;
}

void Q(float *x, int m, float *y, int n, float alpha, float critical_value)
{
	float z_av = 0, s = 0, Q;
	int i;

	for(i = 0; i < m; i++)
		z_av += x[i] - y[i];

	z_av /= m;

	for(i = 0; i < m; i++)
		s += (x[i] - y[i] - z_av) * (x[i] - y[i] - z_av);

	s /= m - 1;

	s = powf(s, 0.5f);

	Q = powf((float) n, 0.5f) * z_av / s;

	printf("\nQ = %f, Qcr = %f", Q, critical_value);
	getch();
	return;
}



float Criterion2(float *x, int nx, float *y, int ny)
{
	int i, j, k, l, n = nx + ny, *p = new int [nx];	// r - массив позиций x
												// в объединенной выборке
	int numr = (n + 1) / 2;						// количество рангов
	float *z = new float [nx + ny], tmp;		// объединенный массив
	float sum, W;
	for(i = 0; i < nx; i++)
	{
        z[i] = x[i];
		p[i] = i;
	}

	for(i = 0; i < ny; i++)
		z[nx + i] = y[i];
	
	for(i = 0; i < n; i++)
		for(j = 0; j < n - 1; j++)
			if(z[j + 1] < z[j])
			{
				for(k = 0; k < nx; k++)
					if(p[k] == j) break;

				for(l = 0; l < nx; l++)
					if(p[l] == j) break;

				if(k != nx && l != nx)
				{
					p[k]++;
					p[l]--;
				}
				if(k != nx && l == nx) p[k]++;
				if(k == nx && l != nx) p[l]--;
							
				tmp = z[j];
				z[j] = z[j + 1];
				z[j + 1] = tmp;
			}

	sum = 0;
	for(i = 0; i < nx; i++)
		if(numr - p[i] > 0) sum += numr - p[i];
		else sum += p[i] - numr + 1;

	if((n / 2) * 2 == n)
		W = (sum - (n + 2) / 4) * sqrt((double) 48 * nx * (n - 1) / (ny * (n - 2) * (n + 2)));
	else
		W = (sum - (n + 1) * (n + 1) / (4 * n) ) *
			sqrt((double) 48 * nx * (n - 1) / (ny * (n + 1) * (n + 3)));

	printf("W = %f", W);
	return W;
}