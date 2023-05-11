#include "Core.h"

int local::count = 0;
int coords::count = 0;
double u = 1e50;

double C[8][8] = {
	{8},
	{4,8},
	{4,2,8},
	{2,4,4,8},
	{4,2,2,1,8},
	{2,4,1,2,4,8},
	{2,1,4,2,4,2,8},
	{1,2,2,4,2,4,4,8}
};
double Bx[8][8] = {
	{4},
	{-4,4},
	{2,-2,4},
	{-2,2,-4,4},
	{2,-2,1,-1,4,4},
	{-2,2,-1,1,-4,4},
	{1,-1,2,-2,2,-2,4},
	{-1,1,-2,2,-2,2,-4,4}
};
double By[8][8] = {
	{4},
	{2,4},
	{-4,-2,4},
	{-2,-4,2,4},
	{2,1,-2,-1,4},
	{1,2,-1,-2,2,4},
	{-2,-1,2,1,-4,-2,4},
	{-1,-2,1,2,-2,-4,2,4}
};
double Bz[8][8] = {
	{4},
	{2,4},
	{2,1,4},
	{1,2,2,4},
	{-4,-2,-2,-1,4},
	{-2,-4,-1,-2,2,4},
	{-2,-1,-4,-2,2,1,4},
	{-1,-2,-2,-4,1,2,2,4}
};
double RightV[8];
double C1[4][4] = {
	{4},
	{2,4},
	{2,1,4},
	{1,2,2,4}
};
double RfV[8][8] = {
	{-4,4,-2,2,-2,2,-1,1},
	{-4,4,-2,2,-2,2,-1,1},
	{-2,2,-4,4,-1,1,-2,2},
	{-2,2,-4,4,-1,1,-2,2}
};

void C1Vector(double* vector, double* result)
{
	for (int i = 0; i < 4; i++)
	{
		result[i] = 0;

		for (int j = 0; j < 4; j++)
		{
			int i1, j1;
			if (i > j)
			{
				i1 = i;
				j1 = j;
			}
			else
			{
				i1 = j;
				j1 = i;
			}
			result[i] += C1[i1][j1] * vector[j];
		}
	}
}
double Func(double* x)
{
	return x[0] + x[1] + x[2];
}
double Scal(double* vector1, double* vector2, int n)
{
	double f = 0;
	for (int i = 0; i < n; i++)
		f += vector1[i] * vector2[i];
	return f;
}
double Norm(double* vector, int n)
{
	return sqrt(Scal(vector, vector, n));
}
void Sum(double* vector1, double* vector2, double* result, int n)
{
	for (int i = 0; i < n; i++)
		result[i] = vector1[i] + vector2[i];
}
void Mul(double a, double* vector, double* result, int n)
{
	for (int i = 0; i < n; i++)
		result[i] = a * vector[i];
}
void Mul2(double* a, double* vector, double* result, int n)
{
	for (int i = 0; i < n; i++)
		result[i] = a[i] * vector[i];
}
void Mov(double* vector, double* result, int n)
{
	for (int i = 0; i < n; i++)
		result[i] = vector[i];
}
void SolutionX(MATRIX* s, double* f, double* x)
{
	int i;
	for (i = 0; i < s->n; i++)
		x[i] = f[i];

	//SY=F 
	for (i = 0; i < s->n; i++)
	{
		for (int j = s->iG[i] - 1; j < s->iG[i + 1] - 1; j++)
			x[i] -= s->GG[j] * x[s->jG[j] - 1];
		x[i] /= s->Di[i];
	}

	//STX=Y
	for (i = s->n - 1; i >= 0; i--)
	{
		x[i] /= s->Di[i];
		for (int j = s->iG[i + 1] - 2; j >= s->iG[i] - 1; j--)
			x[s->jG[j] - 1] -= s->GG[j] * x[i];
	}
}

FEM::FEM(string LocalNum, string Coords, string bounds1, string bounds2, string bounds3)
{
	CoordsIn(Coords);
	LocalIn(LocalNum);
	AllocateMem();
	MakeProfile();
	MakeMatrix();
	//Boundary2(bounds2);
	//Boundary3(bounds3);
	Boundary1(bounds1);
	MSG();
	Output();
}
void FEM::CVector(double* vector, double* result)
{
	for (int i = 0; i < 8; i++)
	{
		result[i] = 0;

		for (int j = 0; j < 8; j++)
		{
			int i1, j1;
			if (i > j)
			{
				i1 = i;
				j1 = j;
			}
			else
			{
				i1 = j;
				j1 = i;
			}
			result[i] += C[i1][j1] * vector[j];
		}
	}
}
void FEM::MVector(double* vector, double* result)
{
	int i;
	for (i = 0; i < M.n; i++)
		result[i] = M.Di[i] * vector[i];

	for (i = 0; i < M.n; i++)
	{
		for (int j = M.iG[i]; j < M.iG[i + 1]; j++)
		{
			result[i] += M.GG[j - 1] * vector[M.jG[j - 1] - 1];
			result[M.jG[j - 1] - 1] += M.GG[j - 1] * vector[i];
		}
	}
}
void FEM::LocalIn(string LocalNum)
{
	ifstream in(LocalNum);
	in >> local::count;

	Matrix = new local[local::count];

	for (int i = 0; i < local::count; i++)
	{
		for (int j = 0; j < 8; j++)
			in >> Matrix[i].mas[j];
		in >> Matrix[i].lambda >> Matrix[i].gamma;
	}
	in.close();
}
void FEM::CoordsIn(string Coords)
{
	ifstream in(Coords);
	in >> coords::count;

	M.n = coords::count;

	Set = new coords[M.n];

	for (int i = 0; i < M.n; i++) {
		for (int j = 0; j < 3; j++)
			in >> Set[i].mas[j];
	}


	in.close();
}
void FEM::MakeProfile()
{
	int** prof = new int* [M.n - 1];

	for (int i = 0; i < M.n - 1; i++)
		prof[i] = new int[i + 1];

	for (int i = 0; i < M.n - 1; i++)
		for (int j = 0; j <= i; j++)
			prof[i][j] = 0;

	int s = 0;//Кол-во ненулевых эл-тов
	for (int k = 0; k < local::count; k++)
	{
		for (int i = 0; i < 8; i++)
		{
			for (int j = i + 1; j < 8; j++)
			{
				int i1 = Matrix[k].mas[i];
				int j1 = Matrix[k].mas[j], k_b;

				if (i1 < j1)
				{
					k_b = i1;
					i1 = j1;
					j1 = k_b;
				}

				if (prof[i1 - 2][j1 - 1] == 0)
				{
					s++;
					prof[i1 - 2][j1 - 1] = 1;
				}
			}
		}
	}

	//Формирование массива iG и jG	
	M.jG = new int[s];
	M.GG = new double[s];
	M.iG[0] = 1;
	M.iG[1] = 1;
	for (int i = 0, d = 0; i < M.n - 1; i++) {
		int k = 0;
		for (int j = 0; j <= i; j++)
			if (prof[i][j] == 1) {
				M.jG[d] = j + 1;
				d++;
				k++;
			}
		M.iG[i + 2] = M.iG[i + 1] + k;
	}
}
void FEM::Add(int i, int j, double x)
{
	int k;
	for (k = M.iG[i] - 1; k < M.iG[i + 1] - 1; k++)
		if (M.jG[k] == j + 1)
			break;
	M.GG[k] += x;
}
void FEM::MakeMatrix()
{
	double h[3];

	for (int i = 0; i < M.n; i++)
	{
		M.Di[i] = 0.0;
		f[i] = 0.0;
	}

	for (int i = 0; i < M.iG[M.n] - 1; i++)
		M.GG[i] = 0.0;

	for (int k = 0; k < local::count; k++)
	{
		//Вычисление шага
		//h----------------------------------------------------------
		for (int i = 0; i < 3; i++) //i-по x,y,z  
		{
			int flag = 1;
			int j;
			for (j = 1; j < 8 && flag; j++) //1 узел фиксируем,пробегаем по остальным 
			{
				flag = 0;
				for (int l = 0; l < 3 && !flag; l++)//проверяем,лежат ли точки на нужном ребре
					if (i != l && Set[Matrix[k].mas[0] - 1].mas[l] != Set[Matrix[k].mas[j] - 1].mas[l])
						flag = 1;
			}
			if (!flag)
				h[i] = fabs(Set[Matrix[k].mas[0] - 1].mas[i] - Set[Matrix[k].mas[j - 1] - 1].mas[i]);
		}
		//------------------------------------------------------------
		//формирование элементов матрицы
		//заполнение мссива GG
		double b_k = Matrix[k].lambda * h[0] * h[1] * h[2] / 36.;
		double c_k = h[0] * h[1] * h[2] / 216.;

		double fr[8];
		//вектор правой части для локал. матрицы
		for (int i = 0; i < 8; i++)
			RightV[i] = Func(Set[Matrix[k].mas[i] - 1].mas);
		CVector(RightV, fr);

		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < i; j++)
			{
				int i1 = Matrix[k].mas[i] - 1;
				int j1 = Matrix[k].mas[j] - 1;
				double s = Matrix[k].gamma * c_k * C[i][j] + b_k / (h[0] * h[0]) * Bx[i][j] +
					b_k / (h[1] * h[1]) * By[i][j] + b_k / (h[2] * h[2]) * Bz[i][j];

				//добавка в gg
				if (i1 < j1)
					Add(j1, i1, s);
				else
					Add(i1, j1, s);
			}
			//добавка в диагональ
			M.Di[Matrix[k].mas[i] - 1] += c_k * Matrix[k].gamma * C[i][i] + b_k / (h[0] * h[0]) * Bx[i][i] +
				b_k / (h[1] * h[1]) * By[i][i] + b_k / (h[2] * h[2]) * Bz[i][i];
			//добавка в правую часть
			f[Matrix[k].mas[i] - 1] += c_k * fr[i];
		}
	}
}
void FEM::AllocateMem()
{
	M.Di = new double[M.n];
	f = new double[M.n];
	x = new double[M.n];
	M.iG = new int[M.n + 1];
}
void FEM::SST(MATRIX* my)
{
	int i, j, k, l;

	my->n = M.n;
	my->Di = new double[M.n];
	my->iG = new int[M.n + 1];

	for (i = 0; i <= M.n; i++)
		my->iG[i] = M.iG[i];


	my->jG = new int[M.iG[M.n] - 1];
	my->GG = new double[M.iG[M.n] - 1];

	for (i = 0; i < M.iG[M.n] - 1; i++)
		my->jG[i] = M.jG[i];

	for (j = 0; j < M.n; j++)
	{
		double Sum = 0;
		for (k = M.iG[j] - 1; k < M.iG[j + 1] - 1; k++)
			Sum += my->GG[k] * my->GG[k];

		my->Di[j] = sqrt(fabs(M.Di[j] - Sum));

		for (i = j + 1; i < M.n; i++)
		{
			int number;
			int flag = 1;
			for (l = M.iG[i] - 1; l < M.iG[i + 1] - 1 && flag; l++)
				if (M.jG[l] == j + 1)
					flag = 0;

			number = l - 1;
			if (flag) continue;

			Sum = 0;
			for (k = M.iG[i] - 1; k < M.iG[i + 1] - 1 && M.jG[k] <= j; k++)
			{
				flag = 1;
				for (l = M.iG[j] - 1; l < M.iG[i + 1] - 1 && flag && M.jG[l] <= M.jG[k]; l++)
					if (M.jG[l] == M.jG[k]) flag = 0;

				l--;

				if (!flag)
					Sum += my->GG[l] * my->GG[k];
			}
			my->GG[number] = (M.GG[number] - Sum) / my->Di[j];
		}
	}
}
void FEM::MSG()
{
	double e = 1e-14;
	int maxiter = 10000;
	MATRIX s;
	SST(&s);

	double* r = new double[s.n];
	double* z = new double[s.n];
	double* temp = new double[s.n];
	double* temp2 = new double[s.n];
	double* temp3 = new double[s.n];

	for (int i = 0; i < M.n; i++)
		x[i] = 0;

	MVector(x, temp);
	Mul(-1, temp, temp, s.n);
	Sum(f, temp, r, s.n);
	SolutionX(&s, r, z);

	double e2 = Norm(r, s.n) / Norm(f, s.n);
	double a, b;
	int k;

	for (k = 1; k<maxiter && e2>e; k++)
	{
		//a
		SolutionX(&s, r, temp2);
		MVector(z, temp);
		a = Scal(temp2, r, s.n) / Scal(temp, z, s.n);

		//r
		Mov(r, temp3, s.n);
		Mul(-a, temp, temp, s.n);
		Sum(r, temp, r, s.n);

		//x
		Mul(a, z, temp, s.n);
		Sum(x, temp, x, s.n);

		//b
		b = 1 / Scal(temp2, temp3, s.n);
		SolutionX(&s, r, temp2);
		b *= Scal(temp2, r, s.n);

		//z
		Mul(b, z, z, s.n);
		Sum(z, temp2, z, s.n);

		e2 = Norm(r, s.n) / Norm(f, s.n);
	}
}
void FEM::Boundary2(string bounds2)
{
	int cr[4];//граничные узлы
	double teta[4];//значения потока в узлах
	ifstream in(bounds2);

	while (!in.eof())
	{
		for (int i = 0; i < 4; i++)
			in >> cr[i];

		for (int i = 0; i < 4; i++)
			in >> teta[i];

		double hx = 0.0,
			hy = 0.0,
			hxyz = 1.0,
			eps = 1e-12;

		for (int i = 0; i < 3; i++)
		{
			double tmp = fabs(Set[cr[0] - 1].mas[i] - Set[cr[1] - 1].mas[i]);

			if (tmp > eps)
				hx += tmp;

			tmp = fabs(Set[cr[0] - 1].mas[i] - Set[cr[2] - 1].mas[i]);

			if (tmp > eps)
				hy += tmp;
		}

		hxyz = hx * hy;

		double vec_tet[4];
		C1Vector(teta, vec_tet);

		for (int i = 0; i < 4; i++)
			f[cr[i] - 1] += hxyz / 36.0 * vec_tet[i];
	}
	in.close();
}
void FEM::Boundary3(string bounds3)
{
	int cr[4];//граничные узлы
	double teta[4];//значения потока в узлах
	double betta;
	ifstream in(bounds3);

	while (!in.eof())
	{

		for (int i = 0; i < 4; i++)
			in >> cr[i];

		for (int i = 0; i < 4; i++)
			in >> teta[i];

		in >> betta;

		double hx = 0.0,
			hy = 0.0,
			hxyz = 1.0,
			eps = 1e-12;

		for (int i = 0; i < 3; i++)
		{
			double tmp = fabs(Set[cr[0] - 1].mas[i] - Set[cr[1] - 1].mas[i]);

			if (tmp > eps)
				hx += tmp;

			tmp = fabs(Set[cr[0] - 1].mas[i] - Set[cr[2] - 1].mas[i]);

			if (tmp > eps)
				hy += tmp;
		}
		hxyz = hx * hy;

		for (int i = 0; i < 4; i++)
		{
			M.Di[cr[i] - 1] += betta * hxyz / 36.0 * C1[i][i];
			for (int j = 0; j < i; j++)
			{
				int ind = M.iG[cr[i] - 1] - 1;

				while (M.jG[ind] != cr[j])
					ind++;

				M.GG[ind] += betta * hxyz / 36.0 * C1[i][j];
			}
		}

		double vec_tet[4];
		C1Vector(teta, vec_tet);

		for (int i = 0; i < 4; i++)
			f[cr[i] - 1] += hxyz / 36.0 * vec_tet[i];
	}
	in.close();
}
void FEM::Boundary1(string bounds1)
{
	int i;
	double q;
	int flag;
	ifstream in(bounds1);

	while (in >> i)
	{
		in >> q;
		i--;
		M.Di[i] = u;
		f[i] = q * u;
	}
	in.close();
}
void FEM::Output()
{
	for (int i = 0; i < M.n; i++)
		cout << fixed << setprecision(14) << x[i] << endl;
}


