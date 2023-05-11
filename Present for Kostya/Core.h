#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <functional>
#include <iostream>
#include <iomanip>

using namespace std;

struct local
{
	static int count;
	int mas[8];
	double lambda;
	double gamma;
};
struct coords
{
	double mas[3];
	static int count;
};
struct MATRIX
{
	int n;
	double* Di;
	double* GG;
	int* iG;
	int* jG;
};
class FEM
{
public:
   FEM(string LocalNum, string Coords, string bounds1, string bounds2, string bounds3);

private:
	local* Matrix;
	coords* Set;
	MATRIX M;

	double* f;
	double* x;

	void MakeProfile();
	void MakeMatrix();
	void Add(int i,int j,double x);
   void LocalIn(string LocalNum);
	void AllocateMem();
	void CoordsIn(string Coords);
	void MSG();
	void SST(MATRIX* my);
	void Output();
	void Boundary2(string bounds2);
	void Boundary3(string bounds3);
	void Boundary1(string bounds1);
	void CVector(double* vector, double* result);
	void MVector(double* vector, double* result);

	// Õ¿’”ﬂ ›“Œ????????
	double GetPsi(int num_fe, int num_basis, double x, double y, double z);
	double U(double nx, double y, double z);
};