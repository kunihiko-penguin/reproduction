#include<iostream>
#include<cmath>
#include<fstream>
#include<random>
#include<vector>
#include<array>
using namespace std;

const int  sweeps = 500000, relax = 30000;
// a = 1.5, b = 1;
const int length = 32, width = 32;


const double H = 0;//external field

//function prototype
double energy_loc(int* arr, int);

int main()
{
	//declare local variables
	int s = length * width;
	int lattice[length * width];

	//generate random number
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distrib(0, s - 1);
	uniform_int_distribution<>distrib_spin(0, 1);
	uniform_real_distribution<> prob(0.0, 1.0);

	ofstream resfile;
	resfile.open("C://Users//35737//Desktop//mc//static//L32//c.txt");

	if (!resfile) {
		cout << "fail to open" << endl;
		return 1;
	}

	for (int a=0;a<9;++a)
	{
		for (double T = 2.2670; T <2.2715 ; T = T + 0.001)
		{
			//some statistical parameters
			double E_tot, M_tot;
			resfile << " \n";
			resfile << T << "\n ";

			//generate random initial configuration
			for (int i = 0; i < s; ++i)
			{
				int rdnum = distrib_spin(gen);
				lattice[i] = rdnum * 2 - 1;

				//visual
				/*for (int i = 0; i < L * L; ++i)
				{
					if (i % L == L-1) {
						cout << lattice[i] << endl;
					}
					else {
						cout << lattice[i] << ' ';
					}
				}*/


			}


			//sweeps
			for (int sweep = 0; sweep < sweeps + relax; ++sweep)
			{
				//evolve the system
				int loc = 0;
				for (int k = 0; k < 5 * s; ++k)
				{
					//choose a random spin
					loc = distrib(gen);

					//cout << loc << ' ' << lattice[loc] << endl;

					/*int r = loc / L, c = loc % L;
					int left, right, above, below;
					left = c == 0 ? L - 1 : c - 1; //c changed
					right = c == L - 1 ? 0 : c + 1;//new c
					above = r == 0 ? L - 1 : r - 1;//r
					below = r == L - 1 ? 0 : r + 1;//r
					cout << lattice[r * L + left] << ' ' << lattice[r * L + right] << ' ' <<lattice[above * L + c] << ' ' <<lattice[below * L + c]<<endl;*/

					//energy chage of operation
					double dE = -2 * energy_loc(lattice, loc);

					//decide whether to flip
					double flip_prob = exp(-dE / T);
					double p = prob(gen);

					//cout << flip_prob << ' ' << p << endl;

					if (flip_prob > p)
					{
						lattice[loc] *= -1;

					}


					//visual
					/*for (int i = 0; i < L * L; ++i)
					{
						if (i % L == L-1) {
							cout << lattice[i] << endl;
						}
						else {
							cout << lattice[i] << ' ';
						}
					}*/

				}


				if (sweep < relax)
				{
					continue;
				}

				else
				{
					//E and M for initial config
					E_tot = 0, M_tot = 0;
					for (int i = 0; i < s; ++i)
					{
						E_tot += energy_loc(lattice, i);
						M_tot += lattice[i];
						//record it!!!

					}
					E_tot *= 0.5;
					resfile << E_tot << " " << M_tot << " ";
				}
			}
			cout << T << endl;
		}
	}
	resfile.close();
	return 0;
}

double energy_loc(int* arr, int loc)
{
	int r = loc / length, c = loc % length;
	int left, right, above, below;
	double E_loc = 0;
	//position r*L+c
	//boundary condition
	left = c == 0 ? length - 1 : c - 1; //c changed
	right = c == length - 1 ? 0 : c + 1;//new c
	above = r == 0 ? width - 1 : r - 1;//r
	below = r == width - 1 ? 0 : r + 1;//r

	E_loc = -arr[loc] * (arr[r * length + left] + arr[r * length + right] + arr[above * length + c] + arr[below * length + c]);
	return E_loc;
}