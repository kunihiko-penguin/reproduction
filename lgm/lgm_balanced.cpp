#include <iostream>
#include <random>
#include <array>
#include <fstream>
#include <algorithm>
using namespace std;

// declare variables
const int L = 6;
const int Length = 6 * L, Width = 2 * L;
const int Size = Length * Width;
const int HalfSize = Size / 2;
const int sweep = 10000, relax = 10000;
int Flag = 0;

// generate random seed
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> distrib(0, HalfSize - 1);
uniform_real_distribution<> prob(0.0, 1.0);

// function prototype
double energy_loc(int *arr, int);
int mass_center(int *arr);

int main()
{
    ofstream resfile;
    resfile.open("C://Users//35737//Desktop//mc//lgm//L6c.txt");

    if (!resfile)
    {
        cout << "fail to open" << endl;
        return 1;
    }   

    ofstream visualfile;
    visualfile.open("C://Users//35737//Desktop//mc//lgm//visualization//L16.txt");

    for(int a=0;a<1;++a)
    {
    

        // declare variable
        int lattice[Length * Width];
        double dE;                // deltaE of changing empty and non-empty grid
        array<int, Size> rdorder; // array to choose order of point to be filled

        // genarate initial configuration
        for (double T = 0.45; T < 1; T += 1) // 0.45 0.725 0.01  0.564 0.577 0.003
        {
            Flag += 1;
            resfile << " \n";
            resfile << T << "\n ";

            // array recording location of 0(empty) and 1(non-empty)
            int lat0[HalfSize];
            int lat1[HalfSize];

            // lattice initialization

            for (int i = 0; i < Size; ++i)
            {
                rdorder[i] = i;
            }
            shuffle(rdorder.begin(), rdorder.end(), gen);

            for (int i = 0; i < HalfSize; ++i)
            {
                int ord1 = rdorder[i]; // location in lattice of non-empty
                lat1[i] = ord1;
                lattice[ord1] = 1;
                int ord0 = rdorder[i + HalfSize]; // location in lattice of empty
                lat0[i] = ord0;
                lattice[ord0] = 0;
            }
            // to this line we generate a lattice of rho=0.5

            // evolve the system
            for (int step = 0; step < sweep + relax; ++step)
            {
                int loc1 = 0, loc0 = 0;
                int order1, order0;
                double rdnum;
                double p_flip;

                for (int k = 0; k < 5 * Size; ++k)
                {
                    order1 = distrib(gen);
                    order0 = distrib(gen);
                    loc1 = lat1[order1];
                    loc0 = lat0[order0];
                    dE = energy_loc(lattice, loc0) - energy_loc(lattice, loc1);
                    rdnum = prob(gen);
                    p_flip = exp(-dE / T);

                    // cout<<"dE="<<dE<<" rdnum="<<rdnum<<" p_flip="<<p_flip<<endl;//debug!!!!!!

                    // decide whether to exchange
                    if (p_flip >= rdnum)
                    {
                        lat1[order1] = loc0;
                        lat0[order0] = loc1;
                        lattice[loc1] = 0;
                        lattice[loc0] = 1;

                        // cout<<"exchange"<<endl;//debug!!!!!!
                    }
                }
                // calculate statistics needed(rho-dense;rho-thin;mass center;total energy;)
                double E_tot = 0;

                // calculate mass centre to decide sampling interval
                // coordinate transformation
                int loc_center = mass_center(lat1);
                // sample inteval
                int row_center = loc_center / Length, col_center = loc_center % Length;
                int row_a[L] = {}, col_dense[L] = {}, row_b[L] = {}, col_thin[L] = {};

                // cout<<"row_center="<<row_center<<" col_center="<<col_center<<endl;//debug!!!!!

                for (int i = 0; i < L; ++i)
                {
                    row_a[i] = row_center - i < 0 ? row_center - i + Width : row_center - i;
                    col_dense[i] = col_center - i + L / 2 < 0 ? col_center - i + Length + L / 2 : col_center - i + L / 2 >= Length ? col_center - i - Length + L / 2
                                                                                                                                : col_center - i + L / 2;
                    row_b[i] = row_center + i + 1 >= Width ? row_center + i + 1 - Width : row_center + i + 1;
                    col_thin[i] = col_center - Length / 2 + L / 2 + 1 - i < 0 ? col_center + Length / 2 + L / 2 + 1 - i : col_center - Length / 2 + L / 2 + 1 - i >= Length ? col_center - 3 * Length / 2 + L / 2 + 1 - i
                                                                                                                                                                            : col_center - Length / 2 + L / 2 + 1 - i;
                }

                // cout<<"above="<<row_a[L-1]<<" below="<<row_b[L-1]<<" dense_left="<<col_dense[L-1]<<" dense_right="
                //<<col_dense[0]<<" thin_left="<<col_thin[L-1]<<" thin_right="<<col_thin[0]<<endl;//debug!!!!!

                if (step >= relax)
                {
                    int N_da = 0, N_db = 0, N_ta = 0, N_tb = 0;
                    // N of dense above
                    int loc_sample;
                    for (int r : row_a)
                    {
                        for (int c : col_dense)
                        {
                            loc_sample = r * Length + c;
                            N_da += lattice[loc_sample];
                        }
                    }
                    // N of dense below
                    for (int r : row_b)
                    {
                        for (int c : col_dense)
                        {
                            loc_sample = r * Length + c;
                            N_db += lattice[loc_sample];
                        }
                    }
                    // N of thin above
                    for (int r : row_a)
                    {
                        for (int c : col_thin)
                        {
                            loc_sample = r * Length + c;
                            N_ta += lattice[loc_sample];
                        }
                    }
                    // N of thin below
                    for (int r : row_b)
                    {
                        for (int c : col_thin)
                        {
                            loc_sample = r * Length + c;
                            N_tb += lattice[loc_sample];
                        }
                    }
                    resfile << N_da << " " << N_db << " " << N_ta << " " << N_tb << " ";

                    // total free energy
                    for (int i : lat1)
                    {
                        E_tot += energy_loc(lattice, i);
                    }
                    E_tot *= 0.5;
                    resfile << E_tot << " ";
                }
            }

            if (Flag == 1 || Flag == 14 || Flag == 21)
            {
                for (int i = 0; i < Size; ++i)
                {
                    visualfile << lattice[i] << " ";
                }

            }

            cout << T << endl;
        }
    }
    resfile.close();
    return 0;
}

// define functions
double energy_loc(int *arr, int loc)
{
    int r = loc / Length, c = loc % Length;
    int left, right, above, below;
    left = c == 0 ? Length - 1 : c - 1;
    right = c == Length - 1 ? 0 : c + 1;
    above = r == 0 ? Width - 1 : r - 1;
    below = r == Width - 1 ? 0 : r + 1;
    return -arr[loc] * (arr[r * Length + left] + arr[r * Length + right] + arr[above * Length + c] + arr[below * Length + c]);
}

int mass_center(int *arr)
{
    int Radius_l = Length / 2 / 3.1415926, Radius_w = Width / 2 / 3.1415926;
    double theta_l = 0, theta_w = 0, theta = 0, x = 0, y = 0, X_bar = 0, Y_bar = 0;

    // center in length direction
    for (int i = 0; i < HalfSize; i++)
    {
        int c = arr[i] % Length;
        theta = 2 * 3.1415926 * c / Length;
        x = Radius_l * cos(theta);
        y = Radius_l * sin(theta);
        X_bar += x;
        Y_bar += y;
    }
    theta_l = atan2(-Y_bar, -X_bar) + 3.1415926;
    X_bar = 0;
    Y_bar = 0;
    // center in width direction
    for (int i = 0; i < HalfSize; i++)
    {
        int r = arr[i] / Length;
        theta = 2 * 3.1415926 * r / Width;
        x = Radius_w * cos(theta);
        y = Radius_w * sin(theta);
        X_bar += x;
        Y_bar += y;
    }
    theta_w = atan2(-Y_bar, -X_bar) + 3.1415926;
    int c_loc = static_cast<int>(theta_l / 2 / 3.1415926 * Length + 0.5);
    int r_loc = static_cast<int>(theta_w / 2 / 3.1415926 * Width + 0.5);

    int loc = c_loc % Length + r_loc % Width * Length;
    return loc;
}