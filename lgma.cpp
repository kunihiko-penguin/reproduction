#include<iostream>
#include<fstream>
#include<random>
#include<algorithm>
#include<vector>
#include<cmath>

using namespace std;

//declare variables
const int L=18;
const int Length=6*L,Width=2*L;
const int Size=Length*Width;
int Flag=0;

const int sweep=10000,relax=1000;

const double SigmaStart=0.25, SigmaEnd=0.346, SigmaStep=0.004;
const int Fill=0.522*Size;
int lat1[Fill];//array to store order of filled points
int latspin[Fill];//array to store spin of filled points,corresponding to the order in lat1

//function prototypes
int move_point(int loc_ini,int spin);//output the location of moved point
int mass_center(int*arr);

int main()
{
    //declare variable
    int loc1=0;//location of the point chosen
    int order1=0;//order in lat1 of the point chosen
    int spin1=0;//spin of the point chosen
    int eta=0;//change of spin
    double double_eta=0;
    double p_move=0;
    int randmove=0;//free difusion
    int loc_move=0;//location of the moved point

    //generate random seed
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(0, Fill - 1);
    uniform_int_distribution<> spin(0, 5);//0 1 2 3 4 5 repressents 0 60 120 180 240 300
    uniform_real_distribution<> prob(0.0, 1.0);


    ofstream resfile;
    resfile.open("C://Users//35737//Desktop//mc//lgmactive//L_18_3.txt");

    if (!resfile) {
        printf("fail to open");
        return 1;}

    ofstream visualfile;
    visualfile.open("C://Users//35737//Desktop//mc//lgmactive//L18visual.txt");


    for (double Sigma=SigmaStart;Sigma<SigmaEnd;Sigma+=SigmaStep)
    {
        Flag+=1;
        resfile << " \n";
		resfile << Sigma << "\n ";

        normal_distribution<double> gaussian(0.0, Sigma);

        //initialize lattice
        vector<int> rdorder(Size);//array to choose order of point to be filled
        int lattice[Size]={0};

        for(int i=0;i<Size;++i)
        {
            rdorder[i]=i;
        }
        shuffle(rdorder.begin(),rdorder.end(),gen);
        for (int i=0;i<Fill;++i)
        {
            lattice[rdorder[i]]=1;
            lat1[i]=rdorder[i];
            latspin[i]=spin(gen);
        }


        for (int step=0;step<sweep+relax;++step)
        {
            for (int k=0;k<5*Size;++k)
            {
                order1=distrib(gen);
                loc1=lat1[order1];

                //change the spin
                double_eta=gaussian(gen);
                eta = static_cast<int>(round(double_eta));
                spin1=latspin[order1];
                spin1+=eta%6;
                spin1=(spin1+6)%6;
                latspin[order1]=spin1;

                p_move=prob(gen);

                if(p_move<1.0/6.0)
                {
                    randmove=spin(gen);
                    loc_move=move_point(loc1,randmove);
                    if (lattice[loc_move]==0)
                    {
                        lattice[loc1]=0;
                        lattice[loc_move]=1;
                        lat1[order1]=loc_move;
                    }

                }
                else
                {
                    loc_move=move_point(loc1,spin1);
                    if (lattice[loc_move]==0)
                    {
                        lattice[loc1]=0;
                        lattice[loc_move]=1;
                        lat1[order1]=loc_move;
                    }
                }

            }

            if (step>=relax)
            {
                //calculate mass centre to decide sampling interval
                //coordinate transformation
                int loc_center=mass_center(lat1);
                //sample inteval
                int row_center=loc_center/Length,col_center=loc_center%Length;
                int row_a[L],col_dense[L],row_b[L],col_thin[L];

                /*if (step == sweep+relax-1) {
                    switch(Flag) {
                        case 1:
                        case 14:
                        case 21:
                            printf("%d %d\n", row_center, col_center);
                            break;
                    }
                }*/
             

                for (int i=0;i<L;++i)
                {
                    row_a[i] = row_center-i < 0 ? row_center-i+Width : row_center-i;
                    col_dense[i] = col_center-i+L/2 < 0 ? col_center-i+Length+L/2 : 
                                col_center-i+L/2 >= Length ? col_center-i-Length+L/2 : col_center-i+L/2;
                    row_b[i] = row_center+i+1 >= Width ? row_center+i+1-Width : row_center+i+1;
                    col_thin[i] = col_center-Length/2+L/2+1-i < 0 ? col_center+Length/2+L/2+1-i :
                                col_center-Length/2+L/2+1-i >= Length ? col_center-3*Length/2+L/2+1-i : col_center-Length/2+L/2+1-i;
                }


                int N_da=0,N_db=0,N_ta=0,N_tb=0;
                //N of dense above
                int loc_sample;
                for(int r:row_a)
                {
                    for(int c:col_dense)
                    {
                        loc_sample=r*Length+c;
                        N_da+=lattice[loc_sample];
                    }
                }
                //N of dense below
                for(int r:row_b)
                {
                    for(int c:col_dense)
                    {
                        loc_sample=r*Length+c;
                        N_db+=lattice[loc_sample];
                    }
                }
                //N of thin above
                for(int r:row_a)
                {
                    for(int c:col_thin)
                    {
                        loc_sample=r*Length+c;
                        N_ta+=lattice[loc_sample];
                    }
                }
                //N of thin below
                for(int r:row_b)
                {
                    for(int c:col_thin)
                    {
                        loc_sample=r*Length+c;
                        N_tb+=lattice[loc_sample];
                    }
                }
                resfile<<N_da<<" "<<N_db<<" "<<N_ta<<" "<<N_tb<<" ";
            } 
        }
        if (Flag==1 ||Flag==16 ||Flag==24)
        {
            for (int i=0;i<Size;++i)
            {
                visualfile<<lattice[i]<<" ";
            }
                
        }

        printf("%lf",Sigma);
    }
    return 0;
}

int move_point(int loc_ini,int spin)
{
    int row=loc_ini/Length;
    int col=loc_ini%Length;
    int row_new=0,col_new=0;
    switch(spin)
    {
        case 0:
            row_new=row;
            col_new=col+1>=Length?col+1-Length:col+1;
            break;
        case 1:
            if(row%2==0)
            {
                col_new=col+1>=Length?col+1-Length:col+1;
            }
            else
            {
                col_new=col;
            }
            row_new=row-1<0?row-1+Width:row-1;
            break;
        case 2:
            row_new=row-1<0?row-1+Width:row-1;
            if (row%2==0)
            {
                col_new=col;
            }
            else
            {
                col_new=col-1<0?col-1+Length:col-1;
            }
            break;
        case 3:
            row_new=row;
            col_new=col-1<0?col-1+Length:col-1;
            break;
        case 4:
            row_new=row+1>=Width?row+1-Width:row+1;
            if(row%2==0)
            {
                col_new=col;
            }
            else
            {
                col_new=col-1<0?col-1+Length:col-1;
            }
            break;
        case 5:
            row_new=row+1>=Width?row+1-Width:row+1;
            if(row%2==0)
            {
                col_new=col+1>=Length?col+1-Length:col+1;
            }
            else
            {
                col_new=col;
            }
            break;
    }
    return row_new*Length+col_new;
}

int mass_center(int*arr)
{
    double r=0,c=0;
    double Radius_l=Length/2/3.1415926,Radius_w=Width/2/3.1415926;
    double theta_l=0,theta_w=0,theta=0,x=0,y=0,X_bar=0,Y_bar=0;
    int int_r=0,int_c=0;

    //center in length direction
    for (int i = 0; i < Fill; i++)
    {
        int_r = arr[i]/Length;
        int_c = arr[i]%Length;
        c = int_c/1.0;

        if(int_r%2==0)
        {
            theta = 2 * 3.1415926 * (c-0.5) / Length;
        }
        else
        {
            theta = 2 * 3.1415926 * c / Length;
        }
        x=Radius_l*cos(theta);
        y=Radius_l*sin(theta);
        X_bar+=x;
        Y_bar+=y;
    }
    theta_l=atan2(-Y_bar,-X_bar)+3.1415926;
    X_bar=0;Y_bar=0;
    //center in width direction
    for (int i = 0; i < Fill; i++)
    {
        int_r = arr[i]/Length;
        r=int_r/1.0;
        theta = 2 * 3.1415926 * r / Width;
        x=Radius_w*cos(theta);
        y=Radius_w*sin(theta);
        X_bar+=x;
        Y_bar+=y;
    }
    theta_w=atan2(-Y_bar,-X_bar)+3.1415926;

    int r1 = static_cast<int>(round(theta_w / 2 / 3.1415926 * Width));
    int c_new=0,c1=0,r_new=0;
    if(r_new%2==0)
    {
        c1 = static_cast<int>(round(theta_l / 2 / 3.1415926 * Length+0.5));
    }
    else
    {
        c1 = static_cast<int>(round(theta_l / 2 / 3.1415926 * Length));
    }

    c_new=c1>=Length?c1-Length:c1;
    r_new=r1>=Width?r1-Width:r1;

    int loc_new = r_new *Length+c_new;
    return loc_new;
}