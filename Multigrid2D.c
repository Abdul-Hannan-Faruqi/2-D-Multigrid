//    Last change:  Abdul Hannan  03 Jan 2019
#include<stdio.h>
#include<math.h>
#include<time.h>

//Global Variables
int prolong(), power(int, int);;
void rest(), add();
double T[300][300], res[300][300], err[300][300], Tcalc, Tn, h, BE, e=1e-10, er, x, y, x2, y2, r, pi=3.14159265359;
int i, j, k, l, w=2, c, m=256, a, init, cyc, iter=1, n=256, f, f2, ln, lmin, level=0;
float lmax;

int main()
{
//Declaration (Local Variables)
float le, b, v;
char bc, sch, choice;
FILE *fp, *fp2, *fp3;

//File opening
fp = fopen("Poisson Log.csv", "a");     //Data logging
if(fp == NULL){
    printf("Couldn't open file\n");
    return 1;
}

fp2 = fopen("Plot.dat", "w");   //For plotting results in Tecplot
if(fp2 == NULL){
    printf("Couldn't open file\n");
    return 1;
}

fp3 = fopen("Residual.csv", "w");   //For recording observations
if(fp3 == NULL){
    printf("Couldn't open file\n");
    return 1;
}

/*Inputs_________________________________________________________________________________________________________________*/
x = 1.0/m;
x2 = x*x;
y = 1.0/n;
y2 = y*y;
printf("x = %f", x);
printf("\ty = %f", y);

x: printf("\n\nChoose an AMG Scheme\nEnter a character: [V]\t[W]\n");
scanf("%c", &sch);
if (sch == 'V')
{
    fprintf(fp, "V,");
}
if (sch=='W')
{
    fprintf(fp, "W,");
    printf("Enter the grid level at W point:\t");
    scanf("%d", &lmin);
    fprintf(fp, "%d,", lmin);
}
printf("\nEnter the no. of levels: ");
scanf("%d", &ln);
fprintf(fp, "%d,", ln);
printf("\nEnter the no. of inner iterations: ");
scanf("%d", &init);
printf("\nEnter the no. of %c cycles: ", sch);
scanf("%d", &cyc);
fprintf(fp, "%d,", cyc);
/*_________________________________________________________________________________________________________________________*/

/*Solution Initialization__________________________________________________________________________________________________*/
for (i=0; i<=n; i++)
{
    for (j=0; j<=m; j++)
    {
        T[i][j]=0;
    }
}
/*_________________________________________________________________________________________________________________________*/

//Solution
clock_t begin = clock();
while (1)
{
    BE=0;
/*Gauss-Seidel____________________________________________________________________________________________________________*/
    for (i= 1; i<n; i++)
    {
        for (j= 1; j<m; j++)
        {
            T[i][j]= (x2*(T[i+1][j] + T[i-1][j]) + y2*(T[i][j+1] + T[i][j-1])+ 2*x2*y2*((1-6*x2*j*j)*y2*i*i*(1-y2*i*i) + (1-(6*y2*i*i))*x2*j*j*(1-x2*j*j)))/(2*(x2+y2));
        }
    }
    
    fprintf(fp3, "%g,", T[n/2][m/2]);

    for (i= 1; i<n; i++)
    {
        for (j= 1; j<m; j++)
        {
            res[i][j] =-2*((1-6*x2*j*j)*y2*i*i*(1-y2*i*i) + (1-6*y2*i*i)*x2*j*j*(1-x2*j*j)) - (T[(i+1)][j] + T[(i-1)][j]- 2*T[i][j])/y2 - (T[i][(j+1)] + T[i][(j-1)]- 2*T[i][j])/x2;
        }
    }
/*_________________________________________________________________________________________________________________________*/

/*Multigrid________________________________________________________________________________________________________________*/
    iter=1;
    if ((sch == 'V') || (sch=='v'))
    {
        while(iter<=cyc)
        {
            k=n/w;
            l=m/w;
            level++;
            rest();
            k=k*w;
            l=l*w;
            while (level>0)
            {
                prolong ();
            }
            add();
            iter++;
        }
    }

    else if ((sch == 'W') || (sch == 'w'))
    {
        while(iter<=cyc)
        {
            k=n/2;
            l=m/2;
            level++;
            rest();
            k=k*2;
            l=l*2;
        while (level>lmin)
        {
            prolong ();
        }
        k=k/2;
        l=l/2;
        rest();
        k=k*2;
        l=l*2;
        while (level>0)
        {
            prolong ();
        }
        add();
        iter++;
        }
    }
/*_________________________________________________________________________________________________________________________*/

    fprintf(fp3, "%g\n", BE);

    if (BE<e)
        break;
}

clock_t end= clock();
double time=(double)(end-begin)/CLOCKS_PER_SEC;

/*Printing the Result______________________________________________________________________________________________________*/
printf("\n\nNo. of Iterations = %d\nResidual = %g\nRun time = %g s\n", iter, BE, time);
fprintf(fp, "%d,%g\n", init, time);

//Writing Tecplot .dat file
fprintf(fp2, "zone\ni=%d, j=%d\n", n+1, m+1);
for (i=0; i<=n; i++)
{
    for (j=0; j<=m; j++)
    {
    Tcalc= x2*y2*i*i*j*j*(1-j*j*x2)*(i*i*y2-1);
    er = Tcalc - T[i][j];
    printf("%g  ", T[i][j]);
    fprintf(fp2, "%f  %f  %g\n", x*j, y*i, T[i][j]);      //Tecplot File
    }
    printf("\n");
}
fclose (fp2);
fclose (fp);
fclose(fp3);
/*_________________________________________________________________________________________________________________________*/
printf("\n\nCLOCKS_PER_SEC = %g\t", CLOCKS_PER_SEC);
printf("\n\n\nRun again with same conditions: [Y]\t[N]\n");
scanf("%*c%c", &choice);
if (choice=='Y')
{
    level++;
    goto x;
}

return 0;
}
/*_________________________________________________________________________________________________________________________*/

void add()
{
    for (i=1; i<n; i++)
    {
        for (j=1; j<m; j++)
        {
            Tn = T[i][j] + err[i][j];
            T[i][j] = Tn;
            if (fabs(res[i][j])>BE)
                BE = fabs(res[i][j]);
        }
    }
}

void rest()
{
    while (level<=ln)
    {
        f = power(w ,level);
        f2= f*f;
        for (a=1; a<init; a++)
        {
            for (i= 1; i<k; i++)
            {
                for (j= 1; j<l; j++)
                {
                    res[f*i][f*j] = (res[f*(i-1)][f*j]+res[f*(i+1)][f*j]+res[f*i][f*(j+1)]+res[f*i][f*(j-1)])/4;
                    err[f*i][f*j] = (f2*y2*(err[f*(i+1)][f*j]+err[f*(i-1)][f*j]) + f2*x2*(err[f*i][f*(j+1)]+err[f*i][f*(j-1)])-f2*f2*x2*y2*res[f*i][f*j])/(2*f2*(x2+y2));
                }
            }
        }
        for (i= 1; i<n; i++)
        {
            for (j= 1; j<m; j++)
            {
                err[i][j]= (err[(i+1)][j]+err[(i-1)][j] + err[i][(j+1)]+err[i][(j-1)])/4;
            }
        }
        k = k/w;
        l=l/w;
        level++;
    }
    lmax=level;
}

int prolong()
{
    level--;
    f = power(w,level);
    for (i= 1; i<k; i++)
    {
        for (j= 1; j<l; j++)
        {
            err[f*i][f*j]= (4*err[f*i][f*j]+err[f*(i+1)][f*j]+err[f*(i-1)][f*j] + err[f*i][f*(j+1)]+err[f*i][f*(j-1)])/8;
        }
    }
    k=w*k;
    l=w*l;
    return level;
}

int power(int w, int level)
{
    int i, f=1;
    for (i=1; i<=level; i++)
        f=f*w;
    i++;
    return f;
}
