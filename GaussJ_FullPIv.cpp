#include <iostream>
#include <math.h>
using namespace std;

/*
    Gauss-Jordan method of solving linear sets of equations with full pivoting.
    
    "n" is the size of the matrix
    "m" is the number of right handside vectors
    "A" is the coefficient matrix that will be transformed into the inverse
    "B" is the matrix containing the m right handside vectors that 
        will be transformed into the solution vectors
*/

void GaussJ(int n, int m, double A[100][100], double B[100][100])
{
    int ipiv, jpiv, is_piv[100] = {0}, indR[100], indC[100];
    double big, auxPiv, line_nr;
    
    for(int p = 0; p < n; p++)
    {
        big = 0;
        for(int i = 0; i < n; i++)
            if(!is_piv[i])
                for(int j = 0; j < n; j++)
                    if(!is_piv[j] && fabs(A[i][j]) >= fabs(big))
                    {
                        ipiv = i;
                        jpiv = j;
                        big = A[i][j];
                    }
        if(ipiv != jpiv)
        {
            for(int i = 0; i < n; i++)
                swap(A[ipiv][i], A[jpiv][i]);

            for(int i = 0; i < m; i++)
                swap(B[ipiv][i], B[jpiv][i]);
        }

        indR[p] = ipiv;
        indC[p] = jpiv;
        is_piv[jpiv] = 1;

        if(!A[jpiv][jpiv])
            throw("system either inconsistent or dependent");

        auxPiv = A[jpiv][jpiv];
        A[jpiv][jpiv] = 1;

        for(int i = 0; i < n; i++)
            A[jpiv][i] /= auxPiv;
        for(int i = 0; i < m; i++)
            B[jpiv][i] /= auxPiv;

        for(int i = 0; i < n; i++)
            if(i != jpiv)
            {
                line_nr = A[i][jpiv];
                A[i][jpiv] = 0;
                for(int j = 0; j < n; j++)
                    A[i][j] -= line_nr * A[jpiv][j];
                for(int j = 0; j < m; j++)
                    B[i][j] -= line_nr * B[jpiv][j];
            }
    }

    for(int i = n-1; i >= 0; i--)
        if(indR[i] != indC[i])
            for(int j = 0; j < n; j++)
                swap(A[j][indR[i]], A[j][indC[i]]);
}

int n, m;
double A[100][100], B[100][100];

int main()
{
    cin>>n>>m;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            cin>>A[i][j];
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            cin>>B[i][j];

    GaussJ(n, m, A, B);
    
    cout<<endl;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout<<A[i][j]<<" ";

        cout<<endl;
    } cout<<endl;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
            cout<<B[i][j]<<" ";
        cout<<endl;
    }
}
