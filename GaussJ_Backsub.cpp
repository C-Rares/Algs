#include <iostream>
#include <math.h>
using namespace std;
/*
    Gauss-Jordan method of solving linear sets of equations using backsubstitution. The idea is that we might encounter zeros on the main diagonal
despite searching for maximum magnitude pivot (because we only form zeros below and not above also) so the first loop is dedicated towards pushing 
the existing zeros below the main diagonal, and only then we start making matrix operations.

    "n" is the size of the matrix
    "A" is the coefficient matrix 
    "B" is the matrix containing the right handside vector (You can solve for multiple right handsides)
*/

void GaussJ(int n, double A[100][100], double B[100])
{
    int ln;
    double big, auxPiv, line_nr, nr, X[100] = {0};
    
    for(int p = 0; p < n; p++)
    {
        ln = -1;
        if(!A[p][p])
            for(int i = 0; i < n; i++)
                if(i != p)
                    if(A[i][p] && A[p][i])
                    {
                        for(int j = 0; j < n; j++)
                            swap(A[p][j], A[i][j]);
                        
                        swap(B[p], B[i]);
                        break;
                    }
    }

    for(int ln = 0; ln < n; ln++)
        for(int i = ln + 1; i < n; i++)
        {
            nr = A[i][ln];
            if(nr)
            {
                for(int j = 0; j < n; j++)
                    A[i][j] = A[ln][ln] * A[i][j] / nr - A[ln][j];
                B[i] = B[i] * A[ln][ln]/nr - B[ln];
            }
        }


    double sum;
    for(int i = n-1; i >= 0; i--)
    {
        sum = 0;
        for(int j = i+1; j < n; j++)
            sum += A[i][j] * X[j];

        X[i] = 1/A[i][i]*(B[i] - sum);
    }

    for(int i = 0; i < n; i++)
        B[i] = X[i];
}

int n;
double A[100][100], B[100];

int main()
{
    cin>>n;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            cin>>A[i][j];

    for(int i = 0; i < n; i++)
            cin>>B[i];

    GaussJ(n, A, B);
    
    cout<<endl;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout<<A[i][j]<<" ";

        cout<<endl;
    } cout<<endl;

    for(int i = 0; i < n; i++)
        cout<<B[i]<<" ";
}