#include <bits/stdc++.h>
#include <math.h>
using namespace std;

void GaussJ(int n, int m, double A[100][100], double B[100][100])
{
    int ipiv, jpiv, is_piv[100] = {0}, indR[100], indC[100];
    double big, auxPiv, line_nr;
    
    for(int p = 0; p < n; p++)
    {
        /* There are always n pivots that we need to find. We look throught the matrix A for the biggest element in magnitude to be our pivot
           for computational reasons */

        big = 0; // Setting variable for a biggest number elligible to be vector
        for(int i = 0; i < n; i++)    // Start searching for pivots on each column
            if(!is_piv[i])    // If we have a pivot on column i (position (i,i)) than we are not interested because the column is computed with 0's (in our case, the inverse elements, we'll see..)  
                for(int j = 0; j < n; j++)    // Else we search the column for a biggest number and store the position
                    if(!is_piv[j])  /* We can t have a pivot* on the row that we have found our next possible pivot because there can be only 1 pivot* per row and by interchanging the rows
                                     we are destroying our past computations that we have done to establish that row's already pivot*  */
                        if(fabs(A[i][j]) >= fabs(big))     //  If there's no pivot we ask if it's bigger than our other findings
                        {               // And then we store the information about it
                            ipiv = i;
                            jpiv = j;
                            big = A[i][j];
                        }   // And so on..
                    
        // Now we have chosen our pivot located on line [ipiv] and column [jpiv]
        // We swap the rows in order to have our pivot on the position [jpiv],[jpiv] on the main diagonal (all pivots are located on the main diagonal)
        if(ipiv != jpiv)   // That is, of course, if it isn't already positioned on the main diagonal
        {
            for(int i = 0; i < n; i++)
                swap(A[ipiv][i], A[jpiv][i]);

            for(int i = 0; i < m; i++)        // Don't forget to swap the lines of the right handside vectors also, so that our solutions don't get lost
                swap(B[ipiv][i], B[jpiv][i]);
        }

        /*  For normal G-J elimination, we would have a seperate identity matrix that would suffer the reduced row echelon operations to become the inverse of A
           And the original matrix A would eventually turn into the identity matrix.
            After we have chosen a pivot* and we used row operations to make 0's above and below our pivot* we move to the next pivot**.
           We realise that no further operations made to compute the pivot** column will ever change the value of our pivot* column.
           (That is because we would substract the pivot* by a number times 0, which is 0. And the 0-elements the same way, that is because we would never substract
           by the pivot* in the future because he is the only pivot on that line) 
            In conclusion, we can use the eternally unchanged computed column to store our inverse of A! That is, after we've made our pivot 1(by dividing the row)
           so that we know the final solution.
            But, by interchanging rows (in the process of choosing the pivot) we also change the rows of what would be the parallel identity matrix so we need not
           to forget that in the end we have to change them back
        */

        indR[p] = ipiv;     // We mark the position of the first biggest pivot by row and column, so we can unscramble the inverse in the end if the rows have been changed
        indC[p] = jpiv;
        is_piv[jpiv] = 1; 

        // If the biggest number that we've found was a 0 it means we have a 0 as pivot, which is nonsense (concluded from translation into equations system)
        if(!A[jpiv][jpiv])
            throw("inconsistent system, no possible solution");

        auxPiv = A[jpiv][jpiv]; // Storing the pivot 
        A[jpiv][jpiv] = 1;    // We make the pivot position 1 so it mimmicks the position of the identity matrix (which has 1's on the main diagonal)

        for(int i = 0; i < n; i++)
            A[jpiv][i] /= auxPiv;   // And we divide the row by the pivot
        for(int i = 0; i < m; i++)
            B[jpiv][i] /= auxPiv;   // And of course, the solutions vectors

        // Now we do the classic row operations of making 0's above and below our pivot
        for(int i = 0; i < n; i++)
            if(i != jpiv) // only below and above the pivot!
            {   //The operations we are doing is { Ri <- Ri - line_nr * Rp }
                line_nr = A[i][jpiv]; // This is the element on our pivot's column
                
                A[i][jpiv] = 0; // We make 0's below and above beforehand because so it mimmicks the identity matrix's poisitions
                                // and we know the row will be left with only 1 non zero number so the solution is untouched

                for(int j = 0; j < n; j++)  // And then we do the operations that leaves our matrix A equivalent to the what would be the transformed identity matrix
                    A[i][j] -= line_nr * A[jpiv][j];
                for(int j = 0; j < m; j++)  // And obviously, the same operations on the solution vectors
                    B[i][j] -= line_nr * B[jpiv][j];
            }
    }

    /* We have reached the final step. We have our solutions in the B vectors but our inverse is scrambled from the row operations 
       so we swap the columns of the rows that we have interchanged, in reverse order */
    for(int i = n-1; i >= 0; i--)
        if(indR[i] != indC[i])
            for(int j = 0; j < n; j++)
                swap(A[j][indR[i]], A[j][indC[i]]);
}

int main()
{
    int n, m;
    double A[100][100], B[100][100];
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

