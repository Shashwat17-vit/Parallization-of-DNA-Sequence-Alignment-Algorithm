// Parallel Smith Waterman Algorithm For Large Sequences

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <omp.h>
#include <time.h>

//Constants
#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3

void SetColor(int ForgC)
 {
     WORD wColor;

      HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
      CONSOLE_SCREEN_BUFFER_INFO csbi;

                       //We use csbi for the wAttributes word.
     if(GetConsoleScreenBufferInfo(hStdOut, &csbi))
     {
                 //Mask out all but the background attribute, and add in the forgournd     color
          wColor = (csbi.wAttributes & 0xF0) + (ForgC & 0x0F);
          SetConsoleTextAttribute(hStdOut, wColor);
     }
     return;
 }


void similarityScore(long long int i, long long int j, int* H, int* P, long long int* maxPos);
int matchMissmatchScore(long long int i, long long int j);
void backtrack(int* P, long long int maxPos);
void printMatrix(int* matrix);
void printPredecessorMatrix(int* matrix);
void generate(void);
long long int nElement(long long int i);
void calcFirstDiagElement(long long int *i, long long int *si, long long int *sj);

/* End of prototypes */

 int m ; //Columns - Size of string a
 int n ;  //Lines - Size of string b

// Define Scoring System
int matchScore = 5;
int missmatchScore = -3;
int gapScore = -4;

// SEQUENCE Array
char *a, *b;

int main() {
    int thread_count;
    int ans;
    printf("Enter the Number of Threads:");
    scanf("%d",&thread_count);
    printf("\nGenerate Random Sequence(1/0):");
    scanf("%d",&ans);
    if(ans==1)
    {printf("\nEnter the Length of Sequence A:");
    scanf("%d",&m);
    printf("\nEnter the Length of Sequence B:");
    scanf("%d",&n);
    a = malloc(m * sizeof(char));
    b = malloc(n * sizeof(char));
    generate();
    }
    else if(ans==0)
    {   printf("Please Calculate Sequence Size and Enter Lengths:");
        scanf("%d",&m);
        scanf("%d",&n);
        a = malloc(m * sizeof(char));
        b = malloc(n * sizeof(char));
        printf("\nEnter Sequence A:");
        for(int i=0;i<m;i++)
        scanf("%s",&a[i]);
        printf("\nEnter Sequence B:");
        for(int j=0;j<n;j++)
        scanf("%s",&b[j]);
    }


    void generate() {
    //Generates the values of a
    long long int i;
    for (i = 0; i < m; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            a[i] = 'A';
        else if (aux == 2)
            a[i] = 'C';
        else if (aux == 3)
            a[i] = 'G';
        else
            a[i] = 'T';
    }
    //Generates the values of b
    for (i = 0; i < n; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            b[i] = 'A';
        else if (aux == 2)
            b[i] = 'C';
        else if (aux == 3)
            b[i] = 'G';
        else
            b[i] = 'T';
    }
}

    //The Scoring Matrix or Similarity Matrix contains "0" valued First row and Columns
    m++;
    n++;

    //Allocates similarity matrix H
    int *H;
    H = calloc(m * n, sizeof(int));

    //Allocates predecessor matrix P
    int *P;
    P = calloc(m * n, sizeof(int));

    //Start position for backtrack
    long long int maxPos = 0;
    //Calculates the similarity matrix
    long long int i, j;

    //Gets Initial time
    double initialTime = omp_get_wtime();

    long long int si, sj, ai, aj;

    //Because now we have zeros ((m-1) + (n-1) - 1)
    long long int nDiag = m +n -3;
    long long int nEle;

    #pragma omp parallel num_threads(thread_count) \
    default(none) shared(H, P, maxPos, nDiag) private(nEle, i, si, sj, ai, aj,j)
    {
        for (i = 1; i <= nDiag; ++i)
        {
            nEle = nElement(i);
            calcFirstDiagElement(&i, &si, &sj);

            for (j = 1; j <= nEle; ++j)
            {
                ai = si - j + 1;
                aj = sj + j - 1;
                similarityScore(ai, aj, H, P, &maxPos);
            }
        }
    }

    backtrack(P, maxPos);


    printf("\nSimilarity Matrix:\n");
    printMatrix(H);
       double finalTime = omp_get_wtime();
    printf("\nElapsed time: %f\n\n", finalTime - initialTime);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P);


// FREE the Variables and Array
    free(H);
    free(P);
    free(a);
    free(b);

    return 0;
}
 // Calculate the number of i-diagonal elements

long long int nElement(long long int i) {
    if (i < m && i < n) {
        //Number of elements in the diagonal is increasing
        return i;
    }
    else if (i < max(m, n)) {
        //Number of elements in the diagonal is stable
        long int min = min(m, n);
        return min - 1;
    }
    else {
        //Number of elements in the diagonal is decreasing
        long int min = min(m, n);
        return 2 * min - i + abs(m - n) - 2;
    }
}

//     Calculate the position of (si, sj)-element

void calcFirstDiagElement(long long int *i, long long int *si, long long int *sj) {
    // Calculate the first element of diagonal
    if (*i < n) {
        *si = *i;
        *sj = 1;
    } else {
        *si = n - 1;
        *sj = *i - n + 2;
    }
}

 //  Calculate and create the Similarity Matrix H And Predecessor Matrix P.

void similarityScore(long long int i, long long int j, int* H, int* P, long long int* maxPos) {

    int up, left, diag;


    long long int index = m * i + j;

    //Get element above
    up = H[index - m] + gapScore;

    //Get element on the left
    left = H[index - 1] + gapScore;

    //Get element on the diagonal
    diag = H[index - m - 1] + matchMissmatchScore(i, j);

    //Calculates the maximum
    int max = NONE;
    int pred = NONE;
    if (diag > max) { //same letter ↖
        max = diag;
        pred = DIAGONAL;
    }

    if (up > max) { //remove letter ↑
        max = up;
        pred = UP;
    }

    if (left > max) { //insert letter ←
        max = left;
        pred = LEFT;
    }
    //Inserts the value in the similarity and predecessor matrixes
    H[index] = max;
    P[index] = pred;

    //Updates maximum score to be used as seed on backtrack
    if (max > H[*maxPos]) {
        #pragma omp critical
        *maxPos = index;
    }

}

//Similarity function on the alphabet for match/missmatch

int matchMissmatchScore(long long int i, long long int j) {
    if (a[j - 1] == b[i - 1])
        return matchScore;
    else
        return missmatchScore;
}

// Function to modify matrix to print, path change from value to PATH - Backtrack from maxpostion to 0;

void backtrack(int* P, long long int maxPos) {
    //hold maxPos value
    long long int predPos;

    //backtrack from maxPos to startPos = 0
    do {
        if (P[maxPos] == DIAGONAL)
            predPos = maxPos - m - 1;
        else if (P[maxPos] == UP)
            predPos = maxPos - m;
        else if (P[maxPos] == LEFT)
            predPos = maxPos - 1;
        P[maxPos] *= PATH;
        maxPos = predPos;
    } while (P[maxPos] != NONE);
}
// print Similarity matrix
void printMatrix(int* matrix) {
    long long int i, j;
    printf("-\t-\t");
    for (j = 0; j < m-1; j++) {
    	printf("%c\t", a[j]);
    }
    printf("\n-\t");
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
        	if (j==0 && i>0) printf("%c\t", b[i-1]);
            printf("%d\t", matrix[m * i + j]);
        }
        printf("\n");
    }

}
//   Print predecessor matrix

void printPredecessorMatrix(int* matrix) {
    long long int i, j, index;
    printf("\t\t");
    for (j = 0; j < m-1; j++) {
    	printf("%c\t", a[j]);
    }
    printf("\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
        	if (j==0 && i>0) printf("%c\t", b[i-1]);
            index = m * i + j;

            if (matrix[index] < 0) {
                    SetColor(3);
                if (matrix[index] == -UP)
                    printf("^\t");
                else if (matrix[index] == -LEFT)
                    printf("<\t");
                else if (matrix[index] == -DIAGONAL)
                    printf("[]\t");
                else
                    printf("-\t");
                    SetColor(15);

            } else {
                if (matrix[index] == UP)
                    printf("^\t");
                else if (matrix[index] == LEFT)
                    printf("<\t");
                else if (matrix[index] == DIAGONAL)
                    printf("[]\t");
                else
                    printf("-\t");
            }
        }
        printf("\n");
    }

}
 //Random generate arrays a and b
void generate() {

    //Generates the values of a
    long long int i;
    for (i = 0; i < m; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            a[i] = 'A';
        else if (aux == 2)
            a[i] = 'C';
        else if (aux == 3)
            a[i] = 'G';
        else
            a[i] = 'T';
    }

    //Generates the values of b
    for (i = 0; i < n; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            b[i] = 'A';
        else if (aux == 2)
            b[i] = 'C';
        else if (aux == 3)
            b[i] = 'G';
        else
            b[i] = 'T';
    }
}
