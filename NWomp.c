#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>
int max(int num1, int num2);
char *substring(char *string, int position, int length);
char *charStrConcat(char c, char *str);

char *refSeq, *alignSeq;

int M[200][200];
#define GAP -4;
#define MATCH 5;
#define MISSMATCH -3;
#define MAX(x, y) (((x) > (y)) ? (x) : (y))


void generate(int c) {
    //Generates the values of a
    long long int i;
    for (i = 0; i < c; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            refSeq[i] = 'C';
        else if (aux == 2)
            refSeq[i] = 'A';
        else if (aux == 3)
            refSeq[i] = 'T';
        else
            refSeq[i] = 'G';
    }
    //Generates the values of b
    for (i = 0; i < c; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            alignSeq[i] = 'A';
        else if (aux == 2)
            alignSeq[i] = 'C';
        else if (aux == 3)
            alignSeq[i] = 'G';
        else
            alignSeq[i] = 'T';
    }
}

void ScoreTable(char seq1[],char seq2 [],int num_proc){
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);
    int i,j;
    #pragma omp parallel num_threads(num_proc) \
    default(none) shared(seq1len,seq2len,seq1,seq2,M) private(i,j)
    for (int i = 1; i < seq1len + 1; i++)
    {
        #pragma omp for
        for (int j = 1; j < seq2len + 1; j++)
        {
            int scoreDiag = 0;

            if (seq1[j - 1] == seq2[i - 1]){
                scoreDiag = M[i - 1][j - 1] + MATCH;
            }
            else{
                scoreDiag = M[i - 1][j - 1] + MISSMATCH;
            }

            int scoreLeft = M[i][j - 1] + GAP;
            int scoreUp =  M[i - 1][j] + GAP;

            int maxScore = MAX(MAX(scoreDiag, scoreLeft), scoreUp);
            M[i][j] = maxScore;
        }
    }
}

//Step 3: PrintTable function
void PrintTable(char seq1[],char seq2 []){
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);

    for (int i = 0; i < seq1len + 1; i++){

        printf("\n");
        #pragma parallel for
        for (int j = 0; j < seq2len + 1; j++){
            printf("%d\t",M[i][j]);
        }

    }

}

int main(){
int numOfChars;
int num_proc=0;

printf("Enter Number of Threads:");
scanf("%d",&num_proc);
printf("\nEnter Number of Character in sequence:");
scanf("%d",&numOfChars);
  //allocate memory for sequences
   refSeq = malloc(numOfChars * sizeof(char));
    alignSeq = malloc(numOfChars * sizeof(char));

  //open file and get the sequences
  //refSeq=["A","T","C","T","G","A","T","C","C","A"];
  //alignSeq=["T","G","A","A","C","C","C","A","T","T"];
  generate(numOfChars);



  int refSeqCnt = strlen(refSeq) + 1;
  int alignSeqCnt = strlen(alignSeq) + 1;
  int scoringMatrix [alignSeqCnt][refSeqCnt];

  int scoreLeft, scoreUp;
  double elapsTime;


  //get start time

  double start = omp_get_wtime();
#pragma omp parallel sections
{

#pragma omp section
{


   ScoreTable(refSeq,alignSeq,num_proc);
   PrintTable(refSeq,alignSeq);
}
  //initialization

  #pragma omp section
{

  for (int i = 0; i < alignSeqCnt; i++){
    for (int j = 0; j < refSeqCnt; j++){
      scoringMatrix[i][j] = 0;
    }
  }


  //fill the matrix
  for (int i = 1; i < alignSeqCnt; i++){
    for (int j = 1; j < refSeqCnt; j++){
      int scoreDiag = 0;
      char *sub1,*sub2;
      char c1 = refSeq[j-1];
      char c2 = alignSeq[i-1];
      if (c1 == c2){
        scoreDiag = scoringMatrix[i - 1][j - 1] + 2;
      } else {
        scoreDiag = scoringMatrix[i - 1][j - 1] + -1;
      }

      scoreLeft = scoringMatrix[i][j - 1] - 2;
      scoreUp = scoringMatrix[i - 1][j] - 2;

      int maxScore = max(max(scoreDiag, scoreLeft), scoreUp);
      scoringMatrix[i][j] = maxScore;
    }
  }
}
}
  //traceback
  char *AlignmentA="";
  char *AlignmentB=" ";

  int m = alignSeqCnt - 1;
  int n = refSeqCnt - 1;

  while (m > 0 || n > 0){
    int scoreDiag = 0;
    if (m == 0 && n > 0){
      AlignmentA = charStrConcat(refSeq[n-1],AlignmentA);
      AlignmentB = charStrConcat('-',AlignmentB);
      n = n - 1;
    } else if (n == 0 && m > 0){
      AlignmentA = charStrConcat('-',AlignmentA);
      AlignmentB = charStrConcat(alignSeq[m-1],AlignmentB);
      m = m - 1;
    } else {
      if (alignSeq[m - 1] == refSeq[n - 1]){
        scoreDiag = 2;
      } else {
        scoreDiag = -1;
      }

      if (m > 0 && n > 0 && scoringMatrix[m][n] == scoringMatrix[m - 1][n - 1] + scoreDiag){
        AlignmentA = charStrConcat(refSeq[n-1],AlignmentA);
        AlignmentB = charStrConcat(alignSeq[m-1],AlignmentB);
        m = m - 1;
        n = n - 1;
      } else if (n > 0 && scoringMatrix[m][n] == scoringMatrix[m][n - 1] - 2) {
        AlignmentA = charStrConcat(refSeq[n-1],AlignmentA);
        AlignmentB = charStrConcat('-',AlignmentB);
        n = n - 1;
      } else {
        AlignmentA = charStrConcat('-',AlignmentA);
        AlignmentB = charStrConcat(alignSeq[m-1],AlignmentB);
        m = m - 1;
      }
    }
  }


  //get end time
  printf("\n");
  printf("\n");

  printf("%s\n", AlignmentA);
  printf("%s\n", AlignmentB);


double end = omp_get_wtime();


  //calculate total time
  elapsTime = (double)(end - start) / CLOCKS_PER_SEC;
  printf("\n\nTotal Time(in Milliseconds) = %lf\n\n",elapsTime * 1000);
}

int max(int num1, int num2) {
  return num1 > num2 ? num1 : num2;
}

char *substring(char *string, int position, int length){
  char *pointer;
  int c;
  pointer = malloc(length+1);

  if (pointer == NULL){
    printf("Unable to allocate memory.\n");
    exit(1);
  }

  for (c = 0 ; c < length ; c++){
    *(pointer+c) = *(string+position-1);
    string++;
  }

  *(pointer+c) = '\0';
  return pointer;
}

char *charStrConcat(char c, char *str){
  char *tempStr;
  char charStr[2];
  tempStr = (char*) malloc(sizeof(char *)* (strlen(str)));
  strcpy(tempStr,str);
  str = (char *) malloc(sizeof(char *) * (strlen(str)+1));
  charStr[0] = c;
  charStr[1] = '\0';
  strcpy(str, charStr);
  strcat(str,tempStr);
  return str;
}

