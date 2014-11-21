/***********************

Conway Game of Life

OpenMp version

************************/

#include <stdio.h>

#include <stdlib.h>

#include <omp.h>


#define ROWS 100    

#define COLS 100

#define NSTEPS 100    


typedef struct m{
    int rows;
    int cols;
    int* mat;
}matrix;


int getPos(int i, int j)
{
    int pos = i * COLS + j;
    return pos;
}

void update(matrix *old,matrix *new,int thread_count)
{
    # pragma omp parallel for num_threads(thread_count) shared(old,new)
    for(int i=0;i<old->rows;i++)
    {
        for(int j=0;j<old->cols;j++)
        {
            old->mat[getPos(i,j)] = new->mat[getPos(i,j)];      
        }
    }
}

void initMatrix(matrix *m, int thread_count)
{
    int rows,cols;
    float x;
    printf("Enter matrix dimensions rows and cols\n");
    scanf("%d %d",&rows,&cols);
    m->rows = rows;
    m->cols = cols;
    m->mat = malloc(rows*cols*sizeof(int));
    # pragma omp parallel for num_threads(thread_count) private(x) shared(m)
    for(int i=1;i<m->rows-1;i++)
    {
        for(int j=1;j<m->cols-1;j++)
        {
            x = rand()/((float)RAND_MAX + 1);
            if(x < 0.5)
                m->mat[getPos(i,j)] = 0;
            else
                m->mat[getPos(i,j)]= 1;
        }
    }

}

void printMatrix(matrix *m, char* filepath)
{
    int i,j;
    FILE* f = fopen(filepath,"w+");
    for(i = 0;i < m->rows;i++)
    {
        for(j=0; j < m->cols;j++)
            fprintf(f," %d ",m->mat[getPos(i,j)]);
        fprintf(f,"\n");
    }
    fclose(f);
}


void fillBorders(matrix *m, int thread_count)
{
    m->mat[0] = m->mat[getPos(m->rows-2,m->cols-2)];
    m->mat[getPos(m->rows-1,m->cols-1)] = m->mat[getPos(1,1)];
    m->mat[getPos(0,m->cols-1)] = m->mat[getPos(m->rows-2,1)];
    m->mat[getPos(m->rows-1,0)] = m->mat[getPos(1,m->cols-2)];

    # pragma omp parallel num_threads(thread_count) shared(m)
    # pragma omp for
    for(int j=1;j<m->cols-1;j++)
    {
        m->mat[getPos(0,j)]=m->mat[getPos(m->rows-2,j)];  
    }
    # pragma omp for
    for(int j=1;j<m->cols-1;j++)
    {      
        m->mat[getPos(m->rows-1,j)] = m->mat[getPos(1,j)];
    }
    # pragma omp for
    for(int i=1;i<m->rows-1;i++)
    {
        m->mat[getPos(i,0)] = m->mat[getPos(i,m->cols-2)];
               
    }
    # pragma omp for
    for(int i=1;i<m->rows-1;i++)
    {     
        m->mat[getPos(i,m->cols-1)] = m->mat[getPos(i,1)];      
    }
}

void calc(matrix *old, matrix *new,int thread_count)
{
    # pragma omp parallel for num_threads(thread_count)
    for(int i=1;i<old->rows-1;i++)
    {
        for(int j=1;j<old->cols-1;j++)
        {
            int nsum=0;
            nsum += old->mat[getPos(i+1,j)];
            nsum += old->mat[getPos(i-1,j)];
            nsum += old->mat[getPos(i,j-1)];
            nsum += old->mat[getPos(i,j+1)];
            nsum += old->mat[getPos(i+1,j+1)];
            nsum += old->mat[getPos(i+1,j-1)];
            nsum += old->mat[getPos(i-1,j-1)];
            nsum += old->mat[getPos(i-1,j+1)];
            printf("nsum = %d\n",nsum);

            switch(nsum)
            {
                case 3:
                    new->mat[getPos(i,j)] = 1;
                    break;
                case 2:
                    new->mat[getPos(i,j)] = old->mat[getPos(i,j)];
                    break;
                default:
                    new->mat[getPos(i,j)] = 0;
                    break;        
            }

        }
    }
}


int main(int argc, char *argv[])
{
    int thread_count = strtol(argv[1],NULL,10);
    int n;
    int global_result=0;
    char* filepath = "/home/luis/CPD/PCP/OpenMP/live/matrix.txt";
    matrix old,new;
    initMatrix(&old,thread_count);
    initMatrix(&new,thread_count);
    

    for(n=0;n<NSTEPS;n++)
    {
        fillBorders(&old,thread_count);
        calc(&old,&new,thread_count);
        printMatrix(&old,filepath);
        update(&old,&new,thread_count);            
    }

    # pragma omp parallel for num_threads(thread_count) reduction(+ : global_result) shared(new)
    for(int i=1;i<new.rows-1;i++)
    {
        for(int j=1;j<new.cols-1;j++)
        {
            global_result += new.mat[getPos(i,j)];
        }
    }
    free(old.mat);
    free(new.mat); 
    printf("\nNumber of live cells = %d\n", global_result);

  return 0;
}