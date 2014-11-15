#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "param.h"
#include "macro.h"

float AlR = 0.2;
float BtR = 0.6;
float MuR = 0.01;
float BtF = 0.6;
float AlF = -1.8;
float MuF = 0.02;
char* s = "Passou";

char* filepathLRabbit = "/home/luis/CPD/PCP/Echo/fileRabbit.txt";
char* filepathLFox = "/home/luis/CPD/PCP/Echo/fileFox.txt";
#define TRUE 1
#define FALSE 0

int getPos(int i, int j)
{
    int pos = i*(WE_Size) + j;
    return pos;
}

void printMatrix(float* matrix,int line_offset,char* filepath,int rank)
{
    int i,j,pos;
    FILE* f = fopen(filepath,"w+");
    fprintf(f, "Matriz local do processo %d\n",rank );
    for (i = 0; i < line_offset; i++)
    {
        fprintf(f, "linha: %d \n", i);
        for (j = 0; j < WE_Size; j++)
        {
            pos=getPos(i,j);
            fprintf(f, " %f ", matrix[pos]);
        }
        fprintf(f, "\n" );
    }
    fclose(f);
}

int SetLand ( float *Rabbit,float *Fox, int offset,int my_rank,int comm_sz);

int nonCriticalLines(float* Rabbit,float* Fox,float* TRabbit,float* TFox,int offset);

int criticalLines(float* lineUpR,float* lineDownR,float* lineUpF,float* lineDownF,
    float* Rabbit,float* Fox,float* TRabbit,float* TFox, int my_rank, int comm_sz,int offset);

int FillBorder(float  *Animal);

int GetPopulation(float *Animal,float *tcount,int offset);





int main(int argc, char *argv[]) 
{
    int rank,comm_size;
    MPI_Comm my_grid;
    int dim[2],period[2],reorder;
    int up,down,right,left,line_offset,err=0,k,j;
    float nbrab,nbfox,totalrab,totalfox;
    MPI_Init(&argc, &argv);
    double MPI_Wtime(void);
    double start, finish;
    start=MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Request reqR[4],reqF[4];
    MPI_Status statusR[4],statusF[4];
    line_offset = (NS_Size)/comm_size;
    dim[0] = 1; 
    dim[1] = comm_size;
    period[0]=TRUE; 
    period[1]=TRUE;
    reorder=FALSE;
    MPI_Cart_create(MPI_COMM_WORLD,2,dim,period,reorder,&my_grid);
    /*
    if(rank==8)
    {
        MPI_Cart_shift(my_grid,0,1,&left,&right);
        MPI_Cart_shift(my_grid,1,1,&up,&down);
        printf("P:%d My neighbors are r: %d d:%d 1:%d u:%d\n",rank,right,down,left,up);
    }
    */
    if(rank == comm_size -1)
    {
        line_offset = line_offset + (int)NS_Size % comm_size;
    }
    //Criação de dados derivados
    MPI_Datatype row;
    MPI_Type_contiguous((WE_Size), MPI_FLOAT, &row);
    MPI_Type_commit(&row);
    MPI_Datatype matrix;
    MPI_Type_contiguous(NS_Size*WE_Size, MPI_FLOAT, &matrix);
    MPI_Type_commit(&matrix);
    MPI_Datatype local_matrix;
    MPI_Type_contiguous(line_offset, row, &local_matrix);
    MPI_Type_commit(&local_matrix);

    float* sendUp_rowR = malloc(WE_Size*sizeof(float));
    float* sendDown_rowR = malloc(WE_Size*sizeof(float));

    float* recvUp_rowR = malloc(WE_Size*sizeof(float));
    float* recvDown_rowR = malloc(WE_Size*sizeof(float));

    float* sendUp_rowF = malloc(WE_Size*sizeof(float));
    float* sendDown_rowF = malloc(WE_Size*sizeof(float));

    float* recvUp_rowF = malloc(WE_Size*sizeof(float));
    float* recvDown_rowF = malloc(WE_Size*sizeof(float));

    float* lmRabbit = malloc(line_offset*WE_Size*sizeof(float));
    float* lmFox = malloc(line_offset*WE_Size*sizeof(float));
    
    float* lmTRabbit = malloc(line_offset*WE_Size*sizeof(float));
    float* lmTFox = malloc(line_offset*WE_Size*sizeof(float));
    
    float *Rabbit,*Fox;

    err = SetLand(lmRabbit,lmFox,line_offset,rank,comm_size);
    if(rank == 0)
    {

        Rabbit = malloc(NS_Size*WE_Size*sizeof(float));
        Fox = malloc(NS_Size*WE_Size*sizeof(float));
    }
    MPI_Barrier(my_grid);

    MPI_Gather(lmRabbit,1,local_matrix,Rabbit,1,matrix,0,my_grid);
    //printf("%s o processo %d\n", s,rank);
    MPI_Gather(lmFox,1,local_matrix,Fox,1,matrix,0,my_grid); 
    //printf("%s o processo %d\n", s,rank);

    for( k=1; k<=NITER; k++) 
    {

        nbrab=0;
        nbfox=0;
        totalrab=0;
        totalfox=0;

        if(rank==0)
        {
            //err = FillBorder(Rabbit); 
            //err = FillBorder(Fox);
        }
        MPI_Scatter(Rabbit,1,local_matrix,lmRabbit,1,local_matrix,0,my_grid);
        MPI_Scatter(Fox,1,local_matrix,lmFox,1,local_matrix,0,my_grid);
        MPI_Cart_shift(my_grid,1,1,&up,&down);

        if (up == -1)
        {
            down = MPI_PROC_NULL;
        }
        if (down == comm_size)
        {
            down = MPI_PROC_NULL;
        }

        //err = nonCriticalLines(lmRabbit,lmFox,lmTRabbit,lmTFox,line_offset);
        int pos;

        //////UP Lines

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(0,j);
            lmRabbit[pos] = sendUp_rowR[j];
        }

        MPI_Isend(sendUp_rowR, 1, row, up, 0, my_grid,reqR);


        MPI_Recv(recvDown_rowR,1, row, down, 0, my_grid, statusR);

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(0,j);
            lmFox[pos] = sendUp_rowF[j];
        }
        MPI_Isend(sendUp_rowF, 1, row, up, 0, my_grid,reqF);

        MPI_Recv(recvDown_rowF,1, row, down, 0, my_grid, statusF);

        //////DOWN Lines

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(line_offset-1,j);
            lmRabbit[pos] = sendDown_rowR[j];
        }
        MPI_Isend(sendDown_rowR, 1, row, down, 0, my_grid,reqR);

        MPI_Recv(recvUp_rowR,1, row, up, 0, my_grid, statusR);

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(line_offset-1,j);
            lmRabbit[pos] = sendDown_rowF[j];
        }
        MPI_Isend(sendDown_rowF, 1, row, down, 0, my_grid,reqF);

        MPI_Recv(recvUp_rowF,1, row, up, 0, my_grid, statusF);


        MPI_Wait(reqR,statusR);
        MPI_Wait(reqF,statusF);

        //err = criticalLines(recvUp_rowR,recvDown_rowR,recvUp_rowF,recvDown_rowF,lmRabbit,lmFox,lmTRabbit,lmTFox,rank,comm_size,line_offset);           

        if( (k % PERIOD) == 1 )
        {
            err = GetPopulation(lmTRabbit,&nbrab,line_offset); 
            err = GetPopulation(lmTFox,&nbfox,line_offset);
            MPI_Reduce(&nbrab,&totalrab,1,MPI_FLOAT,MPI_SUM,0,my_grid);
            MPI_Reduce(&nbfox,&totalfox,1,MPI_FLOAT,MPI_SUM,0,my_grid);
            printf("In the year %d, the number of rabbits is %f and the number of foxes is %f. \n", k,totalrab,totalfox);
        }
        MPI_Barrier(my_grid);
        MPI_Gather(lmRabbit,1,local_matrix,Rabbit,1,matrix,0,my_grid);
        MPI_Gather(lmFox,1,local_matrix,Fox,1,matrix,0,my_grid);
        printf("%s o processo %d\n", s,rank);

    }
    if(rank==0)
    {
        finish=MPI_Wtime();
        printf("tempo decorrido: %f\n", finish - start);
    }
    MPI_Finalize();
    return 0;
}


int SetLand ( float *Rabbit,float *Fox, int offset,int my_rank,int comm_sz)
{
    int err,pos;
    int gi, gj;

    err = 1;
    if(my_rank==0)
    {
        for( gi=1;  gi < offset; gi++) 
        {
            for( gj=1; gj<=WE_Size-2; gj++) 
            {   
                pos=getPos(gi,gj);
                Rabbit[pos] =
                    128.0*(gi-1)*(NS_Size-gi)*(gj-1)*(WE_Size-gj) /
                    (float)(NS_Size*NS_Size*WE_Size*WE_Size);
                Fox[pos] =
                    8.0*(gi/(float)(NS_Size)-0.5)*(gi/(float)(NS_Size)-0.5)+ 
                    8.0*(gj/(float)(WE_Size)-0.5)*(gj/(float)(WE_Size)-0.5);
            }
        }
    }
    if((my_rank+1) == comm_sz)
    {
        for( gi=0;  gi < (offset-1); gi++) 
        {
            for( gj=1; gj<=WE_Size-2; gj++) 
            {
                pos=getPos(gi,gj);
                Rabbit[pos] =
                    128.0*(gi-1)*(NS_Size-gi)*(gj-1)*(WE_Size-gj) /
                    (float)(NS_Size*NS_Size*WE_Size*WE_Size);
                Fox[pos] =
                    8.0*(gi/(float)(NS_Size)-0.5)*(gi/(float)(NS_Size)-0.5)+ 
                    8.0*(gj/(float)(WE_Size)-0.5)*(gj/(float)(WE_Size)-0.5);
            }
        }
    }
    else
    {
        for( gi=0;  gi < offset; gi++) 
        {
            for( gj=1; gj<=WE_Size-2; gj++) 
            {
                pos=getPos(gi,gj);
                Rabbit[pos] =
                    128.0*(gi-1)*(NS_Size-gi)*(gj-1)*(WE_Size-gj) /
                    (float)(NS_Size*NS_Size*WE_Size*WE_Size);
                Fox[pos] =
                    8.0*(gi/(float)(NS_Size)-0.5)*(gi/(float)(NS_Size)-0.5)+ 
                    8.0*(gj/(float)(WE_Size)-0.5)*(gj/(float)(WE_Size)-0.5);
            }
        }
    }        
    return(err);
}


int nonCriticalLines(float* Rabbit,float* Fox,float* TRabbit,float* TFox,int offset)
{
    int gi,gj;
    int err = 0,pos;

    for( gi=1; gi < offset-1; gi++)
    {
        for( gj=1; gj<=WE_Size-3; gj++)
        {
            pos=getPos(gi,gj);
            TRabbit[pos] = (1.0+AlR-4.0*MuR)*Rabbit[pos] +
                                    BtR*Fox[pos] +
                                    MuR*(Rabbit[pos-1]+Rabbit[pos+1]+
                                    Rabbit[getPos(gi-1,gj)]+Rabbit[getPos(gi+1,gj)]);

            TFox[pos] = AlF*Rabbit[pos] +
                                 (1.0+BtF-4.0*MuF)*Fox[pos] +
                                 MuF*(Fox[pos-1]+Fox[pos+1]+
                                 Fox[getPos(gi-1,gj)]+Fox[getPos(gi+1,gj)]);
        }
    }
    return (err);
}

int criticalLines(float* lineUpR,float* lineDownR,float* lineUpF,float* lineDownF,
    float* Rabbit,float* Fox,float* TRabbit,float* TFox, int my_rank, int comm_sz,int offset)
{
    int gj,pos;
    int err = 0;
    if(my_rank==0)
    {

        for( gj=1; gj<=WE_Size-2; gj++)
        {
            pos=getPos(offset-1,gj);
            TRabbit[pos] = (1.0+AlR-4.0*MuR) * Rabbit[pos] +
                                    BtR * Fox[pos] +
                                    MuR*(Rabbit[pos-1] + Rabbit[pos+1] +
                                    Rabbit[getPos(offset-2,gj)] + lineDownR[gj]);                     
            TFox[pos] = AlF * Rabbit[pos] +
                                (1.0+BtF-4.0*MuF) * Fox[pos] +
                                MuF * (Fox[pos-1] + Fox[pos+1] +
                                Fox[getPos(offset-2,gj)] + lineDownF[gj]);
        }
    }

    if ((my_rank + 1) == comm_sz)
    {
        for( gj=1; gj<=WE_Size-2; gj++)
        {
            pos=getPos(0,gj);
            TRabbit[pos] = (1.0+AlR-4.0*MuR)*Rabbit[pos] +
                                    BtR*Fox[pos] +
                                    MuR*(Rabbit[pos-1]+Rabbit[pos+1]+
                                    lineUpR[gj]+Rabbit[getPos(1,gj)]);

            TFox[pos] = AlF*Rabbit[pos] +
                                    (1.0+BtF-4.0*MuF)*Fox[pos] +
                                    MuF*(Fox[pos-1]+Fox[pos+1]+
                                    lineUpF[gj]+Fox[getPos(1,gj)]);
        }
    }
    else
    {
        for( gj=1; gj<=WE_Size-2; gj++)
        {
            pos=getPos(0,gj);
            TRabbit[pos] = (1.0+AlR-4.0*MuR)*Rabbit[pos] +
                                    BtR*Fox[pos] +
                                    MuR*(Rabbit[pos-1]+Rabbit[pos+1]+
                                    lineUpR[gj]+Rabbit[getPos(1,gj)]);

            TFox[pos] = AlF*Rabbit[pos] +
                                (1.0+BtF-4.0*MuF)*Fox[pos] +
                                MuF*(Fox[pos-1]+Fox[pos+1]+
                                lineUpF[gj]+Fox[getPos(1,gj)]);
        }
        for( gj=1; gj<=WE_Size-2; gj++)
        {
            pos=getPos(offset-1,gj);
            TRabbit[pos] = (1.0+AlR-4.0*MuR)*Rabbit[pos] +
                                BtR*Fox[pos] +
                                MuR*(Rabbit[pos-1] + Rabbit[pos+1]
                                + Rabbit[getPos(offset-2,gj)]+lineDownR[gj]);

            TFox[pos] = AlF*Rabbit[pos] +
                                (1.0+BtF-4.0*MuF)*Fox[pos] +
                                MuF*(Fox[pos-1]+Fox[gj*pos+1]+
                                Fox[getPos(offset-2,gj)]+lineDownF[gj]);
        }                
    }
    return(err);
}


int FillBorder(float  *Animal)
{
    int        err;
    int        i, j;
    err = 0;

    for( i=0; i <=NS_Size-2; i++) 
    {
        Animal[getPos(i,0)] = Animal[getPos(i,WE_Size-2)];
        Animal[getPos(i,WE_Size-1)] = Animal[getPos(i,1)];
    }

    for( j=0; j<=WE_Size-2; j++) 
    {
        Animal[getPos(0,j)] = Animal[getPos(NS_Size-2,j)];
        Animal[getPos(NS_Size-1,j)] = Animal[getPos(1,j)];
    }
    return(err);
}


int GetPopulation(float *Animal,float *tcount,int offset)
{
    int   err;
    int   i, j;
    float p;

    err = 0;
    p = 0.0;
    for( i=1; i < offset; i++)
      for( j=1;j<=WE_Size-2; j++)
            p = p + Animal[getPos(i,j)];

    *tcount = p;
    return(err);
}
