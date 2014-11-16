#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "param.h"
#include "macro.h"

float AlR = -0.2;
float BtR = 0.6;
float MuR = 0.01;
float BtF = 0.6;
float AlF = -1.8;
float MuF = 0.02;
char* s = "Passou";

char* filepathLRabbit = "/home/luis/git/Echo/fileRabbit.txt";
char* filepathLFox = "/home/luis/git/Echo/fileFox.txt";


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

void update(float* matrix1,float* matrix2,int line_offset)
{
    int i,j,pos;
    for (i = 0; i < line_offset; i++)
    {
        for (j = 0; j < WE_Size; j++)
        {
            pos=getPos(i,j);
            matrix1[pos] = matrix2[pos];
        }
    }
}

int SetLand ( float *Rabbit,float *Fox, int offset,int my_rank,int comm_sz);

int nonCriticalLines(float* Rabbit,float* Fox,float* TRabbit,float* TFox,int offset);

int criticalLines(float* lineUpR,float* lineDownR,float* lineUpF,float* lineDownF,
    float* Rabbit,float* Fox,float* TRabbit,float* TFox, int my_rank, int comm_sz,int offset);

int FillBorder(float  *Animal,int line_offset, float* lineUp, float* lineDown, int rank, int comm_size);

int GetPopulation(float *Animal,float *tcount,int offset,int rank, int comm_size);





int main(int argc, char *argv[]) 
{
    int rank,comm_size;
    MPI_Comm my_grid;
    int dim[2],period[2],reorder;
    int up,down,right,left,line_offset,err=0,k,j,pos;
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

    if(rank == comm_size -1)
    {
        line_offset = line_offset + (int)NS_Size % comm_size;
    }

    //Criação de dados derivados
    MPI_Datatype row;
    MPI_Type_contiguous((WE_Size), MPI_FLOAT, &row);
    MPI_Type_commit(&row);

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
    

    err = SetLand(lmRabbit,lmFox,line_offset,rank,comm_size);

    if(rank == comm_size - 5)
    {
        printMatrix(lmRabbit,line_offset,filepathLRabbit,rank);
        printMatrix(lmFox,line_offset,filepathLFox,rank);
    }

    for( k=1; k<=NITER; k++) 
    {

        nbrab=0;
        nbfox=0;
        totalrab=0;
        totalfox=0;

        if(rank == 0)
        {
            for(j=0; j <WE_Size;j++)
            {
                sendUp_rowR[j] = lmRabbit[getPos(1,j)];
            }

            MPI_Isend(sendUp_rowR, 1, row, comm_size - 1, 0, my_grid,reqR);

            MPI_Recv(recvDown_rowR,1, row, comm_size - 1, 0, my_grid, statusR);

            for(j=0; j <WE_Size;j++)
            {
                sendUp_rowF[j] = lmFox[getPos(1,j)];
            }

            MPI_Isend(sendUp_rowF, 1, row, comm_size - 1, 0, my_grid,reqF);

            MPI_Recv(recvDown_rowF,1, row, comm_size - 1, 0, my_grid, statusF);

            MPI_Wait(reqR,statusR);
            MPI_Wait(reqF,statusF);

        }

        if(rank == comm_size-1)
        {
            for(j=0; j <WE_Size;j++)
            {
                sendDown_rowR[j] = lmRabbit[getPos(line_offset-1,j)];
            }

            MPI_Isend(sendDown_rowR, 1, row, 0, 0, my_grid,reqR);

            MPI_Recv(recvUp_rowR,1, row, 0, 0, my_grid, statusR);

            for(j=0; j <WE_Size;j++)
            {
                sendDown_rowF[j] = lmFox[getPos(line_offset-1,j)];
            }

            MPI_Isend(sendDown_rowF, 1, row, 0, 0, my_grid,reqF);

            MPI_Recv(recvUp_rowF,1, row, 0, 0, my_grid, statusF);

            MPI_Wait(reqR,statusR);
            MPI_Wait(reqF,statusF);                        
        }


        err = FillBorder(lmRabbit,line_offset,recvUp_rowR,recvDown_rowR,rank,comm_size); 
        err = FillBorder(lmFox,line_offset,recvUp_rowF,recvDown_rowF,rank,comm_size);


        MPI_Cart_shift(my_grid,1,1,&up,&down);

        if (up == -1)
        {
            down = MPI_PROC_NULL;
        }
        if (down == comm_size)
        {
            down = MPI_PROC_NULL;
        }

        //////UP Lines

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(0,j);
            sendUp_rowR[j] = lmRabbit[pos];
        }

        MPI_Isend(sendUp_rowR, 1, row, up, 0, my_grid,reqR);


        MPI_Recv(recvDown_rowR,1, row, down, 0, my_grid, statusR);

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(0,j);
            sendUp_rowF[j] = lmFox[pos];
        }
        MPI_Isend(sendUp_rowF, 1, row, up, 0, my_grid,reqF);

        MPI_Recv(recvDown_rowF,1, row, down, 0, my_grid, statusF);

        //////DOWN Lines

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(line_offset-1,j);
            sendDown_rowR[j] = lmRabbit[pos];
        }
        MPI_Isend(sendDown_rowR, 1, row, down, 0, my_grid,reqR);

        MPI_Recv(recvUp_rowR,1, row, up, 0, my_grid, statusR);

        for(j = 0; j < (WE_Size); j++)
        {
            pos=getPos(line_offset-1,j);
            sendDown_rowF[j] = lmRabbit[pos];
        }
        MPI_Isend(sendDown_rowF, 1, row, down, 0, my_grid,reqF);

        MPI_Recv(recvUp_rowF,1, row, up, 0, my_grid, statusF);

        err = nonCriticalLines(lmRabbit,lmFox,lmTRabbit,lmTFox,line_offset);

        MPI_Wait(reqR,statusR);
        MPI_Wait(reqF,statusF);


        err = criticalLines(recvUp_rowR,recvDown_rowR,recvUp_rowF,recvDown_rowF,lmRabbit,lmFox,lmTRabbit,lmTFox,rank,comm_size,line_offset);           

        if( (k % PERIOD) == 1 )
        {
            err = GetPopulation(lmTRabbit,&nbrab,line_offset,rank,comm_size); 
            err = GetPopulation(lmTFox,&nbfox,line_offset,rank,comm_size);
            MPI_Reduce(&nbrab,&totalrab,1,MPI_FLOAT,MPI_SUM,0,my_grid);
            MPI_Reduce(&nbfox,&totalfox,1,MPI_FLOAT,MPI_SUM,0,my_grid);
            if(rank == 0)
                printf("In the year %d, the number of rabbits is %d and the number of foxes is %d. \n", k,(int)totalrab,(int)totalfox);
        }
        update(lmRabbit,lmTRabbit,line_offset);
        update(lmFox,lmTFox,line_offset);

    }
    free(sendUp_rowR);
    free(sendDown_rowR);
    free(sendUp_rowF);
    free(sendDown_rowF);
    free(recvUp_rowR);
    free(recvDown_rowR);
    free(recvUp_rowF);
    free(recvDown_rowF);
    free(lmRabbit);
    free(lmFox);
    free(lmTRabbit);
    free(lmTFox);
        
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
                            MuR*(Rabbit[pos-1]+Rabbit[pos+1]+Rabbit[getPos(gi-1,gj)]+Rabbit[getPos(gi+1,gj)]);

            TFox[pos] = AlF*Rabbit[pos] +
                         (1.0+BtF-4.0*MuF)*Fox[pos] +
                         MuF*(Fox[pos-1]+Fox[pos+1]+Fox[getPos(gi-1,gj)]+Fox[getPos(gi+1,gj)]);
        }
    }

    for( gi=1; gi < offset-1; gi++)
    {
        for(gj=1; gj<=WE_Size-3; gj++ )
        {
            pos=getPos(gi,gj);
            TRabbit[pos] = MAX( 0.0, TRabbit[pos]);
            TFox[pos] = MAX( 0.0, TFox[pos]);
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
            TRabbit[pos] = MAX( 0.0, TRabbit[pos]);
            TFox[pos] = MAX( 0.0, TFox[pos]);

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
            TRabbit[pos] = MAX( 0.0, TRabbit[pos]);
            TFox[pos] = MAX( 0.0, TFox[pos]);

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

            TRabbit[pos] = MAX( 0.0, TRabbit[pos]);
            TFox[pos] = MAX( 0.0, TFox[pos]);
        }
        for( gj=1; gj<=WE_Size-2; gj++)
        {
            pos=getPos(offset-1,gj);
            TRabbit[pos] = (1.0+AlR-4.0*MuR)*Rabbit[pos] + BtR*Fox[pos] +
                                MuR*(Rabbit[pos-1] + Rabbit[pos+1]
                                + Rabbit[getPos(offset-2,gj)]+lineDownR[gj]);

            TFox[pos] = AlF*Rabbit[pos] +
                                (1.0+BtF-4.0*MuF)*Fox[pos] +
                                MuF*(Fox[pos-1]+Fox[pos+1]+
                                Fox[getPos(offset-2,gj)]+lineDownF[gj]);
            TRabbit[pos] = MAX( 0.0, TRabbit[pos]);
            TFox[pos] = MAX( 0.0, TFox[pos]);                    
        }               
    }
    return(err);
}


int FillBorder(float  *Animal,int line_offset, float* lineUp, float* lineDown, int rank, int comm_size)
{
    int        err;
    int        i, j;
    err = 0;

    if(rank == 0)
    {

        for( j=0; j<=WE_Size-2; j++) 
        {
            Animal[getPos(0,j)] = lineUp[j];
        }
    }
    if (rank +1 == comm_size)
    {
        for( j=0; j<=WE_Size-2; j++) 
        {
            Animal[getPos(line_offset-1,j)] = lineDown[j];
        }
    }    
    for( i=0; i < line_offset; i++) 
    {
        Animal[getPos(i,0)] = Animal[getPos(i,WE_Size-2)];
        Animal[getPos(i,WE_Size-1)] = Animal[getPos(i,1)];
    }
    return(err);
}


int GetPopulation(float *Animal,float *tcount,int offset,int rank,int comm_size)
{
    int   err;
    int   i, j;
    float p;

    err = 0;
    p = 0.0;
    if(rank == 0)
    {
        for( i=1; i < offset; i++)
            for( j=1;j<=WE_Size-2; j++)
                p = p + Animal[getPos(i,j)];
    }
    if(rank == comm_size-1)
    {
        for( i=0; i < offset -1; i++)
            for( j=1;j<=WE_Size-2; j++)
                p = p + Animal[getPos(i,j)];
    }
    else
    {
        for( i=0; i < offset; i++)
            for( j=1;j<=WE_Size-2; j++)
                p = p + Animal[getPos(i,j)];
    }
    *tcount = p;
    return(err);
}
