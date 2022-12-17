/// Nome dos integrantes do grupo
// Artur Auresco Damasio
// Isaac Santos Romão
// Rodrigo Augusto Alves Dias

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <pthread.h>

#define N 2048
#define REQ 0
#define RES 1

int processId;   /* rank do processo */
int noProcesses; /* Número de processos */
int **gridA, **gridB;

int GetOffset()
{
    int linesPerProccess = N / noProcesses;
    return processId * linesPerProccess;
}

void FillGlider(int **grid, int startLine, int endLine)
{
    // GLIDER
    int lin = 1, col = 1;

    if (lin >= startLine && lin <= endLine)
    {
        grid[lin - GetOffset()][col + 1] = 1;
    }

    if (lin + 1 >= startLine && lin + 1 <= endLine)
    {
        grid[lin + 1][col + 2] = 1;
    }

    if (lin + 2 >= startLine && lin + 2 <= endLine)
    {
        grid[lin + 2 - GetOffset()][col] = 1;
        grid[lin + 2 - GetOffset()][col + 1] = 1;
        grid[lin + 2 - GetOffset()][col + 2] = 1;
    }
}

void FillRPentonimo(int **grid, int startLine, int endLine)
{
    // R-pentomino
    int lin = 10;
    int col = 30;

    if (lin >= startLine && lin <= endLine)
    {
        grid[lin - GetOffset()][col + 1] = 1;
        grid[lin - GetOffset()][col + 2] = 1;
    }

    if (lin + 1 >= startLine && lin + 1 <= endLine)
    {
        grid[lin + 1 - GetOffset()][col] = 1;
        grid[lin + 1 - GetOffset()][col + 1] = 1;
    }

    if (lin + 2 >= startLine && lin + 2 <= endLine)
    {
        grid[lin + 2 - GetOffset()][col + 1] = 1;
    }
}

void PrintGrid(int **grid)
{
    int i, j;
    for (i = 0; i < 50; i++)
    {
        for (j = 0; j < 50; j++)
        {
            char c = grid[i][j] == 1 ? '*' : '.';
            printf("%c", c);
        }

        printf("\n");
    }
    printf("\n\n\n");
}

int GetAliveResponse(int gridIndex, int i, int j)
{
    int req[] = {gridIndex, i, j};
    int response;
    int targetProccess = i / (N / noProcesses);

    printf("Proccess %d requiring info from procces %d\n", processId, targetProccess);

    MPI_Status status;
    MPI_Send(&req, 3, MPI_INT, targetProccess, REQ, MPI_COMM_WORLD);
    MPI_Recv(&response, 1, MPI_INT, targetProccess, RES, MPI_COMM_WORLD, &status);

    return response;
}

int GetAliveNeighbors(int **grid, int line, int column, int startLine, int endLine, int gridIndex)
{
    int alive = 0;
    int j = 0;
    int current = 0;

    // obter índice de linha acima
    int above = line == 0 ? N - 1 : line - 1;

    // checar linha de cima
    for (j = column - 1; j <= column + 1; j++)
    {
        // checa borda infinita
        current = j;
        if (current < 0)
            current = N - 1;

        if (current > N - 1)
            current = 0;

        if (above >= startLine && above <= endLine)
        {
            if (grid[above][current] == 1)
                alive++;
        }
        else
        {
            alive += GetAliveResponse(gridIndex, above, current);
        }
    }

    // obter linha abaixo
    int below = (line + 1) % N;

    // checar linha de baixo
    for (j = column - 1; j <= column + 1; j++)
    {
        // checa borda infinita
        current = j;
        if (current < 0)
            current = N - 1;

        if (current > N - 1)
            current = 0;

        if (below >= startLine && below <= endLine)
        {
            if (grid[below][current] == 1)
                alive++;
        }
        else
        {
            alive += GetAliveResponse(gridIndex, below, current);
        }
    }

    // checar esquerda
    int left = column > 0 ? column - 1 : N - 1;
    if (grid[line][left] == 1)
        alive++;

    // checar direita
    int right = column < N - 1 ? column + 1 : 0;
    if (grid[line][right] == 1)
        alive++;

    return alive;
}

int GetNewState(int **grid, int line, int column, int startLine, int endLine, int gridIndex)
{
    int neighbors = GetAliveNeighbors(grid, line, column, startLine, endLine, gridIndex);

    // A célula está viva?
    if (grid[line][column] == 1)
    {
        // mantém-se viva
        if (neighbors == 2 || neighbors == 3)
            return 1;

        // deve morrer
        return 0;
    }

    // célula está morta
    // tem vizinhos suficientes para viver?
    if (neighbors == 3)
        return 1;

    return 0;
}

int GetSurvivors(int **grid, int startLine, int endLine, int gridIndex)
{
    int alive = 0;
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i >= startLine && i <= endLine)
            {
                alive += grid[i][j];
            }
            else
            {
                int i_j[] = {gridIndex, i, j};
                int response;
                int targetProccess = i / (N / noProcesses);
                MPI_Status status;

                MPI_Send(&i_j, 3, MPI_INT, targetProccess, REQ, MPI_COMM_WORLD);
                MPI_Recv(&response, 1, MPI_INT, targetProccess, RES, MPI_COMM_WORLD, &status);

                alive += response;
            }
        }
    }

    return alive;
}

int **GetCurrentGrid(int **gridA, int **gridB, int iteration)
{
    if (iteration % 2 == 0)
        return gridA;
    else
        return gridB;
}

int **GetNextGrid(int **gridA, int **gridB, int iteration)
{
    if (iteration % 2 == 0)
        return gridB;
    else
        return gridA;
}

void PlayGameOfLife(int **gridA, int **gridB, int iterations, int startLine, int endLine)
{
    int i, j, k;

    for (k = 0; k < iterations; k++)
    {
        int **nextGrid = GetNextGrid(gridA, gridB, k);
        int **currentGrid = GetCurrentGrid(gridA, gridB, k);

        for (i = 0; i < N; i++)
        {
            if (i >= startLine && i <= endLine)
            {
                for (j = 0; j < N; j++)
                {
                    nextGrid[i - GetOffset()][j] = GetNewState(currentGrid, i, j, startLine, endLine, k % 2);
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void *ReveiceMessage(void *reqProccess)
{
    int *pointer = (int *)reqProccess;
    int proccess = *pointer;
    free(pointer);
    while (1)
    {
        int response[] = {0, 0, 0};

        MPI_Status status;
        MPI_Recv(&response, 3, MPI_INT, proccess, REQ, MPI_COMM_WORLD, &status);

        printf("Proccess [%d] receive message from proccess [%d]\n", processId, proccess);

        int gridIndex = response[0];
        int i = response[1];
        int j = response[2];

        int aliveResponde = GetCurrentGrid(gridA, gridB, gridIndex)[i - GetOffset()][j];

        MPI_Send(&aliveResponde, 1, MPI_INT, proccess, RES, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[])
{
    int nameSize; /* Tamanho do nome */
    char computerName[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    MPI_Get_processor_name(computerName, &nameSize);

    int linesPerProccess = N / noProcesses;
    int startLine = processId * linesPerProccess;
    int endLine = startLine + linesPerProccess - 1;

    gridA = (int **)malloc(linesPerProccess * sizeof(int *));
    gridB = (int **)malloc(linesPerProccess * sizeof(int *));

    double start, end;
    pthread_t receive_thread;

    int i = 0, j = 0;
    for (i = 0; i < linesPerProccess; i++)
    {
        gridA[i] = (int *)malloc(N * sizeof(int));
        gridB[i] = (int *)malloc(N * sizeof(int));

        for (j = 0; j < N; j++)
        {
            gridA[i][j] = 0;
            gridB[i][j] = 0;
        }
    }

    FillGlider(gridA, startLine, endLine);
    FillRPentonimo(gridA, startLine, endLine);

    /*for (i = 0; i < noProcesses; i++)
    {
        if (i != processId)
        {
            int *id = (int *)malloc(sizeof(int));
            *id = i;
            pthread_create(&receive_thread, NULL, ReveiceMessage, (void *)id);
        }
    }*/

    printf("Proccess %d before barrier\n", processId);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Proccess %d after barrier\n", processId);

    if (processId == 0)
    {
        printf("*** Game of Life (OPEN MP)\n");
        printf("Condição inicial: %d\n", GetSurvivors(gridA, startLine, endLine, 0));
    }
    else
    {
        for (i = 0; i < noProcesses; i++)
        {
            if (i != processId)
            {
                int *id = (int *)malloc(sizeof(int));
                *id = i;
                pthread_create(&receive_thread, NULL, ReveiceMessage, (void *)id);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*start = omp_get_wtime();
    PlayGameOfLife(gridA, gridB, 2001, startLine, endLine);
    end = omp_get_wtime();*/

    if (processId == 0)
    {
        printf("Última geração (2000 iterações): %d\n", GetSurvivors(gridB, startLine, endLine, 1));
        printf("Tempo execução: %f\n", end - start);
    }
    else
    {
        pthread_exit(NULL);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}