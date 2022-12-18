/// Nome dos integrantes do grupo
// Artur Auresco Damasio
// Isaac Santos Romão
// Rodrigo Augusto Alves Dias

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define N 2048
#define REQL 0
#define SURV 1

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

int GetAliveNeighbors(int **grid, int line, int column, int startLine, int endLine, int *above, int *below)
{
    int alive = 0;
    int j = 0;
    int current = 0;

    // checar linha de cima
    for (j = column - 1; j <= column + 1; j++)
    {
        // checa borda infinita
        current = j;
        if (current < 0)
            current = N - 1;

        if (current > N - 1)
            current = 0;

        alive += above[current];
    }

    // checar linha de baixo
    for (j = column - 1; j <= column + 1; j++)
    {
        // checa borda infinita
        current = j;
        if (current < 0)
            current = N - 1;

        if (current > N - 1)
            current = 0;

        alive += below[current];
    }

    // checar esquerda
    int left = column > 0 ? column - 1 : N - 1;
    if (grid[line - GetOffset()][left] == 1)
        alive++;

    // checar direita
    int right = column < N - 1 ? column + 1 : 0;
    if (grid[line - GetOffset()][right] == 1)
        alive++;

    return alive;
}

int GetNewState(int **grid, int line, int column, int startLine, int endLine, int *above, int *below)
{
    int neighbors = GetAliveNeighbors(grid, line, column, startLine, endLine, above, below);

    // A célula está viva?
    if (grid[line - GetOffset()][column] == 1)
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

int GetSurvivors(int **grid, int startLine, int endLine)
{
    int alive = 0;
    int i, j;

    for (i = startLine; i <= endLine; i++)
    {
        for (j = 0; j < N; j++)
        {
            alive += grid[i - GetOffset()][j];
        }
    }

    for (int i = 1; i < noProcesses; i++)
    {
        int response;
        MPI_Status status;

        MPI_Recv(&response, 1, MPI_INT, i, SURV, MPI_COMM_WORLD, &status);

        alive += response;
    }

    return alive;
}

void SendSurvivors(int **grid, int startLine, int endLine)
{
    int alive = 0;
    int i, j;

    for (i = startLine; i <= endLine; i++)
    {
        for (j = 0; j < N; j++)
        {
            alive += grid[i - GetOffset()][j];
        }
    }

    MPI_Send(&alive, 1, MPI_INT, 0, SURV, MPI_COMM_WORLD);
}

void SendLine(int *line, int targetProccess)
{
    int lineToSend[N];
    int i;
    for (i = 0; i < N; i++)
    {
        lineToSend[i] = line[i];
    }

    MPI_Send(&lineToSend, N, MPI_INT, targetProccess, REQL, MPI_COMM_WORLD);
}

int *ReceiveLine(int proccess)
{
    int line[N];
    MPI_Status status;

    MPI_Recv(&line, N, MPI_INT, proccess, REQL, MPI_COMM_WORLD, &status);

    int *lineToReturn = (int *)malloc(sizeof(int) * N);
    int i;
    for (i = 0; i < N; i++)
    {
        lineToReturn[i] = line[i];
    }

    return lineToReturn;
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

        int *above, *below;

        if (processId == 0)
        {
            SendLine(currentGrid[endLine - GetOffset()], processId + 1);
            above = ReceiveLine(noProcesses - 1);

            SendLine(currentGrid[startLine - GetOffset()], noProcesses - 1);
            below = ReceiveLine(processId + 1);
        }
        else
        {
            above = ReceiveLine(processId - 1);
            SendLine(currentGrid[endLine - GetOffset()], (processId + 1) % noProcesses);

            below = ReceiveLine((processId + 1) % noProcesses);
            SendLine(currentGrid[startLine - GetOffset()], processId - 1);
        }

        for (i = startLine; i <= endLine; i++)
        {
            int *lineAbove = i == startLine ? above : currentGrid[i - GetOffset() - 1];
            int *lineBelow = i == endLine ? below : currentGrid[i - GetOffset() + 1];
            for (j = 0; j < N; j++)
            {
                nextGrid[i - GetOffset()][j] = GetNewState(currentGrid, i, j, startLine, endLine, lineAbove, lineBelow);
            }
        }

        free(above);
        free(below);
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

    if (processId == 0)
    {
        printf("*** Game of Life (OPEN MP)\n");
        printf("Condição inicial: %d\n", GetSurvivors(gridA, startLine, endLine));
    }
    else
    {
        SendSurvivors(gridA, startLine, endLine);
    }

    start = omp_get_wtime();
    PlayGameOfLife(gridA, gridB, 2001, startLine, endLine);
    end = omp_get_wtime();

    if (processId == 0)
    {
        printf("Última geração (2000 iterações): %d\n", GetSurvivors(gridB, startLine, endLine));
        printf("Tempo execução: %f\n", end - start);
    }
    else
    {
        SendSurvivors(gridB, startLine, endLine);
    }

    MPI_Finalize();

    return 0;
}