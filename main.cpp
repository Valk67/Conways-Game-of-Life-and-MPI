#include <mpi.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

using namespace std;

//ended up not using objects but still had my first grid an array of objects
struct cell {
  int corner = 0;
  int top = 0;
  int bottom = 0;
  int value;
};

// 1. Each cell with one or no neighbors dies, as if by solitude
// 2. Each cell with four or more neighbors dies, as if by overpopulation
// 3. Each cell with two or three neighbors survives.
// 4. Each cell with three neighbors becomes populated.
int evaluate(int CurCellValue, int TotalNeighbors) {
  if (TotalNeighbors < 2) 
    return 0;
  if (TotalNeighbors > 3) 
    return 0;
  if (TotalNeighbors == 2) 
    return CurCellValue;
  if (TotalNeighbors == 3) 
    return 1;

  return 99;
}

int main(int argc, char *argv[]) {
  if (argc != 2)
    cout << "usage: " << argv[0] << "<filename>\n";
  else {
    ifstream file(argv[1]);
    if (!file.is_open())
      cout << "Could not open file\n";
    else {
      clock_t start;
      double totalTime;
      start = clock();
      ofstream out("output.txt");

      int N = 0, G = 0, O = 0;
      int count = 0;
      string NGO;
      file >> N;
      file >> G;
      file >> O;

      cell **grid = new cell *[N];
      for (int i = 0; i < N; i++) {
        grid[i] = new cell[N];
      }

      // cout << endl;
      while (file >> NGO) {
        for (int i = 0; i < (int)NGO.size(); i++) {
          grid[count][i].value = (NGO[i] - '0');
        }
        count++;
      }
      // cout << "print of original grid" << endl;
      // for (int i = 0; i < N; i++) {
      //   for (int j = 0; j < N; j++) {
      //     cout << grid[i][j].value;
      //   }
      //   cout << '\n';
      // }
      // cout << '\n';


      // MPI initializers, etc from Dixons code.
      MPI_Init(NULL, NULL);
      // Get the number of processes
      int world_size;
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);




      /***** ERROR checking ******/

      if (N % world_size != 0) {
        cerr << "N is not divisible by num processors" << endl;
        exit(1);
      }





      // Get the rank of the process
      int world_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

      int numberSlices = N / world_size;
      // cout << "HEY ITS NUM SLICES = " << numberSlices << endl;
      int slices[numberSlices][N];


      // Declare the buffer and attach it

      MPI_Request request;


      if (world_rank == 0) {
        
        int boardArguements[4];
        boardArguements[0] = N;
        boardArguements[1] = G;
        boardArguements[2] = O;
        boardArguements[3] = numberSlices;

        for (int i = 0; i < world_size; i++) {
          MPI_Isend(&boardArguements, 4, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
          // MPI_Wait(&request, MPI_STATUS_IGNORE);

        }
        for (int i = 0; i < world_size; i++) {
          // MPI_Wait(&request, MPI_STATUS_IGNORE);

        }

        for (int i = 0; i < world_size; i++) {
          for (int j = 0; j < numberSlices; j++) {
            for (int k = 0; k < N; k++) {
              slices[j][k] = grid[j + (i * numberSlices)][k].value;	
            }
          }
          MPI_Isend(&slices, N*numberSlices, MPI_INT, i, 0, MPI_COMM_WORLD, &request);	
          // MPI_Wait(&request, MPI_STATUS_IGNORE);

	      }
          // for (int i = 0; i < world_size; i++) {

          //   MPI_Wait(&request, MPI_STATUS_IGNORE);
          // } 



      }  // end of world_rank == 0 for loop

      // Recieve Section
      int recivedInfo[4];
      MPI_Recv(&recivedInfo, 4, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
      int curSlice[recivedInfo[1]][recivedInfo[0]];
      MPI_Recv(&curSlice, recivedInfo[0] * recivedInfo[1], MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      N = recivedInfo[0];
      G = recivedInfo[1];
      O = recivedInfo[2];
      numberSlices = recivedInfo[3];

      int sendDown[N];
      int sendUp[N];
      int reciveDown[N];
      int reciveUp[N];

      //generation loop
      for (int generation = 1; generation < G; generation++) {  

        if (world_rank > 0) {
          for (int j = 0; j < N; j++) {
            sendUp[j] = curSlice[0][j];
          }
          MPI_Isend(&sendUp, N, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, &request);
          // MPI_Wait(&request, MPI_STATUS_IGNORE);
          MPI_Recv(&reciveUp, N, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        } else {
          for (int z = 0; z < N; z++) {
            reciveUp[z] = 0;
          }
        }
        if (world_rank != world_size - 1) {
          for (int j = 0; j < N; j++) {
            sendDown[j] = curSlice[numberSlices - 1][j];
          }
          MPI_Isend(&sendDown, N, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD, &request);
          // MPI_Wait(&request, MPI_STATUS_IGNORE);
          MPI_Recv(&reciveDown, N, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        } else {
          for (int k = 0; k < N; k++) {
            reciveDown[k] = 0;
          }
        }

        // cout << "\n";
        // cout << "printing recieveUP of " << world_rank << endl;
        // for (int i = 0; i < 8; i++) {
        //   cout << reciveUp[i] << " ";
        // }
        //   cout << "recieve down printed" << endl;
        //           cout << "\n";

        
        int neighbors = 0;  // sum of neighbours
        int tmpSlice[numberSlices][N];
        //colums are N rows are numberslices...
        // cout << "slices = " << numberSlices << " N = " << N << endl;
        for (int x = 0; x < numberSlices; x++) {
          // cout << "below it = " << reciveDown[x] << " " << reciveDown[x+1] << endl;
          for (int y = 0; y < N; y++)  {
            //CORNERS FIRST
            
            //top left
            if (x == 0 && y == 0) {
              neighbors = curSlice[x][y + 1] + reciveUp[y] + reciveUp[y + 1] + reciveDown[y] + reciveDown[y + 1];
              // cout << "sum worldrank "<< world_rank << "= " << sum << endl;
              // cout << "below it = " << reciveDown[0] << " " << reciveDown[1] << endl;
              // cout << "sum = " << sum << endl;

            }
            //top right of slice
            else if (x == 0 && y == N - 1) {  
              neighbors = curSlice[x][N - 2] + reciveDown[N - 1] + reciveDown[N - 2] + reciveUp[N - 1] + reciveUp[N - 2];
            }
            //bottom left of slice
            else if (x == numberSlices - 1 && y == 0) { 
              neighbors = curSlice[x][y + 1] + reciveUp[0] + reciveUp[1] + reciveDown[0] + reciveDown[1];
              // cout << "sum = " << sum << endl;
            }

            //bottom right of slice
            else if (x == numberSlices - 1 && y == N - 1) {
              neighbors = curSlice[x][N - 2] + reciveUp[N - 1] + reciveUp[N - 2] + reciveDown[N - 1] + reciveDown[N - 2];
            }
            else {

              //left of slice
              if (y == 0)  
                neighbors = reciveUp[0] + reciveUp[1] + reciveDown[0] + reciveDown[1] + curSlice[x][y + 1];

              //right of slice
              else if (y == N - 1) {
                neighbors = reciveUp[N - 1] + reciveUp[N - 2] + reciveDown[N - 1] + reciveDown[N - 2] + curSlice[x][N - 2];
              }

              //cells that are surrounded, or top and bottom
              else  {
                neighbors = reciveUp[y - 1] + reciveUp[y] + reciveUp[y + 1] +
                            reciveDown[y - 1] + reciveDown[y] + reciveDown[y + 1] +
                            curSlice[x][y - 1] + curSlice[x][y + 1];
              }
            }
            // cout << "sum = " << sum << endl;

            tmpSlice[x][y] = evaluate(curSlice[x][y], neighbors);
            neighbors = 0;

          }
        }
        
        for (int i = 0; i < numberSlices; i++) {
          for (int j = 0; j < N; j++) { 
            curSlice[i][j] = tmpSlice[i][j];
          }
        }        
        //print every O generation
        if (generation % O == 0) { 
          if (world_rank == 0) {
            
            int grid2[numberSlices][N];
            
            out << "\n";
            out << "Generation " << generation << endl;
            for (int i = 0; i < numberSlices; i++) {
              for (int j = 0; j < N; j++) {
                out << curSlice[i][j];
              }
             out << endl;
            }
            for (int i = 1; i < world_size; i++) {
              MPI_Recv(&grid2, N * numberSlices, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
              for (int j = 0; j < numberSlices; j++) {
                for (int k = 0; k < N; k++) {
                  out << grid2[j][k];
                }
                out << endl;
              }
            }

          } else {
            MPI_Isend(&curSlice, N * numberSlices, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            // MPI_Wait(&request, MPI_STATUS_IGNORE);
          }
        }

      }  //end of G loop

      
      //close mpi
      MPI_Finalize();
      totalTime = (clock() - start) / (double)CLOCKS_PER_SEC;
      cout << totalTime << endl;
    }
  }
  return 0;
}






//top 
              // else if (x == 0) 
              //   neighbors = reciveDown[y - 1] + reciveDown[y] + reciveDown[y + 1] + reciveUp[y - 1] + reciveUp[y] + reciveUp[y + 1] +
              //               curSlice[x][y - 1] + curSlice[x][y + 1];

              // //bottom
              // else if (x == numberSlices - 1)  
              //   neighbors = reciveUp[y - 1] + reciveUp[y] + reciveUp[y + 1] + reciveDown[y - 1] + reciveDown[y] + reciveDown[y + 1] +
              //               curSlice[x][y - 1] + curSlice[x][y + 1];

  