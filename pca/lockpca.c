/* this version of the cellular automata currently stores ints,
* just a test version to get it up and running, will convert to
* chars later, but most importantly also runs in parallel
*
* known issues: randomness is not properly initialized
* when finding cells of same distance, chooses in an order
* 
*/
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <omp.h>
#include <string.h>

// will not be needed later
#define GRIDWIDTH 20
#define GRIDHEIGHT 20
#define GRIDBUFFER 1
#define TOTALWIDTH (GRIDWIDTH + 2 * GRIDBUFFER)
#define TOTALHEIGHT (GRIDHEIGHT + 2 * GRIDBUFFER)
#define MAXGRIDVAL 100

typedef struct{
  int val;
  int x;
  int y;
  char flag;
} cell_t;


int get_valid_random_number(short unsigned int* thread_data, int low, int high){
// issue with assigning val and some sort of overwriting happening
//  double val;
//  #pragma omp critical
  double val = erand48(thread_data);
  val = val * (high - low) + low;
  return (int) val;
  //return (int) arc4random_uniform(high - low) + low;
}

cell_t* init_grid(int width, int height, short unsigned int** data){
  cell_t* grid_start = malloc(width * height * sizeof(cell_t));
  #pragma omp parallel for
  for (int i = 0; i < width * height; i++){
    cell_t cell;
    if ((i % width) == 0 || (i / width) == 0 || (i % width) == (width - 1) || (i / width) >= (width - 1))
      cell.val = 0;
    else {
      cell.val = get_valid_random_number(*(data + omp_get_thread_num()), 1, MAXGRIDVAL);
    }
    cell.flag = 'f';
    cell.x = i / width;
    cell.y = i % width;
    *(grid_start + i) = cell;
  }
  return grid_start;
}



// the following methods will never be called on an edge cell
// so don't need to perform any edge checking

// returns the cell above the current cell in the grid
cell_t* get_above_cell(cell_t* cell, int width){
  return (cell - width); 
}

// returns the cell below the current cell in the grid
cell_t* get_below_cell(cell_t* cell, int width){
  return (cell + width); 
}

// returns the cell to the left of the current cell in the grid
cell_t* get_left_cell(cell_t* cell){
  return (cell - 1); 
}

// returns the cell to the right of the current cell in the grid
cell_t* get_right_cell(cell_t* cell){
  return (cell + 1); 
}

// selects a random cell inside of the grid using a random number generator
cell_t* select_random_cell(cell_t* grid, int width, int height, short unsigned int* thread_data){
  cell_t* cell = (grid + get_valid_random_number(thread_data, 1, width * height));
  while (cell->val == 0){
    cell = (grid + get_valid_random_number(thread_data, 1, width * height));
  }
  return cell;
}


// returns the difference between the value of the cell passed in
// and the cell above it in the grid
int above_cell_distance(cell_t* cell, int width){
  cell_t* above_cell = get_above_cell(cell, width);
  if (above_cell->val == 0) return INT_MAX;
  return abs(cell->val - above_cell->val); 
}

// returns the difference between the value of the cell passed in
// and the cell below it in the grid
int below_cell_distance(cell_t* cell, int width){
  cell_t* below_cell = get_below_cell(cell, width);
  if (below_cell->val == 0) return INT_MAX;
  return abs(cell->val - below_cell->val);
}

// returns the difference between the value of the cell passed in
// and the cell to the left of it in the grid
int left_cell_distance(cell_t* cell){
  cell_t* left_cell = get_left_cell(cell);
  if (left_cell->val == 0) return INT_MAX;
  return abs(cell->val - left_cell->val);
}

// returns the difference between the value of the cell passed in
// and the cell to the right of it in the grid
int right_cell_distance(cell_t* cell){
  cell_t *right_cell = get_right_cell(cell);
  if (right_cell->val == 0) return INT_MAX;
  return abs(cell->val - right_cell->val);
}

cell_t* get_random_neighbor(cell_t* curr_cell, int width){
  int rand = arc4random_uniform(4);
  if (rand == 0) return get_above_cell(curr_cell, width);
  if (rand == 1) return get_below_cell(curr_cell, width);
  if (rand == 2) return get_left_cell(curr_cell);
  if (rand == 3) return get_right_cell(curr_cell);
}

// can parallelize using sections ???
// returns the neighbor with the smallest difference in value with the given cell
// if two neighbors are the same, it doesn't effectively choose one yet
cell_t* get_best_transition_neighbor(cell_t* curr_cell, int width){
  int min = MAXGRIDVAL;
  cell_t* neighbor;
  int same_dist = 0;

  if (above_cell_distance(curr_cell, width) < min && get_above_cell(curr_cell, width)->val != 0){
    min = above_cell_distance(curr_cell, width);
    neighbor = get_above_cell(curr_cell, width);
  }
  if (below_cell_distance(curr_cell, width) < min && get_below_cell(curr_cell, width)->val != 0){
    min = below_cell_distance(curr_cell, width);
    neighbor = get_below_cell(curr_cell, width);
  } else if (below_cell_distance(curr_cell, width) == min && get_below_cell(curr_cell, width)->val != 0){
    same_dist = 1;
  }
  if (left_cell_distance(curr_cell) < min && get_left_cell(curr_cell)->val != 0){
    min = left_cell_distance(curr_cell);
    neighbor = get_left_cell(curr_cell);
    same_dist = 0;
  } else if (left_cell_distance(curr_cell) == min && get_left_cell(curr_cell)->val != 0){
    same_dist = 1;
  }
  if (right_cell_distance(curr_cell) < min && get_right_cell(curr_cell)->val != 0){
    min = right_cell_distance(curr_cell);
    neighbor = get_right_cell(curr_cell);
    same_dist = 0;
  }  else if (right_cell_distance(curr_cell) == min && get_right_cell(curr_cell)->val != 0){
    same_dist = 1;
  }
  if (same_dist) {
    neighbor = get_random_neighbor(curr_cell, width);
    while (neighbor->val == 0)  neighbor = get_random_neighbor(curr_cell, width);
  }
  return neighbor;
}

// swaps the value of the current cell and the neighbor it should transition with
void swap_cells_with_display(cell_t* curr_cell, int width , cell_t swaps[TOTALWIDTH][TOTALHEIGHT]){
  cell_t* neighbor = get_best_transition_neighbor(curr_cell, width);
  cell_t* temp = malloc(sizeof(cell_t));

  while (1){
    char c, n;
    #pragma omp atomic read
      c = curr_cell->flag;
    #pragma omp atomic read
      n = neighbor->flag;
    if (c == 'f' && n == 'f') break;
    // flush?
  }
  #pragma atomic write
  {
    curr_cell->flag = 'u';
    neighbor->flag = 'u';
  }
  memmove(temp, curr_cell, sizeof(cell_t));
  memmove(curr_cell, neighbor, sizeof(cell_t));
  memmove(neighbor, temp, sizeof(cell_t));
  #pragma atomic write
  {
    curr_cell->flag = 'f';
    neighbor->flag = 'f';
  }
  free(temp);
  swaps[curr_cell->x][curr_cell->y].val++;
  swaps[neighbor->x][neighbor->y].val++;
}

// swaps the value of the current cell and the neighbor it should transition with
void swap_cells(cell_t* curr_cell, int width){
  cell_t* neighbor = get_best_transition_neighbor(curr_cell, width);
  cell_t* temp = malloc(sizeof(cell_t));

  while (1){
    char c, n;
    #pragma omp atomic read
      c = curr_cell->flag;
    #pragma omp atomic read
      n = neighbor->flag;
    if (c == 'f' && n == 'f') break;
    // flush?
  }
  #pragma atomic write
  {
    curr_cell->flag = 'u';
    neighbor->flag = 'u';
  }
//  #pragma omp flush((curr_cell->flag), (neighbor->flag))
  memmove(temp, curr_cell, sizeof(cell_t));
  memmove(curr_cell, neighbor, sizeof(cell_t));
  memmove(neighbor, temp, sizeof(cell_t));
  #pragma atomic write
  {
    curr_cell->flag = 'f';
    neighbor->flag = 'f';
  }
//  #pragma omp flush(curr_cell->flag)
//  #pragma omp flush(neighbor->flag)
  free(temp);
}


void display_swaps(cell_t grid[TOTALWIDTH][TOTALHEIGHT]){
  int row_sums = 0;
  for (int i = 0; i < TOTALWIDTH; i++){
    int sum = 0;
    for (int j = 0; j < TOTALHEIGHT; j++){
      printf(" %d ", grid[i][j].val);
      sum += grid[i][j].val;
    }
    row_sums += sum;
    printf("     = %d \n", sum);
  }

  for (int i = 0; i < TOTALHEIGHT; i++){
    int sum = 0;
    for (int j = 0; j < TOTALWIDTH; j++){
      sum += grid[j][i].val;
    }
    printf(" %d ", sum);
  }
  printf("\n");
  printf("row sums: %d\n", row_sums);
}


void display_grid(cell_t* start, int width, int height){
  for(int i = 0; i < (width * height); i++){
    if (i % width == 0) printf("\n");
    cell_t cell = (*(start + i));
    if (cell.val != 0) printf(" %d  ", cell.val);
  }
  
}



// simple cellular automata that is executed in parallel with 8 (default) threads
// main method initializes the array (c doesn't like passing arrays around so there
// will be another file that uses pointers instead of a 2d array later) 
// and also executes a given number of iterations for the cellular automata to go through
// the serial version only does one swap at a time (like rescal), so the execution is split
// up into the threads based on the number of iterations it needs to go through
// the next thing that needs to be implemented is how to better use the threads 
// i.e. split up based on region of grid for swapping
int main(int argc, const char* argv[]){
  double start_time = omp_get_wtime();
  int seed = (int) arc4random();
  int width, height, num_iter;
  
  if (argc < 4) {
    width = 20;
    height = 20;
    num_iter = 1000;
  } else {
    width = atoi(argv[1]);
    height = atoi(argv[2]);
    num_iter = atoi(argv[3]);
  }


  short unsigned int** erand_thread_data;
  #pragma omp parallel
  {
    #pragma omp single
    {
      erand_thread_data = malloc(omp_get_num_threads() * sizeof(short unsigned int *));
    }
    
    //double* result;
    short unsigned int* data = calloc(3, sizeof(short unsigned int));
    data[0] = arc4random(); 
    data[1] = arc4random(); 
    data[2] = arc4random(); 
    *(erand_thread_data + omp_get_thread_num()) = data;
  }
  cell_t* grid_start = init_grid(width, height, erand_thread_data);
  
  //display_grid(grid_start, width, height);


  cell_t swaps[TOTALWIDTH][TOTALHEIGHT];
  #pragma omp parallel for
  for (int i = 0; i < TOTALWIDTH; i++){
    for (int j = 0; j < TOTALHEIGHT; j++){
      cell_t newcell;
      newcell.x = i;
      newcell.y = j;
      newcell.val = 0;
      swaps[i][j] = newcell;
    }
  }
  
  
   
// BUGFIX - threads stepping on each other when swapping cells and/or selecting random #'s
  #pragma omp parallel for
  for (int i = 0; i < num_iter; i++){
    cell_t* cell = select_random_cell(grid_start, width, height, 
          *(erand_thread_data + omp_get_thread_num())); 
    swap_cells_with_display(cell, width , swaps);
//    swap_cells(cell, width);
  }
//  display_grid(grid_start, width, height);
  display_swaps(swaps);

/*  
  int* count = calloc(MAXGRIDVAL, sizeof(int)); 
  for (int i = 0; i < width * height; i++){
    (*(count + ((grid_start + i)->val)))++;
  }
  for (int i = 0; i < MAXGRIDVAL; i++){
    printf("%d: %d \n", i, *(count + i) - 4);
  }
*/

  printf("Parallel: Time to completion for %d iterations of a size %d x %d grid: %f \n", 
      num_iter, width, height, (omp_get_wtime() - start_time));
  free(grid_start);
  
  return 0;
}
