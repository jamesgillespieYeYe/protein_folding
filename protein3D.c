#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>

//Options
//#define C1
//#define C2
//#define WATER_PENALTY
#define DIM 3

//Bound defines the range of starting positions
#if (DIM / 2) < (DIM - DIM / 2)
#   define BOUND DIM - DIM / 2
#else
#   define BOUND DIM / 2
#endif

//WP is the Water Penalty Amount
#define WP 2.5

//Values for contacts
// #define Hhh -2
// #define Hhp -1
// #define Hpp -0.1

#define Hhh -2.3
#define Hhp -1
#define Hpp 0



//#define NUM_ACIDS 20
#define NUM_ACIDS 20
#define NUM_HYDROPHILIC 7

struct acid 
{
    char name[4];
    bool hydrophilic;
    struct acid * prev;
    struct acid * next;
};
typedef struct acid acid;

struct coordinate
{
    int x;
    int y;
    int z;
    double diagonal;
    int valid;
};
typedef struct coordinate coordinate;

struct thread_arg
{
    int x;
    int y;
    int z;
    acid * pgrid[DIM][DIM][DIM];
    acid * acids_list[NUM_ACIDS];
    long num_valid_configurations;
    double min_energy;
    acid * optimal_configuration[DIM][DIM][DIM];
    long num_recursive_calls;
    int num_new_mins;
};
typedef struct thread_arg thread_arg;





coordinate * get_positions(int x, int y, int z);
coordinate * Spatial_Database[DIM][DIM][DIM];
void print_grid(acid * pgrid[DIM][DIM][DIM]);


#ifdef WATER_PENALTY
struct edge
{
    //0: interior 
    //1: exposed on one side 
    //2: exposed on two sides 
    //3: exposed on 3 sides
    int edge_type;
};
typedef struct edge edge;

//Give External WP based on edge type
double EXT_WP_Lookup_Table[4] = 
{
    0.0, //Type 0 - all adjacencies are within the grid
    WP + 4*WP/2 + 4*WP/3, //Type 1 : 1 direct adjacency, 4 C1 adjaencies, 4 C2 adjacencies
    2*WP + 7*WP/2 + 6*WP/3,                     //Type 2:  2 direct, 9 C1 diagonals, 4 C2 diagonals
    3*WP + 9*WP/2 + 7*WP/3  //3, 10, 6
};

edge * Edge_Database[DIM][DIM][DIM];

edge * construct_edge(int x, int y, int z)
{
    edge * ret = calloc(sizeof(edge), 1);

    int coors[3] = {x, y, z};
    for (int i = 0; i < 3; i++)
    {
        int curr_coor = coors[i];
        if (curr_coor == 0 || curr_coor == DIM - 1)
        {
            ret->edge_type++;
        }
    }

    return ret;
}

double Water_Penalty(acid * pgrid[DIM][DIM][DIM], int x, int y, int z)
{
    //Get external penalty
    edge * currEdge = Edge_Database[x][y][z];
    //printf("Edge type: %d\n", currEdge->edge_type);
    double external_penalty = EXT_WP_Lookup_Table[currEdge->edge_type];
    //printf("External penalty: %f\n", external_penalty);
    //Get interal penalty 
    coordinate * positions = Spatial_Database[x][y][z];
    coordinate * tmp = positions;
    double internal_penalty = 0;
    while (tmp->valid == 1 && tmp->diagonal == 0)
    {
        //printf("(%d, %d, %d)\n", tmp->x, tmp->y, tmp->z);
        if (pgrid[tmp->x][tmp->y][tmp->z] == NULL)
        {
            internal_penalty += WP;
        }
        tmp++;
    }
    while (tmp->valid == 1 && tmp->diagonal == 1)
    {
        //printf("(%d, %d, %d)\n", tmp->x, tmp->y, tmp->z);
        if (pgrid[tmp->x][tmp->y][tmp->z] == NULL)
        {
            internal_penalty += WP / 2;
        }
        tmp++;
    }
    while (tmp->valid == 1 && tmp->diagonal == 2)
    {
        //printf("(%d, %d, %d)\n", tmp->x, tmp->y, tmp->z);
        if (pgrid[tmp->x][tmp->y][tmp->z] == NULL)
        {
            internal_penalty += WP / 3;
        }
        tmp++;
    }



    

    double total_penalty = internal_penalty + external_penalty;
    

    return total_penalty;



}

double Total_Water_Penalty(acid * pgrid[DIM][DIM][DIM])
{
    double total;;
    for (int x = 0; x < DIM; x++)
    {
        for (int y = 0; y < DIM; y++)
        {
            for (int z = 0; z < DIM; z++)
            {
                acid * currAcid = pgrid[x][y][z];
                if (currAcid != NULL && currAcid->hydrophilic == false)
                {
                    total += Water_Penalty(pgrid, x, y, z);
                }
            }
        }
    }

    return total;
}
#endif


/**
 * Give the energy between the acid pair (a, b)
*/
double Energy(acid * a, acid * b)
{
    if (a->hydrophilic == false && b->hydrophilic == false)
    {
        return Hhh;
    }
    else if (a->hydrophilic == true && b->hydrophilic == false)
    {
        return Hhp;
    }
    else if (a->hydrophilic == false && b->hydrophilic == true)
    {
        return Hhp;
    }
    else 
    {
        return Hpp;
    }
}

/**
 * Compute the total energy of a configuration by examining energies of 
 * all adjacencies
*/
double Total_Energy(acid * pgrid[DIM][DIM][DIM])
{
    double ret = 0;
    acid * grid[DIM][DIM][DIM];
    memcpy(grid, pgrid, sizeof(acid*)*DIM*DIM*DIM);
    for (int x = 0; x < DIM; x++)
    {
        for (int y = 0; y < DIM; y++)
        {
            for (int z = 0; z < DIM; z++)
            {
                acid * curr = grid[x][y][z];
                if (curr != NULL)
                {
                    
                    coordinate * positions = Spatial_Database[x][y][z];
                    coordinate * tmp = positions;
                    while (tmp->valid == 1 && tmp->diagonal == 0)
                    {
                        acid * adjacent = grid[tmp->x][tmp->y][tmp->z];
                        if (adjacent != NULL)
                        {
                            ret += Energy(curr, adjacent);
                        }
                        tmp++;
                    }
#                   ifdef C1
                    while (tmp->valid == 1 && tmp->diagonal == 1)
                    {
                        acid * adjacent = grid[tmp->x][tmp->y][tmp->z];
                        if (adjacent != NULL)
                        {
                            ret += Energy(curr, adjacent) / 2;
                        }
                        tmp++;
                    }
#                   endif
#                   ifdef C2
                    while (tmp->valid == 1 && tmp->diagonal == 2)
                    {
                        acid * adjacent = grid[tmp->x][tmp->y][tmp->z];
                        if (adjacent != NULL)
                        {
                            ret += Energy(curr, adjacent) / 3;
                        }
                        tmp++;
                    }
#                   endif
                    grid[x][y][z] = NULL;

                }
            }
        }
    }
#   ifdef WATER_PENALTY
    ret += Total_Water_Penalty(pgrid);
#   endif
    return ret;
}



char * Sequence[NUM_ACIDS] = {"ASN", "LEU", "TYR", "ILE", "GLN", "TRP", "LEU", "LYS", "ASP", "GLY", "GLY", "PRO", "SER", "SER", "GLY", "ARG", "PRO", "PRO", "PRO", "SER"};
//char * Sequence[NUM_ACIDS] = {"ASN", "LEU", "TYR"};
char * Hydrophilic[NUM_HYDROPHILIC] = {"ASP", "ARG", "LYS", "GLY", "ASN", "SER", "GLN"};



#define PDB_FILENAME "out.pdb"
void write_PDB(acid * pgrid[DIM][DIM][DIM])
{
    FILE * out = fopen(PDB_FILENAME, "w");
    int label = 1;
    for (int x = 0; x < DIM; x++)
    {
        for (int y = 0; y < DIM; y++)
        {
            for (int z = 0; z < DIM; z++)
            {
                acid * curr = pgrid[x][y][z];
                if (curr != NULL)
                {
                    fprintf(out, "ATOM     %d CA   %s     %d      %d.00     %d.00     %d.00  0.00  0.00           C\n", label, "CA", label, 5*x, 5*y, 5*z);
                    label++;
                }
            }
        }
    }
    fclose(out);
}


void write_PDB_TWO(acid * pgrid[DIM][DIM][DIM], acid * acids[NUM_ACIDS])
{
    FILE * out = fopen(PDB_FILENAME, "w");
    for (int index = 0; index < NUM_ACIDS; index++)
    {
        acid * curr = acids[index];
        for (int x = 0; x < DIM; x++)
        {
            for (int y = 0; y < DIM; y++)
            {
                for (int z = 0; z < DIM; z++)
                {
                    if (curr == pgrid[x][y][z])
                    {
                        int label = index + 1;
                        fprintf(out, "ATOM     %d CA   %s     %d      %d.00     %d.00     %d.00  0.00  0.00           C\n", label, curr->name, label, 1*x, 1*y, 1*z);
                    }
                }
            }
        }
    }
    // for (int num = 1; num < NUM_ACIDS; num++)
    // {
    //     fprintf(out, "CONNECT %d %d\n", num, num + 1);
    // }
    fclose(out);
}


void print_grid(acid * pgrid[DIM][DIM][DIM])
{
    for (int z = 0; z < DIM; z++)
    {
        printf("Slice: z = %d\n", z);
        for (int x = 0; x < DIM; x++)
        {
            for (int y = 0; y < DIM; y++)
            {
                acid * curr = pgrid[x][y][z];
                if (curr == NULL)
                {
                    printf("NULL ");
                }
                else
                {
                    printf("%s  ", curr->name);
                }
            }
            printf("\n");
        }
    }
}






/**
 * Recursive callback
 * Look at all positions adjacent to (lastX, lastY, lastZ), and try placing next acid there if possible
*/
void insert(acid * pgrid[DIM][DIM][DIM], acid * acids_list[NUM_ACIDS], int index, int lastX, int lastY, int lastZ, thread_arg * args)
{
    args->num_recursive_calls++;
    if (index == NUM_ACIDS)
    {
        args->num_valid_configurations++;
        if (args->num_valid_configurations % 10000000 == 0)
        {
            printf("(Thread %ld) Base case number: %ld\n", pthread_self(), args->num_valid_configurations);
        }
        double energy = Total_Energy(pgrid);
        if (energy < args->min_energy)
        {
            args->min_energy = energy;
            memcpy(args->optimal_configuration, pgrid, sizeof(acid*)*DIM*DIM*DIM);
            args->num_new_mins++;
        }
        return;
    }

    acid * next_acid = acids_list[index];
    coordinate * moves = Spatial_Database[lastX][lastY][lastZ];
    coordinate * pmove = moves;

    while (pmove->valid == 1 && pmove->diagonal == 0)
    {   
        
        if (pgrid[pmove->x][pmove->y][pmove->z] == NULL)
        {
            acid * next_grid[DIM][DIM][DIM];
            memcpy(next_grid, pgrid, sizeof(acid*)*DIM*DIM*DIM);
            next_grid[pmove->x][pmove->y][pmove->z] = acids_list[index];

            insert(next_grid, acids_list, index + 1, pmove->x, pmove->y, pmove->z, args);
        }
        pmove++;
    }
}

#define GRID_FILENAME "grid.out.txt"
void save_grid(acid * pgrid[DIM][DIM][DIM], acid * acids[NUM_ACIDS])
{
    FILE * f = fopen(GRID_FILENAME, "w");
    fprintf(f, "acid,x,y,z,hydrophilic\n");
    for (int index = 0; index < NUM_ACIDS; index++)
    {
        acid * curr_acid = acids[index];
        for (int x = 0; x < DIM; x++)
        {
            for (int y = 0; y < DIM; y++)
            {
                for (int z = 0; z < DIM; z++)
                {
                    if (curr_acid == pgrid[x][y][z])
                    {
                        if (curr_acid->hydrophilic == true)
                        {
                            fprintf(f, "%s,%d,%d,%d,True\n", curr_acid->name, x, y, z);
                        }
                        else
                        {
                            fprintf(f, "%s,%d,%d,%d,False\n", curr_acid->name, x, y, z);
                        }
                        
                    } 
                }
            }
        }
    }
    fclose(f);
}




void* thread_func(void* varg)
{
    thread_arg * arg = (thread_arg*) varg;
    printf("New thread: %ld: %d, %d, %d\n", pthread_self(), arg->x, arg->y, arg->z);
    arg->pgrid[arg->x][arg->y][arg->z] = arg->acids_list[0];
    insert(arg->pgrid, arg->acids_list, 1, arg->x, arg->y, arg->z, arg);
    return NULL;
}



int main(int argc, char** argv)
{
#   ifdef WATER_PENALTY
    for (int x = 0; x < DIM; x++)
    {
        for (int y = 0; y < DIM; y++)
        {
            for (int z = 0; z < DIM; z++)
            {
                Edge_Database[x][y][z] = construct_edge(x, y, z);
            }
        }
    }
#   endif
    
    acid * grid[DIM][DIM][DIM];
    memset((void*) grid, 0, sizeof(acid*)*DIM*DIM*DIM);
    print_grid(grid);

    //Assemble Spatial Database
    for (int x = 0; x < DIM; x++)
    {
        for (int y = 0; y < DIM; y++)
        {
            for (int z = 0; z < DIM; z++)
            {
                Spatial_Database[x][y][z] = get_positions(x, y, z);
            }
        }
    }

    //Assemble acids
    acid * acids[NUM_ACIDS];

    acid * prev = NULL;

    for (int i = 0; i < NUM_ACIDS; i++)
    {
        acids[i] = calloc(sizeof(acid), 1);
        strcpy(acids[i]->name, Sequence[i]);
        acids[i]->hydrophilic = false;
        if (prev != NULL)
        {
            acids[i]->prev = prev;
            prev->next = acids[i];
        }
        for (int j = 0; j < NUM_HYDROPHILIC; j++)
        {
            if (strcmp(acids[i]->name, Hydrophilic[j]) == 0)
            {
                acids[i]->hydrophilic = true;
            }
        }
        if (acids[i]->hydrophilic == true)
        {
            printf("%s - hydrophilic\n", acids[i]->name);
        }
        else
        {
            printf("%s\n", acids[i]->name);
        }
        prev = acids[i];
    }


    printf("Using BOUND = %d\n", BOUND);

    



    thread_arg * ThreadData[BOUND][BOUND][BOUND];
    pthread_t threads[BOUND][BOUND][BOUND];

    for (int x = 0; x < BOUND; x++)
    {
        for (int y = 0; y < BOUND; y++)
        {
            for (int z = 0; z < BOUND; z++)
            {
                ThreadData[x][y][z] = calloc(sizeof(thread_arg), 1);
                thread_arg * arg = ThreadData[x][y][z];
                arg->x = x;
                arg->y = y;
                arg->z = z;
                arg->min_energy = 99999;
                memcpy(arg->pgrid, grid, sizeof(acid*)*DIM*DIM*DIM);
                memcpy(arg->acids_list, acids, sizeof(acid*)*NUM_ACIDS);
                int ret = pthread_create(&threads[x][y][z], NULL, thread_func, arg);
            }
        }
    }

    long total_recursive_calls = 0;
    long total_valid_configurations = 0;
    double true_min_energy = 9999;
    acid * true_optimal_grid[DIM][DIM][DIM];
    for (int x = 0; x < BOUND; x++)
    {
        for (int y = 0; y < BOUND; y++)
        {
            for (int z = 0; z < BOUND; z++)
            {
                pthread_join(threads[x][y][z], NULL);
                printf("Summary of thread %ld:\n", threads[x][y][z]);
                thread_arg * arg = ThreadData[x][y][z];
                printf("Recursive Calls: %ld, Valid Configurations: %ld, # New Mins: %d\n", 
                    arg->num_recursive_calls, arg->num_valid_configurations, arg->num_new_mins);
                printf("Minimum Energy: %f\n", arg->min_energy);
                //printf("Optimal Configuration:\n");
                //print_grid(arg->optimal_configuration);
                total_recursive_calls += arg->num_recursive_calls;
                total_valid_configurations += arg->num_valid_configurations;
                if (arg->min_energy < true_min_energy)
                {
                    memcpy(true_optimal_grid, arg->optimal_configuration, sizeof(acid*)*DIM*DIM*DIM);
                    true_min_energy = arg->min_energy;
                }
                
            }
        }
    }
    printf("\n\n");
    printf("Total Recursive Calls: %ld\n", total_recursive_calls);
    printf("Total valid configurations: %ld\n", total_valid_configurations);
    printf("True Minimum Energy: %f\n", true_min_energy);
    printf("True optimal state: \n");
    print_grid(true_optimal_grid);

    
    
    // //Write to PDB file
    // //write_PDB(optimal_configuration);
    // write_PDB_TWO(optimal_configuration, acids);
    save_grid(true_optimal_grid, acids);

    // //Cleanup
    for (int x = 0; x < DIM; x++)
    {
        for (int y = 0; y < DIM; y++)
        {
            for (int z = 0; z < DIM; z++)
            {
                free(Spatial_Database[x][y][z]);
            }
        }
    }
    for (int x = 0; x < BOUND; x++)
    {
        for (int y = 0; y < BOUND; y++)
        {
            for (int z = 0; z < BOUND; z++)
            {
                free(ThreadData[x][y][z]);
            }
        }
    }
    for (int i = 0; i < NUM_ACIDS;i++)
    {
        free(acids[i]);
    }

    

}

/**
 * Return Adjacent, C1 Diagonals, C2 Diagonals
 * Adjacent <==> coordinate is 1 unit from source (6)
 * C1 Diagonal <==> coordinate is sqrt(2) away from source (12)
 * C2 Diagonal <==> coordinate is sqrt(3) away from source (8)
*/
coordinate * get_positions(int x, int y, int z)
{
    coordinate * ret = calloc(sizeof(coordinate), 27);
    coordinate * tmp = ret;


    //Directly adjacent
    if (0 <= x - 1)
    {
        tmp->x = x - 1;
        tmp->y = y;
        tmp->z = z;
        tmp->valid = 1;
        tmp++;
    }
    if (x + 1 < DIM)
    {
        tmp->x = x + 1;
        tmp->y = y;
        tmp->z = z;
        tmp->valid = 1;
        tmp++;
    }

    if (0 <= y - 1)
    {
        tmp->x = x;
        tmp->y = y - 1;
        tmp->z = z;
        tmp->valid = 1;
        tmp++;
    }
    if (y + 1 < DIM)
    {
        tmp->x = x;
        tmp->y = y + 1;
        tmp->z = z;
        tmp->valid = 1;
        tmp++;
    }

    if (0 <= z - 1)
    {
        tmp->x = x;
        tmp->y = y;
        tmp->z = z - 1;
        tmp->valid = 1;
        tmp++;
    }
    if (z + 1 < DIM)
    {
        tmp->x = x;
        tmp->y = y;
        tmp->z = z + 1;
        tmp->valid = 1;
        tmp++;
    }

#   ifdef C1
    //C1 Diagonals in x = +/- 1 planes
    //Up to 8 coordinates
    
    int x_pm_1[2] = {x - 1, x + 1};
    int y_pm_1[2] = {y - 1, y + 1};
    int z_pm_1[2] = {z - 1, z + 1};
    for (int i = 0; i < 2; i++)
    {
        int currX = x_pm_1[i];
        if (0 <= currX && currX < DIM)
        {
            //y = +/- 1, z = z
            
            for (int j = 0; j < 2; j++)
            {
                int currY = y_pm_1[j];
                if (0 <= currY && currY < DIM)
                {
                    tmp->x = currX;
                    tmp->y = currY;
                    tmp->z = z;
                    tmp->valid = 1;
                    tmp->diagonal = 1;
                    tmp++;
                }
            }

            
            for (int j = 0; j < 2; j++)
            {
                int currZ = z_pm_1[j];
                if (0 <= currZ && currZ < DIM)
                {
                    tmp->x = currX;
                    tmp->y = y;
                    tmp->z = currZ;
                    tmp->valid = 1;
                    tmp->diagonal = 1;
                    tmp++;
                }
            }
        }
    }

    //C1 diagonals in x = x plane
    //Up to 4 coordinates

    for (int i = 0; i < 2; i++)
    {
        int currY = y_pm_1[i];
        if (0 <= currY && currY < DIM)
        {
            for (int j = 0; j < 2; j++)
            {
                int currZ = z_pm_1[j];
                if (0 <= currZ && currZ < DIM)
                {
                    tmp->x = x;
                    tmp->y = currY;
                    tmp->z = currZ;
                    tmp->valid = 1;
                    tmp->diagonal = 1;
                    tmp++;
                }
            }
        }
    }

#   endif
#   ifdef C2

    //C2 diagonals - up to 8 coordinates

    for (int i = 0; i < 2; i++)
    {
        int currX = x_pm_1[i];
        if (0 <= currX && currX < DIM)
        {
            for (int j = 0; j < 2; j++)
            {
                int currY = y_pm_1[j];
                if (0 <= currY && currY < DIM)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        int currZ = z_pm_1[k];
                        if (0 <= currZ && currZ < DIM)
                        {
                            tmp->x = currX;
                            tmp->y = currY;
                            tmp->z = currZ;
                            tmp->valid = 1;
                            tmp->diagonal = 2;
                            tmp++;
                        }
                    }
                }
            }
        }
    }

#   endif

    return ret;
}