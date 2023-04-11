#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

//Options
#define C1
#define C2
#define WATER_PENALTY
#define DIM 3

#define WP 2.5


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




long num_recursive_calls;


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
        return -2;
    }
    else if (a->hydrophilic == true && b->hydrophilic == false)
    {
        return -1;
    }
    else if (a->hydrophilic == false && b->hydrophilic == true)
    {
        return -1;
    }
    else 
    {
        return -0.1;
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




long num_valid_configurations = 0;
double min_energy = 999;
long num_new_mins = 0;
acid * optimal_configuration[DIM][DIM][DIM];

/**
 * Recursive callback
 * Look at all positions adjacent to (lastX, lastY, lastZ), and try placing next acid there if possible
*/
void insert(acid * pgrid[DIM][DIM][DIM], acid * acids_list[NUM_ACIDS], int index, int lastX, int lastY, int lastZ)
{
    num_recursive_calls++;
    if (index == NUM_ACIDS)
    {
        num_valid_configurations++;
        if (num_valid_configurations % 10000000 == 0)
        {
            printf("Base case number: %ld\n", num_valid_configurations);
        }
        double energy = Total_Energy(pgrid);
        if (energy < min_energy)
        {
            min_energy = energy;
            memcpy(optimal_configuration, pgrid, sizeof(acid*)*DIM*DIM*DIM);
            num_new_mins++;
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

            insert(next_grid, acids_list, index + 1, pmove->x, pmove->y, pmove->z);
        }
        pmove++;
    }
}

#define GRID_FILENAME "grid.out.txt"
void save_grid(acid * pgrid[DIM][DIM][DIM], acid * acids[NUM_ACIDS])
{
    FILE * f = fopen(GRID_FILENAME, "w");

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
                        fprintf(f, "%s,%d,%d,%d\n", curr_acid->name, x, y, z);
                    } 
                }
            }
        }
    }
    fclose(f);
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


    int bound = DIM / 2;
    if (bound < (DIM - bound))
    {
        bound = DIM - bound;
    }
    printf("Using bound of %d\n", bound);


    for (int x = 0; x < bound; x++)
    {
        for (int y = 0; y < bound; y++)
        {
            for (int z = 0; z < bound; z++)
            {
                printf("Starting point: (%d, %d, %d)\n", x, y, z);
                acid * tmp_grid[DIM][DIM][DIM];
                memcpy(tmp_grid, grid, sizeof(acid*)*DIM*DIM*DIM);
                tmp_grid[x][y][z] = acids[0];
                insert(tmp_grid, acids, 1, x, y, z);
            }
        }
    }

    printf("Num Recursive Calls: %ld\n", num_recursive_calls);
    printf("Valid configurations tested: %ld\n", num_valid_configurations);
    printf("Num new mins: %ld\n", num_new_mins);
    printf("Min energy: %f\n", min_energy);
    printf("Optimal configuration:\n");
    print_grid(optimal_configuration);

    //Write to PDB file
    //write_PDB(optimal_configuration);
    write_PDB_TWO(optimal_configuration, acids);
    save_grid(optimal_configuration, acids);


    //Cleanup
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