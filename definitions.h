#ifndef DEFINITIONS_H
#define DEFINITIONS_H
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
#ifndef C1
#define TRACK
#endif
#define NUM_TRACK 8

//Bound defines the range of starting positions
// #if (DIM / 2) < (DIM - DIM / 2)
// #   define BOUND DIM - DIM / 2
// #else
// #   define BOUND DIM / 2
// #endif

#define BOUND 2

//WP is the Water Penalty Amount
#define WP 2.5

//Values for contacts
// #define Hhh -2
// #define Hhp -1
// #define Hpp -0.1

#define Hhh -2.3
#define Hhp -1
#define Hpp 0



#define NUM_ACIDS 20
//#define NUM_ACIDS 27
#define NUM_HYDROPHILIC 7

struct acid 
{
    char name[4];
    bool hydrophilic;
    struct acid * prev;
    struct acid * next;
    int index;
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

struct contact_map
{
    int map[NUM_ACIDS][NUM_ACIDS];
};
typedef struct contact_map contact_map;

struct entry
{
    struct contact_map * map;
    double energy;
    acid * grid[DIM][DIM][DIM];
};
typedef struct entry entry;


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
    entry best_energies[NUM_TRACK];
    entry worst_energies[NUM_TRACK];
};
typedef struct thread_arg thread_arg;



char * Sequence[NUM_ACIDS] = {"ASN", "LEU", "TYR", "ILE", "GLN", "TRP", "LEU", "LYS", "ASP", "GLY", "GLY", "PRO", "SER", "SER", "GLY", "ARG", "PRO", "PRO", "PRO", "SER"};
//char * Sequence[NUM_ACIDS] = {"ASN", "LEU", "TYR", "ILE", "GLN", "TRP", "LEU", "LYS", "ASP", "GLY", "GLY", "PRO", "SER", "SER", "GLY", "ARG", "PRO", "PRO", "PRO", "SER","SER", "GLY", "ARG", "PRO", "PRO", "PRO", "SER"};
char * Hydrophilic[NUM_HYDROPHILIC] = {"ASP", "ARG", "LYS", "GLY", "ASN", "SER", "GLN"};

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
 * Give the energy between the acid pair (a, b)
*/
double Energy(acid * a, acid * b)
{
    if (a->prev == b || a->next == b)
    {
        return 0;
    }
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

double compute_energy_map(contact_map *map, acid * acids[NUM_ACIDS]);
void print_contact_map(contact_map * map);
contact_map * gen_contact_map(acid *pgrid[DIM][DIM][DIM], coordinate * spatial_database[DIM][DIM][DIM]);

double compute_energy_map(contact_map *map, acid * acids[NUM_ACIDS])
{
    double ret = 0;
    for (int i = 0; i < NUM_ACIDS; i++)
    {
        for (int j = 0; j < NUM_ACIDS; j++)
        {
            if (map->map[i][j] == 1)
            {
                ret += Energy(acids[i], acids[j]);
            }
        }
    }
    return ret / 2;
}

void print_contact_map(contact_map * map)
{
    for (int i = 0; i < NUM_ACIDS; i++)
    {
        printf("--");
    }
    printf("\n\n");
    for (int i = 0; i < NUM_ACIDS; i++)
    {
        for (int j = 0; j < NUM_ACIDS; j++)
        {
            printf("%d ", map->map[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (int i = 0; i < NUM_ACIDS; i++)
    {
        printf("--");
    }
    printf("\n\n");
}
contact_map * gen_contact_map(acid *pgrid[DIM][DIM][DIM], coordinate * spatial_database[DIM][DIM][DIM])
{
    contact_map * ret = calloc(sizeof(contact_map), 1);
    for (int x = 0; x < DIM; x++)
    {
        for (int y = 0; y < DIM; y++)
        {
            for (int z = 0; z < DIM; z++)
            {
                acid * curr = pgrid[x][y][z];
                if (curr != NULL)
                {
                    //printf("name: %s\n", curr->name);
                    coordinate * adjacent = spatial_database[x][y][z];
                    coordinate *tmp = adjacent;
                    while (tmp->valid == 1 && tmp->diagonal == 0)
                    {
                        acid * neighbor = pgrid[tmp->x][tmp->y][tmp->z];
                        if (neighbor != NULL)
                        {
                            if (curr->next != neighbor && curr->prev != neighbor)
                            {
                                ret->map[curr->index][neighbor->index] = 1;
                            }
                            
                        }
                        tmp++;
                    }
                }
            }
        }
    }
    return ret;
}

#endif