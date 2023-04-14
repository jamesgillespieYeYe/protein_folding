// #include <stdio.h>
// #include <string.h>
// #include <stdlib.h>
// #include <stdbool.h>
// #include <pthread.h>
#include "definitions.h"







//coordinate * get_positions(int x, int y, int z);
coordinate * Spatial_Database[DIM][DIM][DIM];
//void print_grid(acid * pgrid[DIM][DIM][DIM]);


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
 * Compute the total energy of a configuration by examining energies of 
 * all adjacencies
*/
double Total_Energy(acid * pgrid[DIM][DIM][DIM], acid * acids_list[NUM_ACIDS])
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
// #   ifndef C_1
// #       ifndef C_2
//     contact_map * map = gen_contact_map(pgrid, Spatial_Database);
//     double ret_energy = compute_energy_map(map, acids_list);
//     if (!(ret_energy - ret < .001 && ret - ret_energy < .001))
//     {
//         printf("DISAGREEMENT\n");
//         printf("%f, %f", ret, ret_energy);
//         exit(1);
//     }
//     free(map);
// #       endif
// #   endif
    return ret;
}







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


bool areDifferent(contact_map * one, contact_map * two)
{
    for (int i = 0; i < NUM_ACIDS; i++)
    {
        for (int j = 0; j < NUM_ACIDS; j++)
        {
            if (one->map[i][j] != two->map[i][j])
            {
                return true;
            }
        }
    }
    return false;
}

#define EPS .01
bool doubleEqual(double one, double two)
{
    if (one - two < EPS && two - one < EPS)
    {
        return true;
    }
    return false;
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
        double energy = Total_Energy(pgrid, acids_list);
        if (energy < args->min_energy)
        {
            args->min_energy = energy;
            memcpy(args->optimal_configuration, pgrid, sizeof(acid*)*DIM*DIM*DIM);
            args->num_new_mins++;
        }
#   ifdef TRACK //Track the top NUM_TRACK best and worst structures
        //BEST Structures
        bool needToSkip = false;
        contact_map * tmp = gen_contact_map(pgrid, Spatial_Database);
        for (int index = 0; index < NUM_TRACK; index++)
        {
            if (doubleEqual(energy, args->best_energies[index].energy) == true)
            {
                //Need to make sure this is a unique entry by
                //checking if the contact maps are the same
                //If they aren't the same, we are ok, but if they are,
                //we break out of the loop
                if (args->best_energies[index].map != NULL)
                {
                    if (!areDifferent(tmp, args->best_energies[index].map))
                    {
                        needToSkip = true;
                        break;
                    }
                }
            }
        }
        free(tmp);
        if (!needToSkip)
        {
            for (int index = 0; index < NUM_TRACK; index++)
            {
                if (energy < args->best_energies[index].energy)
                {
                    //This structure is unique AND is better than the current structure
                    //So we want want to copy this data to this location,
                    //and shift the others down
                    entry tmp;
                    memcpy((void*) &tmp, (void*) &(args->best_energies[index]), sizeof(entry));
                    args->best_energies[index].energy = energy;
                    args->best_energies[index].map = gen_contact_map(pgrid, Spatial_Database);
                    memcpy(args->best_energies[index].grid, pgrid, sizeof(acid*)*DIM*DIM*DIM);
                    for (int sub_index = index + 1; sub_index < NUM_TRACK; sub_index++)
                    {
                        entry storage;
                        memcpy((void*) &storage, (void*) &(args->best_energies[sub_index]), sizeof(entry));
                        memcpy((void*) &(args->best_energies[sub_index]), (void*) &tmp, sizeof(entry));
                        memcpy((void*) &tmp, (void*) &storage, sizeof(entry));
                    }
                    //Free the dynamic memory associated with the last entry that got 
                    //"popped" off the list
                    free(tmp.map);


                    break;
                }
            }
        }
        //WORST Structures
        needToSkip = false;
        tmp = gen_contact_map(pgrid, Spatial_Database);
        for (int index = 0; index < NUM_TRACK; index++)
        {
            if (doubleEqual(energy, args->worst_energies[index].energy) == true)
            {
                //Need to make sure this is a unique entry by
                //checking if the contact maps are the same
                //If they aren't the same, we are ok, but if they are,
                //we break out of the loop
                if (args->worst_energies[index].map != NULL)
                {
                    if (!areDifferent(tmp, args->worst_energies[index].map))
                    {
                        needToSkip = true;
                        break;
                    }
                }
                
            }
        }
        free(tmp);
        if (!needToSkip)
        {
            for (int index = 0; index < NUM_TRACK; index++)
            {
                if (energy > args->worst_energies[index].energy)
                {
                    //This structure is unique AND is worse than the current structure
                    //So we want want to copy this data to this location,
                    //and shift the others down
                    entry tmp;
                    memcpy((void*) &tmp, (void*) &(args->worst_energies[index]), sizeof(entry));
                    args->worst_energies[index].energy = energy;
                    args->worst_energies[index].map = gen_contact_map(pgrid, Spatial_Database);
                    memcpy(args->worst_energies[index].grid, pgrid, sizeof(acid*)*DIM*DIM*DIM);
                    for (int sub_index = index + 1; sub_index < NUM_TRACK; sub_index++)
                    {
                        entry storage;
                        memcpy((void*) &storage, (void*) &(args->worst_energies[sub_index]), sizeof(entry));
                        memcpy((void*) &(args->worst_energies[sub_index]), (void*) &tmp, sizeof(entry));
                        memcpy((void*) &tmp, (void*) &storage, sizeof(entry));
                    }
                    //Free the dynamic memory associated with the last entry that got 
                    //"popped" off the list
                    free(tmp.map);


                    break;
                }
            }
        }
        
#   endif
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
    printf("New thread: %ld: Initial Position: (%d, %d, %d)\n", pthread_self(), arg->x, arg->y, arg->z);
    arg->pgrid[arg->x][arg->y][arg->z] = arg->acids_list[0];
    insert(arg->pgrid, arg->acids_list, 1, arg->x, arg->y, arg->z, arg);
    return NULL;
}


void print_acids(acid * first)
{
    acid * curr = first;
    while (curr != NULL)
    {
        printf("%s--", curr->name);
        curr = curr->next;
    }
    printf("\n");
    curr = first;
    while (curr != NULL)
    {
        if (curr->hydrophilic == true)
        {
            printf("P----");
        }
        else
        {
            printf("H----");
        }
        curr = curr->next;
    }
    printf("\n");
}
void print_acids_backwards(acid * last)
{
    acid * curr = last;
    while (curr != NULL)
    {
        printf("%s--", curr->name);
        curr = curr->prev;
    }
    printf("\n");
    curr = last;
    while (curr != NULL)
    {
        if (curr->hydrophilic == true)
        {
            printf("P----");
        }
        else
        {
            printf("H----");
        }
        curr = curr->prev;
    }
    printf("\n");
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
        acids[i]->index = i;
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

    print_acids(acids[0]);
    print_acids_backwards(acids[NUM_ACIDS - 1]);
    
    printf("Using %d as bound on starting position\n", BOUND);
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
#               ifdef TRACK
                for (int index = 0; index < NUM_TRACK; index++)
                {
                    arg->worst_energies[index].energy = -999;
                }
#               endif
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
    printf("Total Recursive Calls: %ld\n", total_recursive_calls);
    printf("Total valid configurations: %ld\n", total_valid_configurations);
    printf("True Minimum Energy: %f\n", true_min_energy);
    printf("True optimal state: \n");
    print_grid(true_optimal_grid);
#   ifdef TRACK
    for (int x = 0; x < BOUND; x++)
    {
        for (int y = 0; y < BOUND; y++)
        {
            for (int z = 0; z < BOUND; z++)
            {
                printf("(%d, %d, %d)   Best   Worst\n", ThreadData[x][y][z]->x, ThreadData[x][y][z]->y, ThreadData[x][y][z]->z);
                for (int index = 0; index < NUM_TRACK; index++)
                {
                    entry currBest = ThreadData[x][y][z]->best_energies[index];
                    entry currWorst = ThreadData[x][y][z]->worst_energies[index];
                    printf("Energy: %f %f\n", currBest.energy, currWorst.energy);
                }
                for (int i = 0; i < NUM_TRACK; i++)
                {
                    for (int j = 0; j < NUM_TRACK; j++)
                    {
                        if (i != j)
                        {
                            if (!areDifferent(ThreadData[x][y][z]->best_energies[i].map, 
                                ThreadData[x][y][z]->best_energies[j].map))
                            {
                                printf("Error with best structures: %d and %d\n", i, j);
                            }
                        }
                    }
                }
                for (int i = 0; i < NUM_TRACK; i++)
                {
                    for (int j = 0; j < NUM_TRACK; j++)
                    {
                        if (i != j)
                        {
                            if (!areDifferent(ThreadData[x][y][z]->worst_energies[i].map, 
                                ThreadData[x][y][z]->worst_energies[j].map))
                            {
                                printf("Error with WORST structures: %d and %d\n", i, j);
                            }
                        }
                    }
                }
            }
        }
    }
#   endif

    
    
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
#               ifdef TRACK
                for (int index = 0; index < NUM_TRACK; index++)
                {
                    free(ThreadData[x][y][z]->best_energies[index].map);
                }
#               endif
                free(ThreadData[x][y][z]);
            }
        }
    }
    for (int i = 0; i < NUM_ACIDS;i++)
    {
        free(acids[i]);
    }

    

}

