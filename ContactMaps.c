#ifndef DEFINITIONS_H
#include "definitions.h"
#endif




int load_grid(acid * pgrid[DIM][DIM][DIM], acid * acids[NUM_ACIDS], char * filename)
{
    FILE * f = fopen(filename, "r");
    char buf[200];
    int acid_count = 0;
    fgets(&buf[0], 200, f);
    while (fgets(&buf[0], 200, f))
    {
        //printf("%s\n", buf);
        int found[3] = {-1, -1, -1};
        int index = 0;
        char * ptr = &buf[0];
        bool cont = true;
        while (cont)
        {
            if (*ptr == ',')
            {
                ptr++;
                found[index] = *ptr - '0';
                index++;
            }
            if (index > 2)
            {
                cont = false;
            }
            ptr++;
        }
        //printf("%d, %d, %d\n", found[0], found[1], found[2]);
        pgrid[found[0]][found[1]][found[2]] = acids[acid_count];
        acid_count++;
    }
    fclose(f);
}


#ifdef MAPS_MAIN
acid * acids[NUM_ACIDS];
coordinate * Spatial_Database[DIM][DIM][DIM];
int main(int argc, char** argv)
{
    printf("hi\n");
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

    acid * grid[DIM][DIM][DIM];
    memset((void*) grid, 0, sizeof(acid*)*DIM*DIM*DIM);
    load_grid(grid, acids, "grid.out.txt");
    print_grid(grid);
    contact_map * map = gen_contact_map(grid, Spatial_Database);
    print_contact_map(map);
    double energy = compute_energy_map(map, acids);
    printf("%f\n", energy);

}

#endif

