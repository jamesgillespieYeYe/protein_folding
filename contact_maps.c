#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>

#define NUM_ACIDS 5


void print_map(int map[NUM_ACIDS][NUM_ACIDS])
{
    for (int x = 0; x < NUM_ACIDS; x++)
    {
        for (int y = 0; y < NUM_ACIDS; y++)
        {
            printf("%d ", map[x][y]);
        }
        printf("\n");
    }
}


long callback(int pmap[NUM_ACIDS][NUM_ACIDS], int row, int index, int count)
{
    int ret = 0;
    int map[NUM_ACIDS][NUM_ACIDS];
    memcpy(map, pmap, sizeof(int)*NUM_ACIDS*NUM_ACIDS);
    if (index < row)
    {
        index = row + 1;
    }
    
    if (count == 4 || index > NUM_ACIDS - 1)
    {
        if (row == NUM_ACIDS - 1)
        {
            return 1;
        }
        else
        {
            int newRow = row + 1;
            int newCount = 0;
            for (int y = 0; y < NUM_ACIDS; y++)
            {
                newCount += map[newRow][y];
            }
            return callback(map, row + 1, 0, newCount);
        }
    }
    if (pmap[row][index] == 1)
    {
        //This position is already filled; can't change it now
        return callback(map, row, index + 1, count);
    }
    if (row == index || row == index - 1 || row == index + 1)
    {
        //We can't be adjacent to ourself, our predecessor, or our successor
        return callback(map, row, index + 1, count);
    }
    //You can leave this slot empty 
    ret += callback(map, row, index + 1, count);
    //Or you can place something here
    map[row][index] = 1;
    map[index][row] = 1;
    ret += callback(map, row, index + 1, count + 1);

}

int main(int argc, char** argv)
{
    int map[NUM_ACIDS][NUM_ACIDS];
    memset((void*) map, 0, sizeof(int)*NUM_ACIDS*NUM_ACIDS);
    int ret = callback(map, 0, 0, 0);
    printf("%d\n", ret);
    print_map(map);
    
    
    
}