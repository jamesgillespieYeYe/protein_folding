
This program the follows procedure outlined in a 1996 article
titled: 
"Emergence of Preferred Structures in a Simple Model of Protein Folding" 
authored by Hao Li, Robert Helling, Chao Tang and Ned Wingreen,
and published by the American Association for the Advancement of Science.
It can be found at:
https://www.jstor.org/stable/2891163
 
The authors sought to analyze the protein folding problem using
a simplified model, representings proteins as "beads (molecules) 
on a chain". The "simplified protein folding" problem is to arrange this
chain of beads in a 2D or 3D lattice in such a way that the energy of
the arrangment is minimized. "Energy" is computed by looking at adjacent
pairs of monomers in the lattice, and an assigning an energy to their 
contact based on their classification as hydrophobic or hydrophilic. 
 
The authors of the paper carried out a full analysis of structures, 
sequences that design them, and the relationship between the models 
they generated and structures in the real world. They enumerated all 
possible structures a 27-monomer protein could take in a 3x3x3 grid, 
and computed the ideal structure for all 2^27 possible configurations of 
hydrophobic-hydrophilic sequences. They carried out a very interesting
and informatative analysis of their results. This is NOT what I am doing.
 
I am simply trying to analyze one structure - the "mini t" protein, which
I am investigating for a course. I would like to obtain an approximation of 
the characteristics this protein might have. The real result will not look
much like the 3D lattices we generate here, but they might provide some 
insight into the general form the structure will take. In particular, 
this experiment might illuminate which monomers are likely to contact each 
other, and what what substructures might appear. 
 
My approach is simple - compute all possible arrangments of the protein, 
compute the energy for each one, and report the optimal arrangement. This is achieved 
through a multithreaded approach. 8 threads are spawned, each with an associated
starting point to place the first acid. The algorithm works recursively - 
on each iterations, it calls itself for each possible position to place the next 
monomer.
 
The design philosophy is also simple - speed and simplicity! Data structures 
are generated before the algorithm runs that significantly reduce complexity.
Dynamic memory usage is minimized. Values are often precomputed. And of course, 
we make good use the compiler :). Any computation that only occurs once can be 
as ugly as inefficient as I want; anything called millions of times, of course, 
must be fast. 
 
There are several modes the program can be run in - they are defined in the header
file. See the writeup for an explanation. A secondary program plots results.
 
compile with gcc -o proteins -ggdb3 -pthread protein3D.c 
 

