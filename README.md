This code was supposed to be private because it's so old that I don't recall many of its functionalities, but I thought it wouldn't hurt if if was public...

The main file is Script.m and all the other files are modularized because that's how MATLAB works.

Basically, it's possible to do two things with this code:

1) You can compute the best covering using only the MATLAB machinery by means of 'uniraioinexrest.m', which also includes verifying whether a given covering is or is not a valid one.
2) It is possible to use this code to generate the constraint functions to be fed to the ALGENCAN implementation of the augmented Lagrangian method (in FORTRAN). Then, ALGENCAN is able to compute the best covering to your input figure as well.

Item 1 was used in my master's thesis, whereas item 2 was used in the paper derived from it: https://link.springer.com/article/10.1007/s40314-019-0991-5

It's very important to mention that the figure to be covered need not be convex, but it has to be described in terms of a bunch of polynomial inequalities (i.e. it must be semialgebraic).

Should this be of anyone's interest, feel free to use my code as you wish and to contact me for whatever reason!
