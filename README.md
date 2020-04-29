# MF_Moments

This Repository contains MATLAB code to calculate the central moments of a matrix Fisher distribution defined on the special orthogonal group SO(3).
The **num** folder contains functions for numerical calculations, and the **sym** folder contains functions for symbolic calculations.

## Guideline
The major function is
```
output = EQ(s,ind)
```
* The input **s** is the concentration parameter for a Matrix Fisher distribution, which is a 1-by-3 or a 3-by-1 vector;
* The input **ind** is a n-by-2 matrix of the form [i1,j1;...;in,jn];
* The **output** is the expectation of Q(i1,j1)*...*Q(in,jn) with Q following a matrix Fisher distribution with parameter S=diag(s).
