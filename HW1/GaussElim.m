% This function solves the equation AX = B using Gaussian elimination with 
% partial pivoting.The program also will calculate the determinant from the
% augmented matrix and number of pivots.

function [X, determ] = GaussElim(A, B)