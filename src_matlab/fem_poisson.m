function [ u, timings ] = fem_poisson( nx, ny )
% Solves uxx + uyy = 1 on a unit square with Q1 bilinear finite elements.
% Returns the computed solution vector in u, and the timimngs vector where:
%
% t_grid   = timings[1] time for grid generation
% t_ptr    = timings[2] time to calculate matrix pointers
% t_asm    = timings[3] time for system matrix assembly
% t_rhs    = timings[4] time for right hand side assembly
% t_bdr    = timings[5] time to set boundary conditions
% t_sparse = timings[6] time to convert to sparse matrix format
% t_sol    = timings[7] time for solution
%
% FEM assembly functions from the FEATool Multiphysics Matlab FEM Toolbox are
% required. Available from https://www.featool.com/

% Copyright 2013-2018 Precise Simulation, Ltd.


tic
i_cub  = 2;
sfun   = 'sf_quad_Q1';
coef_a = [1 1];
coef_f = 1;

if( nargin<1 )
  nx = 16;
end
if( nargin<2 )
  ny = nx;
end
grid = rectgrid(nx,ny);
p    = grid.p;
c    = grid.c;
a    = [];
b    = grid.b;
t_grid = toc;
clear grid


tic
[vRowInds,vColInds,vAvals,n_rows,n_cols] = ...
assemblea( [2 3;2 3], {sfun;sfun}, coef_a, i_cub, p, c, a );
t_asm = toc;


tic
ind_c  = b(1,:);
i_edge = b(2,:);
j_edge = mod( i_edge, size(c,1) ) + 1;
ix = sub2ind( size(c), [i_edge j_edge], [ind_c ind_c] );
indrow = unique( c(ix) );

pp = sparse( max(vRowInds), 1 ) ;
pp(indrow) = 1;
itmp = 1:length(vRowInds);
indr = itmp(logical(pp(vRowInds)));  % Index to entries in matrix corresponding to boundary rows.

vAvals(indr) = 0;   % Zero out bc rows.
indr = indr(1:length(indrow));
vRowInds(indr) = indrow;
vColInds(indr) = indrow;
vAvals(indr) = 1;   % Set diagonals to 1.
t_bdr = toc;
clear b


tic
A = sparse( vRowInds, vColInds, vAvals, n_rows, n_cols );
t_sparse = toc;
clear vRowInds vColInds vAvals n_rows n_cols


tic
f = assemblef( 1, {sfun}, coef_f, i_cub, p, c, a );
t_rhs = toc;
clear p c a


tic
f(indrow) = 0;
t_bdr = t_bdr + toc;


tic
for i=1:100
  tmp = A*f;
end
t_sparse_mv = toc/100;


tic
u = A\f;
t_sol = toc;


T_PTR   = -1;
timings = [ t_grid, T_PTR, t_asm, t_rhs, t_bdr, t_sparse, t_sparse_mv, t_sol ];
