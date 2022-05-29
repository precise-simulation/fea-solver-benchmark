% FEM Benchmark result processing script.

% Copyright 2013-2018 Precise Simulation, Ltd.


clear all
close all
clc


lw = 2;
mz = 15;


outdir = 'output';
addpath(outdir)
if( ~isempty(lastwarn) )
  outdir = ['..',filesep,'output'];
  addpath(outdir)
end
res_o  = load( [outdir,filesep,'output_octave.txt'] );
res_m  = load( [outdir,filesep,'output_matlab.txt'] );
res_jl = load( [outdir,filesep,'output_julia.txt'] );
res_f  = zeros(size(res_o));
for i=1:numel(res_m(:,1))
  res_f_i = load( [outdir,filesep,'results_mgfeat_nx=',num2str(res_m(i,1)),'.log'] );
  res_f(i,:) = sum( res_f_i(2:end-1,:), 1 )/(size(res_f_i,1)-2);
end


figure, hold on
title('Matrix Assembly')
plot( res_o(:,1),  res_o(:,4),  '.r-', 'linewidth', lw, 'markersize', mz )
plot( res_m(:,1),  res_m(:,4),  '.g-', 'linewidth', lw, 'markersize', mz )
plot( res_jl(:,1), res_jl(:,4), '.b-', 'linewidth', lw, 'markersize', mz )
plot( res_o(:,1),  res_f(:,4),  '.m-', 'linewidth', lw, 'markersize', mz )
set( gca, 'xtick', res_m(:, 1) )
ylabel('cpu [s]')
xlabel('1/h')
grid on
legend('Octave','Matlab','Julia','Fortran','location','northwest')
print('-r300','-djpeg',[outdir,filesep,'fig_assembly.jpg'])
print('-r300','-depsc',[outdir,filesep,'fig_assembly.eps'])


figure, hold on
title('Total time (without solver)')
plot( res_o(:,1),  sum(res_o(:,2:7),2),  '.r-', 'linewidth', lw, 'markersize', mz )
plot( res_m(:,1),  sum(res_m(:,2:7),2),  '.g-', 'linewidth', lw, 'markersize', mz )
plot( res_jl(:,1), sum(res_jl(:,2:7),2), '.b-', 'linewidth', lw, 'markersize', mz )
plot( res_o(:,1),  sum(res_f(:,2:7),2),  '.m-', 'linewidth', lw, 'markersize', mz )
set( gca, 'xtick', res_m(:, 1) )
ylabel('cpu [s]')
xlabel('1/h')
grid on
legend('Octave','Matlab','Julia','Fortran','location','northwest')
print('-r300','-djpeg',[outdir,filesep,'fig_total.jpg'])
print('-r300','-depsc',[outdir,filesep,'fig_total.eps'])


figure, hold on
title('Solver')
plot( res_o(:,1),  res_o(:,9),  '.r-', 'linewidth', lw, 'markersize', mz )
plot( res_m(:,1),  res_m(:,9),  '.g-', 'linewidth', lw, 'markersize', mz )
plot( res_jl(:,1), res_jl(:,9), '.b-', 'linewidth', lw, 'markersize', mz )
plot( res_o(:,1),  res_f(:,9),  '.m-', 'linewidth', lw, 'markersize', mz )
set( gca, 'xtick', res_m(:, 1) )
ylabel('cpu [s]')
xlabel('1/h')
grid on
legend('Octave','Matlab','Julia','Fortran','location','northwest')
print('-r300','-djpeg',[outdir,filesep,'fig_solve.jpg'])
print('-r300','-depsc',[outdir,filesep,'fig_solve.eps'])






fid = fopen([outdir,filesep,'tables.txt'],'w');
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
fprintf(fid,['|  Octave |',repmat('           |',1,9),'\n']);
fprintf(fid,'|     1/h | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', '   t_grid', '    t_ptr', '  t_asm_A', '  t_asm_f', '    t_bdr', ' t_sparse', '    t_tot', '   t_spmv', '  t_solve' );
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
for i=1:size(res_jl,1)
  fprintf(fid,'|   %5i | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e |\n', [res_o(i,1:7), sum(res_o(i,2:7)), res_o(i,8:9)] );
end
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);

fprintf(fid,'\n');

fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
fprintf(fid,['|  Matlab |',repmat('           |',1,9),'\n']);
fprintf(fid,'|     1/h | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', '   t_grid', '    t_ptr', '  t_asm_A', '  t_asm_f', '    t_bdr', ' t_sparse', '    t_tot', '   t_spmv', '  t_solve' );
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
for i=1:size(res_jl,1)
  fprintf(fid,'|   %5i | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e |\n', [res_m(i,1:7), sum(res_m(i,2:7)), res_m(i,8:9)] );
end
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);

fprintf(fid,'\n');

fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
fprintf(fid,['|   Julia |',repmat('           |',1,9),'\n']);
fprintf(fid,'|     1/h | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', '   t_grid', '    t_ptr', '  t_asm_A', '  t_asm_f', '    t_bdr', ' t_sparse', '    t_tot', '   t_spmv', '  t_solve' );
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
for i=1:size(res_jl,1)
  fprintf(fid,'|   %5i | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e |\n', [res_jl(i,1:7), sum(res_jl(i,2:7)), res_jl(i,8:9)] );
end
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);

fprintf(fid,'\n');

fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
fprintf(fid,['| Fortran |',repmat('           |',1,9),'\n']);
fprintf(fid,'|     1/h | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', '   t_grid', '    t_ptr', '  t_asm_A', '  t_asm_f', '    t_bdr', ' t_sparse', '    t_tot', '   t_spmv', '  t_solve' );
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
for i=1:size(res_f,1)
  fprintf(fid,'|   %5i | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e | %8.2e |\n', [res_f(i,1:7), sum(res_f(i,2:7)), res_f(i,8:9)] );
end
fprintf(fid,['|---------+',repmat('-----------+',1,8),'-----------|\n']);
fclose(fid);


% Plotly
data = { [res_o(:,1) res_o(:,2) res_m(:,2) res_jl(:,2) res_f(:,2)] ;
         [res_o(:,1) res_o(:,3) res_m(:,3) res_jl(:,3) res_f(:,3)] ;
         [res_o(:,1) res_o(:,4) res_m(:,4) res_jl(:,4) res_f(:,4)] ;
         [res_o(:,1) res_o(:,5) res_m(:,5) res_jl(:,5) res_f(:,5)] ;
         [res_o(:,1) res_o(:,6) res_m(:,6) res_jl(:,6) res_f(:,6)] ;
         [res_o(:,1) res_o(:,7) res_m(:,7) res_jl(:,7) res_f(:,7)] ;
         [res_o(:,1) sum(res_o(:,2:7),2) sum(res_m(:,2:7),2) sum(res_jl(:,2:7),2) sum(res_f(:,2:7),2)] ;
         [res_o(:,1) res_o(:,8) res_m(:,8) res_jl(:,8) res_f(:,8)] ;
         [res_o(:,1) res_o(:,9) res_m(:,9) res_jl(:,9) res_f(:,9)] };
title = 'Results and Timings';
filename = [outdir,filesep,'results.html'];
legends = {'Octave','Matlab','Julia','Fortran'};
titles = { 'Grid' 'Matrix Pointers' 'Matrix Assembly' 'RHS Assembly' 'Boundary Conditions' 'Sparse Allocation' 'Total Time (without solver)', 'Sparse MV', 'Solver' };
plotly_chart( title, titles, data, filename, legends )


% exit
