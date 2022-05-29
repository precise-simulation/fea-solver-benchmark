% FEM Benchmark Matlab and Octave test run script.

% Copyright 2013-2022 Precise Simulation, Ltd.


clear all

path_to_featool = which('featool.m');
path_to_featool = 'C:\featool\featool.m';
if( ~isempty(path_to_featool) )
  path_to_featool = fileparts( path_to_featool );
else
  path_to_featool = input('Please enter the path to your FEATool toolbox installation: ','s')
end

if( exist('OCTAVE_VERSION','builtin') )
  warning( 'off', 'Octave:shadowed-function' )
end
addpath( genpath(path_to_featool) );
addpath( 'src_matlab' );


% n0    = 32
% nlev  = 6
% nruns = 20
params = load( 'testrun_param.txt' );
n0    = params(1);
nlev  = params(2);
nruns = params(3);


fem_poisson(n0);
timings = zeros(nlev,9);
for ilev = 1:nlev
  n = n0*2^(ilev-1);
  timings_ilev = zeros(1,8);

  for i = 1:max(3,nruns)
    [~,timings_i] = fem_poisson( n );
    if i ~= 1 && i ~= nruns
      timings_ilev = timings_ilev + timings_i;
    end
  end
  timings_ilev = timings_ilev/(nruns - 2);

  timings(ilev,1) = n0*2^(ilev-1);
  timings(ilev,2:end) = timings_ilev;
end


if( exist('OCTAVE_VERSION','builtin') )
  save( 'output/output_octave.txt', 'timings', '-ascii' )
else
  save( 'output/output_matlab.txt', 'timings', '-ascii' )
end


exit
