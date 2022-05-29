function plotly_chart( title, titles, data, filename, legends )
% Create and output plotly line plot diagrams.

% Copyright 2013-2018 Precise Simulation, Ltd.


plotlysrc = 'https://cdn.plot.ly/plotly-latest.min.js';

cols = [1 0 0;0 1 0;0 0 1;1 0 1;1 1 0;0 1 1];

if( ~iscell(data) )
  data = { data };
end

sm1 = @(s) s(1:end-1);

s = [ '<!DOCTYPE html><html><head><meta charset="UTF-8"><title>',title,'</title>', ...
      '<script src="',plotlysrc,'"></script></head><body>' ];

for i = 1:numel(data)
  s = [ s, '<div id="div',num2str(i),'" style="width:800px;height:600px;"></div>' ];
end

for i = 1:numel(data)
  d = data{i};
  s = [ s,'<script>'];
  x = d(:,1);
  sx = sprintf('%i,',d(:,1));

  for j = 2:size(d,2)
    s = [ s, plotlyvar( j-1, x, d(:,j), legends{j-1}, cols(mod(j-2,size(cols,1))+1,:) ) ];
  end

  s = [ s, 'var layout={title:''',titles{i},''',xaxis:{title:''1/h'',tickvals:[',sx(1:end-1),']},yaxis:{title:''cpu [s]''}};'];

  s = [ s, 'data=[',sm1(sprintf('tr%i,',1:size(d,2)-1)),'];Plotly.newPlot(''div',num2str(i),''', data, layout);</script>' ];
end

s = [ s, '</body></html>' ];

fid = fopen(filename,'w');
fprintf(fid,'%s',s);
fclose(fid);

try
  if( ispc )
    system( ['cmd.exe /c rundll32 url.dll, FileProtocolHandler "',filename,'"'] );
  elseif( isunix  )
    system( [system_dependent('getpref', 'HTMLSystemBrowserOptions'),' file://',filename] );
  elseif( ismac  )
    system( ['open ',filename] );
  end
catch
  disp(['Could not automatically open ',filename,' in browser.'])
end


function [ s ] = plotlyvar( i, x, y, name, color )

linewidth = 2;

sx = sprintf('%f,',x);
sy = sprintf('%f,',y);
s = ['var tr',num2str(i),'={', ...
     'x:[',sx(1:end-1),'],', ...
     'y:[',sy(1:end-1),'],', ...
     'type:''scatter'',', ...
     'name:''',name,''',', ...
...     'fill:''toself'',', ...
...     'fillcolor:''rgb(%i,%i,%i)'',', ...
     'line:{color:''',sprintf('rgb(%i,%i,%i)',round(255*(color))), ...
...     'marker:{opacity:0},', ...
     ''',width:',int2str(linewidth),'}', ...
     '};'];
