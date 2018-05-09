function evol_nb_cell(state,time_flag)
% EVOL_NB_CELL(CELL_TYPE,TIME_FLAG)
%   Extract the number of cells of type CELL_TYPE (string) from ReSCAL logs
%   and plot the evolution.
%   If TIME_FLAG is set, read time data in t0 unit from TIME.log file.

if ~exist('time_flag','var')
  time_flag=0;
end

%% read header of CELL.log
fprintf('reading CELL.log\n');
fid=fopen('CELL.log');
if (fid==-1)
  fprintf('file not found\n');
  return
end
textscan(fid,'# CELL STATES');
C=textscan(fid,'NB_STATES = %d');
n=C{1};
CTypes=textscan(fid,'ST(%d): %s',n);
id_cell=-1;
for i=1:n
  if strcmp(CTypes{2}(i),state)
    celltype=CTypes{1}(i);
    id_cell=celltype+1;
    fprintf('id of state %s: %d\n', state, celltype);
  end
end

if (id_cell < 0)
  fprintf('state %s not found\n', state);
  return
end

%% skip line of text
C=textscan(fid,'%[^\n]s');
%celldisp(C)

%% read cell data
C=textscan(fid,'%d','delimiter',':');
data0=cell2mat(C);
l=size(data0,1);
nb=l/(n+1);
fprintf('number of records: %d\n', nb);
data=reshape(data0,n+1,nb)';
celldata=data(:,id_cell+1);
fclose(fid);

%% read time data
if time_flag>0
  fprintf('reading TIME.log\n');
  fid=fopen('TIME.log');
  if (fid==-1)
    fprintf('file not found\n');
    return
  end
  C=textscan(fid,'%[^\n]s');
  %celldisp(C)
  C=textscan(fid,'%f','delimiter',':');
  data0=cell2mat(C);
  l=size(data0,1);
  c=l/nb;
  data=reshape(data0,c,nb)';
  timedata=data(:,4);
  fclose(fid);
end

%% plot data
figure
if time_flag==0
  plot(celldata,'linewidth',2)
  xlabel('computation time (min.)')
else
  plot(timedata,celldata,'linewidth',2)
  xlabel('time (t_{0})')
end
title(sprintf('Evolution of the number of %s cells', state))
ylabel('number of cells')

return

