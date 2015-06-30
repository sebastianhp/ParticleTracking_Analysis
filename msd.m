% MSD: calculates mean square displacement from an ensemble of trajectory 
% data of single particles.
%
% List of inputs:
%
% * tracksCell (cell array): 1xnumTracks cell array containing a single
% * trajectory in each cell
% * divs (int): number of times each trajectory will be divided
% * divsize (int): how many timeframes for each division
% * dim (int): dimension along which to take msd. 1 is along x, 2 along y
% 
% The program expects the tracked files to be in the format:
% Column 1: Position in x
% Column 2: Position in y
% Column 3: Timestamp
% Each row its own timestep


function [times, msd] = msd(tracksCell, divs, divsize, dim)

marker = false;
sz = size(tracksCell);
numTracks = sz(2);

if divs == -1
    divs = 1;
end

splitTracksCell = cell(numTracks,divs+1); % Rows are different runs
                                        % Columns are different divisions
                                        % The last column is the time step
                                        % for all the divisions in the run.
                                        % Therefore we assume that the
                                        % timesteps are all the same.
sizes = zeros(1,numTracks);

for i=1:numTracks
    if divsize == -1
        s = size(tracksCell{i});
        sizes(i) = s(1);
        divsize = s(1);
        marker = true;
    end
    trackedRun = tracksCell{i};
    for j=1:divs
        temp = trackedRun( (1+(j-1)*divsize):j*divsize,dim) - trackedRun(1,dim);
        splitTracksCell{i,j} = temp;
    end
    splitTracksCell{i,divs+1} = trackedRun(1:divsize,3) - trackedRun(1,3); % time column
    if marker
        divsize = -1;
    end
end

if marker
    divsize = min(sizes);
end

msd = zeros(divsize,1);

for i=1:divsize
    xvals = zeros(numTracks,divs);
    deltaXsq = zeros(numTracks,divs);
    for j=1:numTracks
        for k=1:divs
            temp = splitTracksCell{j,k};
            xvals(j,k) = temp(i,dim);
        end
    end
    avgx = mean(mean(xvals));
    
    for j=1:numTracks
        for k=1:divs
            deltaXsq(j,k) = (xvals(j,k) - avgx)^2;
        end
    end
    msd(i,1) = mean(mean(deltaXsq));
end

temp = splitTracksCell{1,divs+1};
times = temp(1:divsize);

end
        
