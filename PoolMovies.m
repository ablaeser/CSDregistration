function [regPooledMovie, regPooledShift] = PoolMovies( movieDirs, fovName, fovDir, varargin )
%Concatenate multiple registered movies from the same field of view into
%one movie and register all frames to the 

IP = inputParser;
addRequired( IP, 'movieDirs', @iscell )
addRequired( IP, 'fovName', @ischar )
addRequired( IP, 'fovDir', @ischar )
addParameter( IP, 'regColor', 1, @isnumeric ) % regColor = 1; %2;
%addParameter( IP, 'bin', 1, @isnumeric )
parse( IP, movieDirs, fovName, fovDir, varargin{:} );
regColor = IP.Results.regColor;
%binSize = IP.Results.bin; %binSize = 10;
RGBOpt = struct('overwrite',true, 'message',true, 'append',false, 'big',true, 'color',true );
GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
mkdir(fovDir)
% Load the individually-registered movies and their metadata
Nmovie = numel(movieDirs);
fprintf('\n\nPooling data: %s  (%i movies)', fovName, Nmovie );
for m = flip(1:Nmovie)
    [metadata(m), T{m}, otherData(m), regMovie{m}] = LoadProcessed( movieDirs{m}, 'reg' );
    shift{m} = otherData(m).shift; % shift = (Nframe x 5 x 3)array of [ xshift, yshift, shift distance, error, phase difference ] 
    %regParam(m) = otherData(m).regParam;
    [~,movieName{m}] = fileparts( metadata(m).root );
    planeProj{m} = metadata(m).mean;
    CSDframes{m} = otherData(m).CSD;
end

% Show mean projections
figure('WindowState','max');
opt = {[0.07,0.01], [0.06,0.04], [0.01,0.01]};  % {[vert, horz], [bottom, top], [left, right] } 
for m = 1:Nmovie
    subtightplot(1,Nmovie,m,opt{:}); imshow( metadata(m).mean(:,:,2), [] );
    title( sprintf('m = %i: %s', m, movieName{m} ), 'Interpreter','none' ); 
end
% Save mean projections to FOV folder as a single RGB tiff stack
PlaneMeanPath = sprintf('%s%s_PlaneMean.tif', fovDir, fovName );
if ~exist(PlaneMeanPath, 'file')  
    fovPlaneProj = cat(4, planeProj{:});
    saveastiff( fovPlaneProj, PlaneMeanPath, RGBOpt );  
end


% Calculate the error and shifts between each pair of projections, and identify the most typical movie
projFFT = cell(1,Nmovie);  errorMat = nan(Nmovie,Nmovie);  shiftMat = nan(Nmovie,Nmovie,4); % [distance, x shift, y shift, phase diff]
for m = 1:Nmovie,  projFFT{m} = fft2( metadata(m).mean(:,:,2) );  end

Nframe = zeros(1,Nmovie);
for m = 1:Nmovie
    Nframe(m) = size(regMovie{m},4);
    for n = 1:Nmovie
        [output, ~] = dftregistration( projFFT{m}, projFFT{n}, 10 ); % output = [error, diffphase, net_row_shift, net_col_shift]
        errorMat(m,n) = output(1);
        shiftMat(m,n,:) = [norm(output(3:4)), output(4), output(3), output(2)];
        %regMovie(:,:,c,z) = abs( ifft2(fftIndReg) );
    end
end
[~, mMin] = min( sum(errorMat,1) ); % which movie has the least total error
firstFrame = [0,cumsum(Nframe)] + 1;
movieFrames = cell(1,Nmovie); %CSDpoolFrames = [];
for m = 1:Nmovie
    movieFrames{m} = firstFrame(m):firstFrame(m+1)-1;
    %CSDpoolFrames = [CSDpoolFrames, CSDframes{m} + firstFrame(m)-1];
end
%goodPoolFrames = 1:firstFrame(end)-1;  goodPoolFrames(CSDpoolFrames) = []; 
%maxIndvShift = [max(max( shiftMat(:,:,2) )), max(max( shiftMat(:,:,3) ))]

figure('WindowState','max');
subplot(1,2,1);
imagesc( errorMat ); colorbar; axis square;
title( sprintf('Error (minimum total error = row %i)', mMin) )
subplot(1,2,2);
imagesc( shiftMat(:,:,1) ); colorbar; axis square;
set(gca, 'Xtick',1:Nmovie, 'Ytick',1:Nmovie);
title('Shift Distance (pixels)')
impixelinfo

% Determine largest pixel shifts in each dimension from the individually-registered movies
pooledShift = cat(1, shift{:});
pooledShift = pooledShift(:,1:2,regColor);
maxIndvShift = (ceil( max( abs(pooledShift) ) )); %ceil( prctile( pooledShift, 99, 1 ) ); % 
if any(maxIndvShift > 10), warning('Max individual shift exceeds 10 pixels in at least one dimension'); end


% Combine movies and  initial edge-trimming
fprintf('\nConcatenating movies... '); tic
pooledMovie = cat(4, regMovie{:}); % Concatenate the individually-registered movies
pooledMovie = pooledMovie(maxIndvShift(2)+1:end-maxIndvShift(2), maxIndvShift(1)+1:end-maxIndvShift(1), :, :); % Trim edges   metadata(1).goodColor
refIm = metadata(mMin).mean( maxIndvShift(2)+1:end-maxIndvShift(2), maxIndvShift(1)+1:end-maxIndvShift(1), : );
toc

% Register the full data and trim the edges
[regPooledMovie, regPooledShift, ~, ~, pooledRegParam] = RegisterMovie(pooledMovie, metadata(1), 'color',regColor, 'refIm', refIm ); % , 'save',saveRoot
%poolShiftZ = zscore( regPooledShift(:,1:3,regColor) );
%find( abs(poolShiftZ(:,2)) > 10 )

% Identify and suppress CSD-related outlier shifts
for m = 1:Nmovie % find(~cellfun( @isempty, CSDframes ))
    corrShift = regPooledShift( movieFrames{m}, :, : );
    CSDframes{m} = find( sum( abs(zscore( corrShift(:,1:2,regColor) )) > 10, 2) );
    if ~isempty( CSDframes{m} ) 
        fprintf('\nm = %i: Found %i CSD-related frames', m, numel(CSDframes{m}) ); 
        goodLocalFrames = 1:size(corrShift,1);  goodLocalFrames(CSDframes{m}) = [];
        for c = regColor
            corrShift(CSDframes{m},1,c) = interp1( goodLocalFrames, corrShift(goodLocalFrames,1,c), CSDframes{m} );
            corrShift(CSDframes{m},2,c) = interp1( goodLocalFrames, corrShift(goodLocalFrames,2,c), CSDframes{m} );
            corrShift(CSDframes{m},3,c) = sqrt( corrShift(CSDframes{m},1,c).^2 + corrShift(CSDframes{m},2,c).^2 );
            corrShift(CSDframes{m},4,c) = NaN;
        end
        CSDpoolFrames = CSDframes{m}+firstFrame(m)-1;
        regPooledShift(CSDpoolFrames,:,:) = corrShift(CSDframes{m},:,:); 
        corrShiftMovie = ApplyDFTshift( pooledMovie(:,:,:,CSDpoolFrames), regPooledShift(CSDpoolFrames,[1,2,5],:));
        regPooledMovie(:,:,:,CSDpoolFrames) = corrShiftMovie;
    end
end

maxPooledShift = ceil( max( abs(regPooledShift(:,1:2,pooledRegParam(1).color)) ) ); % ceil( prctile( regPooledShift(:,1:2,pooledRegParam(1).color), 99, 1 ) ); %
if any(maxPooledShift > 10), warning('Max pooled shift exceeds 10 pixels in at least one dimension'); end
trimPooledMovie = regPooledMovie(maxPooledShift(2)+1:end-maxPooledShift(2), maxPooledShift(1)+1:end-maxPooledShift(1), :, :); % Trim edges 

figure('WindowState','max');
subplot(2,2,1);
plot( regPooledShift(:,1,regColor) ); xlim([1,Inf]);
set(gca,'Xtick', firstFrame );
ylabel('X shift');
subplot(2,2,3);
plot( regPooledShift(:,2,regColor) ); xlim([1,Inf]); 
ylabel('Y shift');
subplot(2,2,[2,4]);
plot( regPooledShift(:,1,regColor), regPooledShift(:,2,regColor) ); 
xlabel('X shift'); ylabel('Y shift');
title( sprintf('All Frames Registered to average of movie %i (color channel %i)', mMin, regColor ) );


% Break the supermovie back into its individual consituents and save those movies to their individual folders
fprintf('\nWriting individual movies... ');

for m = 1:Nmovie
    tempMovie = trimPooledMovie(:,:,:,movieFrames{m});
    for c = metadata(m).goodColor
        savePath = sprintf('%s_pooledReg_chn%i.tif',metadata(m).root, c);
        saveastiff( squeeze(tempMovie(:,:,c,:)), savePath, GrayOpt ); 
        fprintf('\n   Wrote %s.  ', savePath ); toc
    end
end

% Save the results
if ~isempty( fovDir )
    % Save the registered supermovie (and other variables) to a mat file
    mkdir(fovDir)
    saveRoot = sprintf('%s%s', fovDir, fovName );
    fix(clock)%datetime('now')
    fprintf('\nSaving pooled/registered data to mat, and writing binned movies to tiff...  ');
    saveMovie = trimPooledMovie(:,:,metadata(1).goodColor,:); % remove unused colors to reduce file size
    save( sprintf('%s_pooled.mat', saveRoot ), 'saveMovie', 'pooledRegParam', 'regPooledShift', 'metadata', 'T', 'errorMat', 'shiftMat', 'mMin', 'Nframe', 'firstFrame', ...
        'maxIndvShift', 'refIm', 'fovDir', 'fovName', 'movieDirs', 'movieName', 'movieFrames', 'maxPooledShift', 'otherData', '-v7.3','-nocompression'); 
    toc
    
    % Write the registered super-movie to tif   
    for c = metadata(1).goodColor
        savePath = sprintf('%s_pooled_chn%i.tif',saveRoot, c);
        saveastiff( squeeze(trimPooledMovie(:,:,c,:)), savePath, GrayOpt ); 
        fprintf('Wrote %s.  ', savePath ); toc
    end
end

end

% Shift all frames of each individual movie to align with the mMin-th movie
%{
shiftedRegMovie = cell(1,Nmovie);
for n = 1:Nmovie
    tempShift = repmat( [shiftMat(mMin,n,2), shiftMat(mMin,n,3), shiftMat(mMin,n,4)], size(regMovie{m},4), 1, 3 ); 
    shiftedRegMovie{n} = ApplyDFTshift( regMovie{m}, tempShift );
end

fprintf('\nConcatenating movies... '); tic
catShiftMovie = cat(4, shiftedRegMovie{:}); % Concatenate the individually-registered movies

fprintf('\nWriting supermovie...\n');
GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
saveRoot = sprintf('%s%s', fovDir, fovName );
for c = metadata(1).goodColor
    savePath = sprintf('%s_shifted_chn%i.tif',saveRoot, c);
    saveastiff( squeeze(catShiftMovie(:,:,c,:)), savePath, GrayOpt ); 
    fprintf('Wrote %s.  ', savePath ); toc
end
%}