function outputMovie = ApplyDFTshift( inputMovie, Shift )
%inputMovie = movie to be shifted
%shift = (Nframe x 3 x 3)array of [ xshift, yshift, phase difference ] 
tic
Nframe = size(inputMovie,4);
[Nrow,Ncol] = size(inputMovie(:,:,1,1)); % adapted from dftregistration
Nr = ifftshift( -fix(Nrow/2):ceil(Nrow/2)-1 ); % adapted from dftregistration
Nc = ifftshift( -fix(Ncol/2):ceil(Ncol/2)-1 ); % adapted from dftregistration
[Nc,Nr] = meshgrid(Nc,Nr); % adapted from dftregistration

% Apply the shift to each frame of the movie
cShift = find( inputMovie(1,1,:,1) > 0 )';
outputMovie = inputMovie;
for z = flip(1:Nframe)
    for c = cShift
        depFT = fft2( inputMovie(:,:,c,z) );
        fftDepReg = depFT.*exp(1i*2*pi*(-Shift(z,2,c)*Nr/Nrow-Shift(z,1,c)*Nc/Ncol)); % adapted from dftregistration
        fftDepReg = fftDepReg*exp(1i*Shift(z,3,c)); % adapted from dftregistration
        outputMovie(:,:,c,z) = abs( ifft2(fftDepReg) );
    end
end
toc
end

