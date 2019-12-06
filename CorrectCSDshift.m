function [shift, CSDframes] = CorrectCSDshift(shift)

% shift = shift values obtained from registration of unregMovie, using dftregistration. organized as (Nframe x 5 x 3)array of [ xshift, yshift, shift distance, error, phase difference ]

% Perform the corrections
corrShift = shift(:,[1:2,5],:); % [x shift, y shift, phase diff]
for c = find( ~isnan( shift(1,1,:) ) )'
    CSDframes = find( sum( abs(zscore( shift(:,1:2,c) )) > 10, 2) );
    fprintf('\nFound %i CSD-related frames', numel(CSDframes) );
    goodFrame = 1:size(shift,1);  goodFrame( CSDframes ) = [];
    % Use interpolation to correct shifts
    corrShift(CSDframes, 1, c) = interp1( goodFrame, shift(goodFrame,1,c), CSDframes, 'spline' ); % x shift
    corrShift(CSDframes, 2, c) = interp1( goodFrame, shift(goodFrame,2,c), CSDframes, 'spline' ); % y shift
    %corrShift(CSDframes, 3, c) = interp1( goodFrame, shift(goodFrame,5,c), CSDframes, 'spline' ); % phase shift - better to leave this one uncorrected

    % Plot the interpolated results vs original
    figure('WindowState','max');
    sp(1) = subplot(3,1,1); 
    h(1) = plot( shift(:,1,c), 'LineWidth', 1.5 ); hold on;
    h(2) = plot( corrShift(:,1,c) ); 
    ylabel('X Shift (pix)');
    legend(h, {'Raw','Corrected'}, 'Location','SouthEast' );
    
    sp(2) = subplot(3,1,2); 
    plot( shift(:,2,c), 'LineWidth', 1.5 ); hold on;
    plot( corrShift(:,2,c) ); 
    ylabel('Y Shift (pix)'); 
    
    sp(3) = subplot(3,1,3); 
    plot( shift(:,5,c), 'LineWidth', 1.5 ); hold on;
    plot( corrShift(:,3,c) ); 
    ylabel('Phase Shift'); 
	xlabel('Frame');
    linkaxes(sp,'x'); 
    xlim([1,Inf]);
    
    % Update values of shift array 
    shift(CSDframes,1:2,c) = corrShift(CSDframes,1:2,c);
    shift(CSDframes,3,c) = sqrt( corrShift(CSDframes,1,c).^2 + corrShift(CSDframes,2,c).^2 );
    shift(CSDframes,4,c) = NaN; % after the correction, we don't know the new image error
end

end

