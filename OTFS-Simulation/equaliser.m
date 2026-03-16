function [eqSig] = equaliser(rxSig,fadedSig,txSig,ofdmSym)

%--------------------------------------------------------------------------
%
%                   Equalises a faded signal 
%
%--------------------------------------------------------------------------
% Input arguments: 
%
% rxSig                         The received faded signal w/ noise applied
% fadedSig                      The received faded signal, NO noise
% txSig                         Transmitted signal prior to channel
% ofdmSym                       Totol ofdm symbols per subframe
% 
%--------------------------------------------------------------------------
% Function returns: 
% 
% eqSig                         The equalised output signal (serial vector)
%
%--------------------------------------------------------------------------
%
% Author: Bradley Bates
% University of Bristol, UK
% email address: bb16177@bristol.ac.uk
% May 2020
%
% Copyright (c) 2020, Bradley Bates
%
%--------------------------------------------------------------------------

% Input parameters
pilotSpacing = 4;                       % Subcarrier sampling period (4 is best)
noSamples = numel(rxSig)/ofdmSym;       % No. of times the channel is sampled
numSC = pow2(floor(log2(noSamples)));   % Calc. no. of OFDM subcarrier

% Reformat input data into a matrix to make sampling easier
fadedSig_matrix = reshape(fadedSig, [noSamples,ofdmSym]);
txSig_matrix = reshape(txSig, [noSamples,ofdmSym]);

% Remove cyclic prefix as do no need to equalise it
fadedSig_matrix(1:(noSamples-numSC),:) = [];   % Remove cp
txSig_matrix(1:(noSamples-numSC),:) = [];   % Remove cp

% Sample signal
h = zeros(size(txSig_matrix));                          % Pre allocate matrix
for i = 1:size(txSig_matrix,1)                          % Loop rows
    for j = 1:size(txSig_matrix,2)                      % Loop Columns
        % Sample every pilotSpacing subcarrier (incuding first and last subcarriers)
        if i==1
             sample(i,j) = fadedSig_matrix(i,j)./txSig_matrix(i,j);
        elseif 0 == mod(i,pilotSpacing)
            sample(i/pilotSpacing+1,j) = fadedSig_matrix(i,j)./txSig_matrix(i,j);
        end
    end    
end

% Linearly Interpolate samples between subcarriers 
interpPoints = ( (1+(1/pilotSpacing)) : (1/pilotSpacing) : size(sample,1) );% Calc. points to interpolate
for j = 1:size(txSig_matrix,2)
    h(:,j) = interp1(sample(:,j),interpPoints);
end

% Keep only 1st, middle and last samples per carrier
midSym = ceil(ofdmSym/2);
h1 = [h(:,1)'; h(:,midSym)'; h(:,ofdmSym)']';

% Linearly Interpolate samples between symbols
seg1Len = midSym;          % number of symbols in first segment (1 to midSym)
seg2Len = ofdmSym - midSym; % number of symbols in second segment (midSym to end)
interpPoints1 = linspace(1, 2, seg1Len);
interpPoints2 = linspace(2, 3, seg2Len + 1);
interpPoints2(1) = [];     % avoid duplicating midSym
for i = 1:size(txSig_matrix,1)
    h_seg1(i,:) = interp1(h1(i,:), interpPoints1);
    h_seg2(i,:) = interp1(h1(i,:), interpPoints2);
end

% Concatenate segments
h = [h_seg1'; h_seg2']';

% Set eq of CP carriers to same as 1st carrier
for k = 1:(noSamples-numSC)
    h = [h(1,:);h];
end

% Convert back into serial
H = reshape(h,[numel(h),1]);
H(isnan(H)) = 1;                  % Set all NaN's to 1

% Equalise the signal 
eqSig = rxSig ./ H;
eqSig(isnan(eqSig))=0;            % Set all NaN's to 0

end