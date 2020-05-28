% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Monte Carlo draws from the distribution of the Kolmogorov 
%              Smirnov Distance Statistic, parameterized by sample sizes.
%              Returns the vector Kolmogorv-Smirnov Statistic between the
%              Kaplan-Meier Estimators: Distance first estimator is above
%              the second estimator; the distance the first estimator is
%              below the second estimator; the fraction of the curve where
%              the first estimator is above the second; and the fraction
%              of the curve where the first estimator is below the second.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Copyright 2010 Aaron Sheldon                                            %
%                                                                         %
% This program is free software: you can redistribute it and/or modify    % 
% it under the terms of the GNU General Public License as published by    %
% the Free Software Foundation version 3 of the License.                  %
%                                                                         %
% This program is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           %
% GNU General Public License for more details.                            %
%                                                                         %
% You should have received a copy of the GNU General Public License       %
% along with this program.  If not, see <http://www.gnu.org/licenses/>.   %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% parameter..: array(float) SamplePoints      Locations of all sample points
% parameter..: scalar(int)  SampleCount       Number of samples in the first cumulative probability distribution
% parameter..: scalar(int)  DrawCount         Number of random draws
% return.....: array(float) DistanceStatistic Random samples of the vector statistic
%
function DistanceStatistic = MonteCarloKolmogorovSmirnov(SamplePoints, SampleCount, DrawCount)

  % Escape with bad data
  if ~ (isnumeric(SamplePoints) && isnumeric(SampleCount) && isnumeric(DrawCount))
    sprintf('%s', 'Data must be numeric.')
    DistanceStatistic = [];
    return
  end

  % Rectify values
  DrawCount   = sum(ceil(abs(DrawCount  (~ (isnan(DrawCount)   | isinf(DrawCount))))));
  SampleCount = sum(ceil(abs(SampleCount(~ (isnan(SampleCount) | isinf(SampleCount))))));
  
  % Coerce the outcome data
  SamplePoints = reshape                                                ...
  (                                                                     ...
    abs  (SamplePoints(~ (isnan(SamplePoints) | isinf(SamplePoints)))), ...
    numel(SamplePoints(~ (isnan(SamplePoints) | isinf(SamplePoints)))), ...
    1                                                                   ...
  );
  
  % Rectify values
  SampleOneCount = SampleCount;
  SampleTwoCount = numel(SamplePoints) - SampleOneCount;  
  SampleCount    = SampleOneCount + SampleTwoCount;

  % Escape with bad data
  if SampleOneCount < 1 || SampleTwoCount < 1 || DrawCount < 1
    sprintf('%s', 'Data out of bounds.')
    DistanceStatistic = [];
    return
  end
      
  % Debug: start the stopwatch
  tic;

  % Sort by event times
  SamplePoints = sortrows(SamplePoints);

  % Calculate time intervals and drop the time column
  TotalInterval = SamplePoints(end, 1) - SamplePoints(1, 1);
  EventInterval =                                          ...
  [                                                        ...
    0;                                                     ...
    SamplePoints(2:end, 1) - SamplePoints(1:(end - 1), 1); ...
    0;                                                     ...
  ];
  
  % Preallocate loop variable
  DistanceStatistic = zeros(DrawCount, 4);
  CompareSample     =                            ...
  [                                              ...
      ones(SampleOneCount, 1) ./ SampleOneCount; ...
    - ones(SampleTwoCount, 1) ./ SampleTwoCount  ...
  ];    
 
  % Debug: display loop count and stopwatch
  sprintf('count(time) = %16.0f(%1.14f)', 0, toc)  
  
  % Run the simulations
  for DrawIndex = 1:DrawCount,

    % Debug: start the stopwatch
    tic;  

    % Shuffle the deck
    CompareSample = CompareSample(randperm(SampleCount), 1);

    % Partial sums for interval lengths and difference between curves
    SignedDifference =      ...
    [                       ...
      0;                    ...
      cumsum(CompareSample) ...
    ];
  
    % Initialize return variables
    AboveMaximum  = max(  SignedDifference);
    BelowMaximum  = max(- SignedDifference);
    AboveFraction = sum((SignedDifference > 0) .* EventInterval) ./ TotalInterval;
    BelowFraction = sum((SignedDifference < 0) .* EventInterval) ./ TotalInterval;

    % Return the result
    DistanceStatistic(DrawIndex, :) = [AboveMaximum, BelowMaximum, AboveFraction, BelowFraction];
    
    % Debug: display loop count and stopwatch
    sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%16.0f(%1.14f)', DrawIndex, toc)     
  end 
end
