% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Calculate the Kolmogorov Smirnov difference statistic of two 
%              samples. Returns the maximum distance that the cumulative 
%              distribution of the first sample is above the cumulative 
%              distribution of the second sample, the maximum distance
%              that the cumulative distribution of the first sample is
%              below the cumulative distribution of the second sample, the
%              fraction of the total interval on which the first cumulative
%              distribution is completely above the second cumulative
%              distribution, and the fraction of the total interval where
%              the first cumulativedistribution is completely below the second
%              cumulative distribution.
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
% parameter..: array(float) SampleOne         First sample of measurements
% parameter..: array(float) SampleTwo         Second sample of measurements
% return.....: array(float) DistanceStatistic Kolmogorov-Smirnov vector statistic
%
function DistanceStatistic = KolmogorovSmirnovStatistic(SampleOne, SampleTwo)
  
  % Escape with bad data
  if ~ (isnumeric(SampleOne) && isnumeric(SampleTwo))
    sprintf('%s', 'Data must be numeric.')
    DistanceStatistic = [];
    return
  end

  % Coerce the first data
  SampleOne = reshape                                          ...
  (                                                            ...
    real (SampleOne(~ (isnan(SampleOne) | isinf(SampleOne)))), ...
    numel(SampleOne(~ (isnan(SampleOne) | isinf(SampleOne)))), ...
    1                                                          ...
  );
  
  % Coerce the second data
  SampleTwo = reshape                                          ...
  (                                                            ...
    real (SampleTwo(~ (isnan(SampleTwo) | isinf(SampleTwo)))), ...
    numel(SampleTwo(~ (isnan(SampleTwo) | isinf(SampleTwo)))), ...
    1                                                          ...
  );

  % Get the total number of data points in each sample
  CountOne = numel(SampleOne);
  CountTwo = numel(SampleTwo);

  % Concatenate the matrices and sort
  SortSample = sortrows                             ...
  (                                                 ...
    [                                               ...
      [SampleOne,   ones(CountOne, 1) ./ CountOne]; ...
      [SampleTwo, - ones(CountTwo, 1) ./ CountTwo]  ...
    ],                                              ...
    1                                               ...
  );

  % Simple intervals for measuing dwell fraction
  TotalInterval = SortSample(end, 1) - SortSample(1, 1);
  EventInterval =                                      ...
  [                                                    ...
    0;                                                 ...
    SortSample(2:end, 1) - SortSample(1:(end - 1), 1); ...
    0;                                                 ...
  ];
  
  % Partial sums for interval lengths and difference between curves
  SignedDifference =         ...
  [                          ...
    0;                       ...
    cumsum(SortSample(:, 2)) ...
  ];
  
  % Initialize return variables
  AboveMaximum  = max(  SignedDifference);
  BelowMaximum  = max(- SignedDifference);
  AboveFraction = sum((SignedDifference > 0) .* EventInterval) ./ TotalInterval;
  BelowFraction = sum((SignedDifference < 0) .* EventInterval) ./ TotalInterval;

  % Return the result
  DistanceStatistic = [AboveMaximum, BelowMaximum, AboveFraction, BelowFraction];
end
