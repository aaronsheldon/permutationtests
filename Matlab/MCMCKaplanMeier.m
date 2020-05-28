% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Markov Chain Monte Carlo draws from the uniform distribution on the
%              combinations of finite sized samples, to simulate all possible 
%              allocations of survival data into to fixed sized samples.
%              Returns the vector Kolmogorv-Smirnov Statistic between the
%              Kaplan-Meier Estimators: Distance first estimator is above
%              the second estimator; the distance the first estimator is
%              below the second estimator; the fraction of the curve
%              where the first estimator is above the second; and the
%              fraction of the curve where the first estimator is below the
%              second.
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
% parameter..: array(float) OutcomePoints     Times of outcome events
% parameter..: array(float) CensorPoints      Times of censor events
% parameter..: scalar(int)  SampleCount       Size of the first sample events
% parameter..: scalar(int)  DrawCount         Number of random draws
% return.....: array(float) DistanceStatistic Random samples of the Kolmogorov-Smirnov vector statistic
%
function DistanceStatistic = MCMCKaplanMeier(OutcomePoints, CensorPoints, SampleCount, DrawCount)  

  % Escape with bad data
  if ~ (isnumeric(OutcomePoints) && isnumeric(CensorPoints) && isnumeric(SampleCount) && isnumeric(DrawCount))
    sprintf('%s', 'Data must be numeric.')
    DistanceStatistic = [];
    return
  end

  % Rectify values
  DrawCount   = sum(ceil(abs(DrawCount  (~ (isnan(DrawCount)   | isinf(DrawCount))))));
  SampleCount = sum(ceil(abs(SampleCount(~ (isnan(SampleCount) | isinf(SampleCount))))));
  
  % Coerce the outcome data
  OutcomePoints = reshape                                                  ...
  (                                                                        ...
    abs  (OutcomePoints(~ (isnan(OutcomePoints) | isinf(OutcomePoints)))), ...
    numel(OutcomePoints(~ (isnan(OutcomePoints) | isinf(OutcomePoints)))), ...
    1                                                                      ...
  );

  % Coerce the censor data
  CensorPoints = reshape                                                ...
  (                                                                     ...
    abs  (CensorPoints(~ (isnan(CensorPoints) | isinf(CensorPoints)))), ...
    numel(CensorPoints(~ (isnan(CensorPoints) | isinf(CensorPoints)))), ...
    1                                                                   ...
  );  
  
  % Total number of data points in each sample
  OutcomeCount   = numel(OutcomePoints);
  CensorCount    = numel(CensorPoints);
  EventCount     = OutcomeCount + CensorCount;
  SampleOneCount = SampleCount;
  SampleTwoCount = EventCount - SampleCount;
  
  % Escape with bad data
  if DrawCount < 1 || EventCount <= SampleCount
    sprintf('%s', 'Data out of bounds.')
    DistanceStatistic = [];
    return
  end
  
  % Debug: start the stopwatch
  tic;

  % Sort by event times 
  EventPoints = sortrows                      ...
  (                                           ...
    [                                         ...
      [OutcomePoints, ones(OutcomeCount, 1)]; ...
      [CensorPoints , zeros(CensorCount, 1)]  ...
    ],                                        ...
    1                                         ...
  );

  % Calculate time intervals and drop the time column
  TotalInterval = EventPoints(end, 1) - EventPoints(1, 1);
  EventInterval =                                         ...
  [                                                       ...
    0;                                                    ...
    EventPoints(2:end, 1) - EventPoints(1:(end - 1), 1); ...
    0;                                                    ...
  ];
  EventPoints(:, 1) = [];

  % Initialize the loop
  DistanceStatistic = zeros(DrawCount , 4);
  CompareSample     = zeros(EventCount, 1);   

  % Seed by shuffling and splitting deck, then sorting the splits
  IndexAll = randperm(EventCount)';
  IndexOne = sort(IndexAll(1:SampleOneCount               , 1));
  IndexTwo = sort(IndexAll((SampleOneCount + 1):EventCount, 1));

  % Probability density of the first hand dealt
  CompareSample(IndexOne, 1) = EventPoints  (IndexOne, 1) ./ (SampleOneCount:-1:1)';
  CompareSample(IndexOne, 1) = CompareSample(IndexOne, 1) .*           ...
  [                                                                    ...
    1;                                                                 ...
    cumprod(1 - CompareSample(IndexOne(1:(SampleOneCount - 1), 1), 1)) ...
  ];
   
  % Probability density of the second hand dealt
  CompareSample(IndexTwo, 1) =   EventPoints  (IndexTwo, 1) ./ (SampleTwoCount:-1:1)';
  CompareSample(IndexTwo, 1) = - CompareSample(IndexTwo, 1) .*         ...
  [                                                                    ...
    1;                                                                 ...
    cumprod(1 - CompareSample(IndexTwo(1:(SampleTwoCount - 1), 1), 1)) ...
  ];

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
  DistanceStatistic(1, :) = [AboveMaximum, BelowMaximum, AboveFraction, BelowFraction];

  % Debug: display loop count and stopwatch
  sprintf('count(time) = %16.0f(%1.14f)', 1, toc)  
  
  % Loop through the simulations
  for DrawIndex = 2:DrawCount

    % Debug: start the stopwatch
    tic;

    % Markov Chain Monte Carlo swapping of indexes in samples
    SwapOne = IndexOne(randi(SampleOneCount, 1), 1);
    SwapTwo = IndexTwo(randi(SampleTwoCount, 1), 1);
    
    % Assign the swapped index in the first hand
    IndexOne =                                                ...
    [                                                         ...
      IndexOne((IndexOne < SwapTwo) & (IndexOne ~= SwapOne)); ...
      SwapTwo;                                                ...
      IndexOne((IndexOne > SwapTwo) & (IndexOne ~= SwapOne))  ...
    ];
    
    % Assign the swapped index in the second hand
    IndexTwo =                                                ...
    [                                                         ...
      IndexTwo((IndexTwo < SwapOne) & (IndexTwo ~= SwapTwo)); ...
      SwapOne;                                                ...
      IndexTwo((IndexTwo > SwapOne) & (IndexTwo ~= SwapTwo))  ...
    ];

    % Add the first hand to the comparison array
    CompareSample(IndexOne, 1) = EventPoints  (IndexOne, 1) ./ (SampleOneCount:-1:1)';
    CompareSample(IndexOne, 1) = CompareSample(IndexOne, 1) .*           ...
    [                                                                    ...
      1;                                                                 ...
      cumprod(1 - CompareSample(IndexOne(1:(SampleOneCount - 1), 1), 1)) ...
    ];
    
    % Add the second hand to the comparison array
    CompareSample(IndexTwo, 1) =   EventPoints  (IndexTwo, 1) ./ (SampleTwoCount:-1:1)';
    CompareSample(IndexTwo, 1) = - CompareSample(IndexTwo, 1) .*         ...
    [                                                                    ...
      1;                                                                 ...
      cumprod(1 - CompareSample(IndexTwo(1:(SampleTwoCount - 1), 1), 1)) ...
    ];

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
