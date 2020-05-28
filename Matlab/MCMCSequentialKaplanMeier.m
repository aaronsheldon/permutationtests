% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: MArkov Chain Monte Carlo draws from the uniform distribution on the
%              combinations of finite sized samples, to simulate all possible 
%              allocations of survival data into to fixed sized samples.
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
% parameter..: array(float) AllEvents         Sequential list of even
% parameter..: scalar(int)  SampleCount       Size of the first sample of events
% parameter..: scalar(int)  DrawCount         Number of random draws
% return.....: array(float) DistanceStatistic Random samples of the Kolmogorov-Smirnov vector statistic
%
function DistanceStatistic = MCMCSequentialKaplanMeier(AllEvents, SampleCount, DrawCount)
  
  % Escape with bad data
  if ~ (isnumeric(AllEvents) && isnumeric(SampleCount) && isnumeric(DrawCount))
    sprintf('%s', 'Data must be numeric.')
    DistanceStatistic = [];
    return
  end

  % Rectify values
  DrawCount   = sum(ceil(abs(DrawCount  (~ (isnan(DrawCount)   | isinf(DrawCount))))));
  SampleCount = sum(ceil(abs(SampleCount(~ (isnan(SampleCount) | isinf(SampleCount))))));
  
  % Coerce the event data
  AllEvents = reshape                     ...
  (                                       ...
    real (AllEvents(~ isinf(AllEvents))), ...
    numel(AllEvents(~ isinf(AllEvents))), ...
    1                                     ...
  );
  
  % Brace events with not a number boundaries
  AllEvents =  ...
  [            ...
    nan;       ...
    AllEvents; ...
    nan        ...
  ];

  % Remove all singleton events, no interval data is available
  SingletonIndex =                                                                               ...
  [                                                                                              ...
    false;                                                                                       ...
    isnan(AllEvents(1:(end - 2))) & (~ isnan(AllEvents(2:(end - 1)))) & isnan(AllEvents(3:end)); ...
    false                                                                                        ...
  ];
  AllEvents(SingletonIndex) = [];

  % Remove duplicate boundaries
  DuplicateIndex =                                          ...
  [                                                         ...
    false;                                                  ...
    isnan(AllEvents(1:(end - 1))) & isnan(AllEvents(2:end)) ...
  ];
  AllEvents(DuplicateIndex) = [];
  
  % Get the indicies of the boundaries
  BoundaryIndex = find(isnan(AllEvents));  
  BlockCount    = numel(BoundaryIndex) - 1;
  BlockOneCount = SampleCount;
  BlockTwoCount = BlockCount - SampleCount;
  
  % Escape with bad data
  if DrawCount < 1 || BlockCount <= SampleCount
    sprintf('%s', 'Data out of bounds.')
    DistanceStatistic = [];
    return
  end
  
  % Debug: start the stopwatch
  tic;
  
  % Upper and lower boundary indices
  LowerInput  = BoundaryIndex(1:(end - 1), 1) + 1;
  UpperInput  = BoundaryIndex(2: end     , 1) - 1;
  OutputCount = UpperInput - LowerInput;
  
  % Preallocate the sample blocks
  EventBlocks = cell(BlockCount, 1);
  
  % Calculate the suvival and censor intervals per block 
  for BlockIndex = 1:BlockCount
    SampleEvents = [zeros(OutputCount(BlockIndex, 1), 1), ones(OutputCount(BlockIndex, 1), 1)];
    
    % Find the time between events and set first and last intervals as censor times
    SampleEvents(:  , 1) = diff(sort(AllEvents(LowerInput(BlockIndex, 1):UpperInput(BlockIndex, 1), 1)));
    SampleEvents(1  , 2) = 0;
    SampleEvents(end, 2) = 0;
    
    % Remove zero length events
    SampleEvents(SampleEvents(:, 1) <= 0, :) = [];
  
    % Sort the intervals
    SampleEvents = sortrows(SampleEvents, 1);
    
    % Store the block of events
    EventBlocks{BlockIndex}=SampleEvents;
  end
      
  % Sort by event times to calculate time intervals
  EventPoints   = sortrows(vertcat(EventBlocks{:, 1}), 1);
  TotalInterval = EventPoints(end, 1) - EventPoints(1, 1);
  EventInterval =                                        ...
  [                                                      ...
    0;                                                   ...
    EventPoints(2:end, 1) - EventPoints(1:(end - 1), 1); ...
    0;                                                   ...
  ];

  % Initialize the loop
  DistanceStatistic = zeros(DrawCount , 4); 

  % Shuffle and split the deck, then sort the splits
  IndexAll = randperm(BlockCount)';
  IndexOne = IndexAll(1:SampleCount               , 1);
  IndexTwo = IndexAll((SampleCount + 1):BlockCount, 1);

  % Get the first set of blocks of events
  SampleOne      = sortrows(vertcat(EventBlocks{IndexOne, :}), 1);
  SampleOneCount = size(SampleOne, 1);
    
  % Probability density of the first hand dealt
  SampleOne(:, 2) = SampleOne(:, 2) ./ (SampleOneCount:-1:1)';
  SampleOne(:, 2) = SampleOne(:, 2) .*                ...
  [                                                   ...
    1;                                                ...
    cumprod(1 - SampleOne(1:(SampleOneCount - 1), 2)) ...
  ];
    
  % Get the second set of blocks of events
  SampleTwo      = sortrows(vertcat(EventBlocks{IndexTwo, :}), 1);
  SampleTwoCount = size(SampleTwo, 1);
   
  % Probability density of the second hand dealt
  SampleTwo(:, 2) =   SampleTwo(:, 2) ./ (SampleTwoCount:-1:1)';
  SampleTwo(:, 2) = - SampleTwo(:, 2) .*                ...
  [                                                     ...
    1;                                                  ...
    cumprod(1 - SampleTwo(1:(SampleTwoCount - 1), 2))   ...
  ];
  
  % Concatenate the matrices and sort by event time
  CompareSample = sortrows ...
  (                        ...
    [                      ...
      SampleOne;           ... 
      SampleTwo            ...
    ],                     ...
    1                      ...
  ); 
  
  % Partial sums for interval lengths and difference between curves
  SignedDifference =            ...
  [                             ...
    0;                          ...
    cumsum(CompareSample(:, 2)) ...
  ];
  
  % Initialize return variables
  AboveMaximum  = max(  SignedDifference);
  BelowMaximum  = max(- SignedDifference);
  AboveFraction = sum((SignedDifference > 0) .* EventInterval) ./ TotalInterval;
  BelowFraction = sum((SignedDifference < 0) .* EventInterval) ./ TotalInterval;

  % Return the result
  DistanceStatistic(1, :) = [AboveMaximum, BelowMaximum, AboveFraction, BelowFraction];
  
  % Debug: display loop count and stopwatch
  sprintf('count(time) = %16.0f(%1.14f)', 0, toc)  
  
  % Loop through the simulations
  for DrawIndex = 1:DrawCount

    % Debug: start the stopwatch
    tic;

    % Markov Chain Monte Carlo swapping of indexes in samples
    SwapOne = IndexOne(randi(BlockOneCount, 1), 1);
    SwapTwo = IndexTwo(randi(BlockTwoCount, 1), 1);
    
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

    % Get the first set of blocks of events
    SampleOne      = sortrows(vertcat(EventBlocks{IndexOne, :}), 1);
    SampleOneCount = size(SampleOne, 1);
    
    % Probability density of the first hand dealt
    SampleOne(:, 2) = SampleOne(:, 2) ./ (SampleOneCount:-1:1)';
    SampleOne(:, 2) = SampleOne(:, 2) .*                ...
    [                                                   ...
      1;                                                ...
      cumprod(1 - SampleOne(1:(SampleOneCount - 1), 2)) ...
    ];
    
    % Get the second set of blocks of events
    SampleTwo      = sortrows(vertcat(EventBlocks{IndexTwo, :}), 1);
    SampleTwoCount = size(SampleTwo, 1);
    
    % Probability density of the second hand dealt
    SampleTwo(:, 2) =   SampleTwo(:, 2) ./ (SampleTwoCount:-1:1)';
    SampleTwo(:, 2) = - SampleTwo(:, 2) .*                ...
    [                                                     ...
      1;                                                  ...
      cumprod(1 - SampleTwo(1:(SampleTwoCount - 1), 2))   ...
    ];
  
    % Concatenate the matrices and sort by event time
    CompareSample = sortrows ...
    (                        ...
      [                      ...
        SampleOne;           ... 
        SampleTwo            ...
      ],                     ...
      1                      ...
    ); 
  
    % Partial sums for interval lengths and difference between curves
    SignedDifference =            ...
    [                             ...
      0;                          ...
      cumsum(CompareSample(:, 2)) ...
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
