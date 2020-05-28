% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Calculate the Kolmogorov Smirnov distance between two Kaplan 
%              Meier estimates of the cumulative distribution, when provide 
%              a list of events, where the first and last events are censor
%              events. Left and right censoring both place a lower bound on
%              the time to an event and are treated equally. When passing 
%              different sequences of events separate the sequences by not %
%              a number (nan). Returns the maximum distance that the cumulative 
%              distribution of the first sample is above the cumulative 
%              distribution of the second sample, the maximum distance
%              that the cumulative distribution of the first sample is
%              below the cumulative distribution of the second sample, the
%              fraction of the total interval on which the first cumulative
%              distribution is completely above the second cumulative
%              distribution, and the fraction of the total interval where
%              the first cumulative distribution is completely below the second
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
% parameter..: array(float) EventsOne         Sequential lists of events in the first sample
% parameter..: array(float) EventsTwo         Sequential lists of events in the second sample
% return.....: array(float) DistanceStatistic Kolmogorov-Smirnov vector statistic
%
function DistanceStatistic = SequentialKaplanMeierStatistic(EventsOne, EventsTwo)
  
  % Escape with bad data
  if ~ (isnumeric(EventsOne) && isnumeric(EventsTwo))
    sprintf('%s', 'Data must be numeric.')
    DistanceStatistic = [];
    return
  end
  
% % Coerce the first event data  % % % % % % % % % % % % % % % % % % % % %
  EventsOne = reshape                     ...
  (                                       ...
    real (EventsOne(~ isinf(EventsOne))), ...
    numel(EventsOne(~ isinf(EventsOne))), ...
    1                                     ...
  );
  
  % Brace events with not a number boundaries
  EventsOne =  ...
  [            ...
    nan;       ...
    EventsOne; ...
    nan        ...
  ];

  % Remove all singleton events, no interval data is available
  SingletonIndex =                                                                               ...
  [                                                                                              ...
    false;                                                                                       ...
    isnan(EventsOne(1:(end - 2))) & (~ isnan(EventsOne(2:(end - 1)))) & isnan(EventsOne(3:end)); ...
    false                                                                                        ...
  ];
  EventsOne(SingletonIndex) = [];

  % Remove duplicate boundaries
  DuplicateIndex =                                          ...
  [                                                         ...
    false;                                                  ...
    isnan(EventsOne(1:(end - 1))) & isnan(EventsOne(2:end)) ...
  ];
  EventsOne(DuplicateIndex) = [];
  
  % Get the indicies of the boundaries
  BoundaryIndex = find(isnan(EventsOne));
  
  % Number of blocks of events and number of intervals between events
  BlockCount = numel(BoundaryIndex) - 1;
  CountOne   = numel(EventsOne(~ isnan(EventsOne))) - BlockCount;
  
  % Upper and lower boundary indices
  LowerInput  = BoundaryIndex(1:(end - 1), 1) + 1;
  UpperInput  = BoundaryIndex(2: end     , 1) - 1;
  LowerOutput = LowerInput - (1:2:(2 .* BlockCount - 1))';
  UpperOutput = UpperInput - (2:2:(2 .* BlockCount))';
  
  % Preallocate the sample
  SampleOne = [zeros(CountOne, 1), ones(CountOne, 1)];

  % Sort each block then find the time between events
  for BlockIndex = 1:BlockCount
    SampleOne(LowerOutput(BlockIndex, 1):UpperOutput(BlockIndex, 1), 1) = diff ...
    (                                                                          ...
      sort                                                                     ...
      (                                                                        ...
        EventsOne(LowerInput(BlockIndex, 1):UpperInput(BlockIndex, 1), 1)      ...
      )                                                                        ...
    );
    
    % First and last intervals are censor times
    SampleOne(LowerOutput(BlockIndex, 1), 2) = 0;
    SampleOne(UpperOutput(BlockIndex, 1), 2) = 0;
  end
  
  % Remove zero and negative intervals
  SampleOne(SampleOne(:, 1) <= 0, :) = [];
  
  % Recount positive intervals
  CountOne = size(SampleOne, 1); 
  
  % Sort the first events by event time
  SampleOne = sortrows(SampleOne, 1);
  
  % First sample Kaplan Meier estimate of the probability density of the
  SampleOne(:, 2) = SampleOne(:, 2) ./ (CountOne:-1:1)';
  SampleOne(:, 2) = SampleOne(:, 2) .*          ...
  [                                             ...
    1;                                          ...
    cumprod(1 - SampleOne(1:(CountOne - 1), 2)) ...
  ];
  
% % Coerce the second event data % % % % % % % % % % % % % % % % % % % % %
  EventsTwo = reshape                     ...
  (                                       ...
    real (EventsTwo(~ isinf(EventsTwo))), ...
    numel(EventsTwo(~ isinf(EventsTwo))), ...
    1                                     ...
  );
  
  % Brace events with not a number boundaries
  EventsTwo =  ...
  [            ...
    nan;       ...
    EventsTwo; ...
    nan        ...
  ];

  % Remove all singleton events, no interval data is available
  SingletonIndex =                                                                               ...
  [                                                                                              ...
    false;                                                                                       ...
    isnan(EventsTwo(1:(end - 2))) & (~ isnan(EventsTwo(2:(end - 1)))) & isnan(EventsTwo(3:end)); ...
    false                                                                                        ...
  ];
  EventsTwo(SingletonIndex) = [];

  % Remove duplicate boundaries
  DuplicateIndex =                                          ...
  [                                                         ...
    false;                                                  ...
    isnan(EventsTwo(1:(end - 1))) & isnan(EventsTwo(2:end)) ...
  ];
  EventsTwo(DuplicateIndex) = [];
  
  % Get the indicies of the boundaries
  BoundaryIndex = find(isnan(EventsTwo));
  
  % Number of blocks of events and number of intervals between events
  BlockCount = numel(BoundaryIndex) - 1;
  CountTwo   = numel(EventsTwo(~ isnan(EventsTwo))) - BlockCount;
  
  % Upper and lower boundary indices
  LowerInput  = BoundaryIndex(1:(end - 1), 1) + 1;
  UpperInput  = BoundaryIndex(2: end     , 1) - 1;
  LowerOutput = LowerInput - (1:2:(2 .* BlockCount - 1))';
  UpperOutput = UpperInput - (2:2:(2 .* BlockCount))';
  
  % Preallocate the sample
  SampleTwo = [zeros(CountTwo, 1), ones(CountTwo, 1)];

  % Sort each block then find the time between events
  for BlockIndex = 1:BlockCount
    SampleTwo(LowerOutput(BlockIndex, 1):UpperOutput(BlockIndex, 1), 1) = diff ...
    (                                                                          ...
      sort                                                                     ...
      (                                                                        ...
        EventsTwo(LowerInput(BlockIndex, 1):UpperInput(BlockIndex, 1), 1)      ...
      )                                                                        ...
    );
    
    % First and last intervals are censor times
    SampleTwo(LowerOutput(BlockIndex, 1), 2) = 0;
    SampleTwo(UpperOutput(BlockIndex, 1), 2) = 0;
  end
  
  % Remove zero and negative intervals
  SampleTwo(SampleTwo(:, 1) <= 0, :) = [];
  
  % Recount positive intervals
  CountTwo = size(SampleTwo, 1);
  
  % Sort the second events by event time
  SampleTwo = sortrows(SampleTwo, 1);
  
  % Second sample Kaplan Meier estimate of the probability density of the
  SampleTwo(:, 2) =   SampleTwo(:, 2) ./ (CountTwo:-1:1)';
  SampleTwo(:, 2) = - SampleTwo(:, 2) .*        ...
  [                                             ...
    1;                                          ...
    cumprod(1 - SampleTwo(1:(CountTwo - 1), 2)) ...
  ];  
  
% % Concatenate the matrices and sort by event time  % % % % % % % % % % %
  SortSample = sortrows ...
  (                     ...
    [                   ...
      SampleOne;        ... 
      SampleTwo         ...
    ],                  ...
    1                   ...
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
