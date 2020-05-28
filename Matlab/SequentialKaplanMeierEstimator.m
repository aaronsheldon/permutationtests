% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Calculate the Kaplan Meier estimator of the cumulative
%              distribution, when provide a list of events, where the first
%              and last events are censor events. Left and right censoring
%              both place a lower bound on the time to an event and are
%              treated equally. When passing Different sequences of events
%              separate the sequences by not a number (nan).
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
% parameter..: array(float) AllEvents             Sequential lists of events
% return.....: array(float) DistributionEstimator Time of events and cumulative frequency
%
function DistributionEstimator = SequentialKaplanMeierEstimator(AllEvents)
  
  % Escape with bad data
  if ~ isnumeric(AllEvents)
    sprintf('%s', 'Data must be numeric.')
    DistributionEstimator = [];
    return
  end
  
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
  
  % Number of blocks of events and number of intervals between events
  BlockCount = numel(BoundaryIndex) - 1;
  EventCount = numel(AllEvents(~ isnan(AllEvents))) - BlockCount;
  
  % Upper and lower boundary indices
  LowerInput  = BoundaryIndex(1:(end - 1), 1) + 1;
  UpperInput  = BoundaryIndex(2: end     , 1) - 1;
  LowerOutput = LowerInput - (1:2:(2 .* BlockCount - 1))';
  UpperOutput = UpperInput - (2:2:(2 .* BlockCount))';
  
  % Preallocate the sample
  SampleEvents = [zeros(EventCount, 1), ones(EventCount, 1)];

  % Sort each block then find the time between events
  for BlockIndex = 1:BlockCount
    SampleEvents (LowerOutput(BlockIndex, 1):UpperOutput(BlockIndex, 1), 1) = diff ...
    (                                                                              ...
      sort                                                                         ...
      (                                                                            ...
        AllEvents(LowerInput (BlockIndex, 1):UpperInput (BlockIndex, 1), 1)        ...
      )                                                                            ...
    );
    
    % First and last intervals are censor times
    SampleEvents(LowerOutput(BlockIndex, 1), 2) = 0;
    SampleEvents(UpperOutput(BlockIndex, 1), 2) = 0;
  end
  
  % Remove zero and negative intervals
  SampleEvents(SampleEvents(:, 1) <= 0, :) = [];
  
  % Recount positive intervals
  EventCount = size(SampleEvents, 1);
  
  % Sort the intervals
  SampleEvents = sortrows(SampleEvents, 1);
  
  % Generate the Kaplan Meier estimate of the cumulative probability function
  SampleEvents(:, 2) = SampleEvents(:, 2) ./ (EventCount:-1:1)';
  SampleEvents(:, 2) = cumsum                          ...
  (                                                    ...
    SampleEvents(:, 2) .*                              ...
    [                                                  ...
      1;                                               ...
      cumprod(1 - SampleEvents(1:(EventCount - 1), 2)) ...
    ]                                                  ...
  );
  
  % Return the results
  DistributionEstimator = SampleEvents;
end
