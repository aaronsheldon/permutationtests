% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Markov Chain Monte Carlo draws from the distribution of the 
%              Kolmogorov Distance Statistic, parameterized by sample size.
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
% parameter..: scalar(int) SampleCount       Number of samples in the cumulative probability distribution
% parameter..: scalar(int) DrawCount         Number of random draws
% return.....: array(int)  DistanceStatistic Random samples of the uniform norm distance
%
function DistanceStatistic = MCMCKolmogorov(SampleCount, DrawCount)

  % Escape with bad data
  if ~ (isnumeric(SampleCount) && isnumeric(DrawCount))
    sprintf('%s', 'Data must be numeric.')
    DistanceStatistic = [];
    return
  end

  % Rectify values
  DrawCount   = sum(ceil(abs(DrawCount  (~ (isnan(DrawCount)   | isinf(DrawCount))))));
  SampleCount = sum(ceil(abs(SampleCount(~ (isnan(SampleCount) | isinf(SampleCount))))));

  % Escape with bad data
  if SampleCount < 1 || DrawCount < 1
    sprintf('%s', 'Data out of bounds.')
    DistanceStatistic = [];
    return
  end
      
  % Debug: start the stopwatch
  tic;

  % The minimum size of the statistic
  MinValue = 1 ./ (2 .* SampleCount);

  % The mid points of the empirical distribution
  MidPoints = MinValue .* (1:2:(2 .* SampleCount - 1))';

  % Preallocate loop variable
  DistanceStatistic = zeros(DrawCount, 1);

  % Seed a random sort list
  SampleRun =                   ...
  [                             ...
    0;                          ...
    sort(rand(SampleCount, 1)); ...
    1                           ...
  ];

  % Find the largest distance
  DistanceStatistic(1, 1) = max(MinValue + abs(SampleRun(2:(SampleCount + 1), 1) - MidPoints));

  % Debug: display loop count and stopwatch
  sprintf('count(time) = %16.0f(%1.14f)', 1, toc)  
  
  % Run the simulations
  for DrawIndex = 2:DrawCount,

    % Debug: start the stopwatch
    tic;  

    % Pick two cards from the deck
    IndexOne = randi(SampleCount + 1, 1);    
    IndexTwo = randi(SampleCount    , 1);
    
    % First index is smaller
    if IndexOne < IndexTwo
      IndexLower = IndexOne;
      IndexUpper = IndexTwo;
    
    % Second index is smaller
    elseif IndexOne > IndexTwo
      IndexLower = IndexTwo;
      IndexUpper = IndexOne;
      
    % The indicies are equal
    else
      IndexLower = IndexOne;
      IndexUpper = SampleCount + 1;
    end
    
    % Build the shift of the intervals
    ShiftWeight   = rand(1, 1);
    IntervalLower = SampleRun(IndexLower + 1, 1) - SampleRun(IndexLower, 1);
    IntervalUpper = SampleRun(IndexUpper + 1, 1) - SampleRun(IndexUpper, 1);
    DeltaInterval = IntervalLower - ShiftWeight .* (IntervalLower + IntervalUpper);
    
    % Build the next sample
    SampleRun((1 + IndexLower):IndexUpper, 1) = SampleRun((1 + IndexLower):IndexUpper, 1) - DeltaInterval;

    % Compute the statistic
    DistanceStatistic(DrawIndex, 1) = max(MinValue + abs(SampleRun(2:(SampleCount + 1), 1) - MidPoints));
    
    % Debug: display loop count and stopwatch
    sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%16.0f(%1.14f)', DrawIndex, toc)    
  end 
end
