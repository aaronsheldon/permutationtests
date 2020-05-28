% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Markov Chain Monte Carlo draws of the Spearman rank correlation
%              coefficient, calculated by permuting a list of ranks generated 
%              from a list of tie counts.
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
% parameter..: array(int)  TiesOne              List of tied scores in the first sample
% parameter..: array(int)  TiesTwo              List of tied scores in the second sample
% parameter..: scalar(int) DrawCount            Number of random draws
% return.....: array(int)  CorrelationStatistic Random samples rank correlation
%
function CorrelationStatistic = MCMCSpearman(TiesOne, TiesTwo, DrawCount)
  
  % Escape with bad data
  if ~ (isnumeric(TiesOne) && isnumeric(TiesTwo) && isnumeric(DrawCount))
    sprintf('%s', 'Data must be numeric.')
    CorrelationStatistic = [];
    return
  end

  % Rectify the data
  DrawCount = sum(ceil(abs(DrawCount(~ (isnan(DrawCount) | isinf(DrawCount))))));
  
  % Coerce the first ties
  TiesOne = reshape                                                         ...
  (                                                                         ...
    ceil(abs(TiesOne(~ (isnan(TiesOne) | isinf(TiesOne) | TiesOne == 0)))), ...
    numel   (TiesOne(~ (isnan(TiesOne) | isinf(TiesOne) | TiesOne == 0))),  ...
    1                                                                       ...
  );

  % Coerce the second ties
  TiesTwo = reshape                                                         ...
  (                                                                         ...
    ceil(abs(TiesTwo(~ (isnan(TiesTwo) | isinf(TiesTwo) | TiesTwo == 0)))), ...
    numel   (TiesTwo(~ (isnan(TiesTwo) | isinf(TiesTwo) | TiesTwo == 0))),  ...
    1                                                                       ...
  );    
  
  % Total sample size
  SampleOneCount = sum(TiesOne);
  SampleTwoCount = sum(TiesTwo);
  
  % Escape with bad data
  if DrawCount < 1 || SampleOneCount ~= SampleTwoCount
    sprintf('%s', 'Data out of bounds.')
    CorrelationStatistic = [];
    return
  end
    
  % Debug: start the stopwatch
  tic;
    
  % Build the list of individual ranks from the list of ties, for the first variable
  RanksOne = zeros(SampleOneCount, 1);
  SizeOne  = numel (TiesOne);
  UpperOne = cumsum(TiesOne);
  LowerOne = 1 +              ...
  [                           ...
    0;                        ...
    UpperOne(1:(SizeOne - 1)) ...
  ];
  
  
  % Assign the mean square error rank of each set of ties
  for IndexOne = 1:SizeOne
    RanksOne(LowerOne(IndexOne):UpperOne(IndexOne), 1) = LowerOne(IndexOne) + (TiesOne(IndexOne) - 1) ./ 2;
  end
  RanksOne = (RanksOne - mean(RanksOne)) ./ sqrt(sum((RanksOne - mean(RanksOne)) .^ 2));

  % Build the list of individual ranks from the list of ties, for the second variable
  RanksTwo = zeros(SampleTwoCount, 1);
  SizeTwo  = numel (TiesTwo);    
  UpperTwo = cumsum(TiesTwo);
  LowerTwo = 1 +              ...
  [                           ...
    0;                        ...
    UpperTwo(1:(SizeTwo - 1)) ...
  ];
  
  % Assign the mean square error rank of each set of ties
  for IndexTwo = 1:SizeTwo
    RanksTwo(LowerTwo(IndexTwo):UpperTwo(IndexTwo), 1) = LowerTwo(IndexTwo) + (TiesTwo(IndexTwo) - 1) ./ 2;
  end
  RanksTwo = (RanksTwo - mean(RanksTwo)) ./ sqrt(sum((RanksTwo - mean(RanksTwo)) .^ 2));

  % Initialize the loop
  CorrelationStatistic = zeros(DrawCount, 1);

  % Seed with a shuffled deck
  RanksTwo = RanksTwo(randperm(SampleTwoCount), 1);

  % Compute the statistic
  CorrelationStatistic(1, 1) = RanksOne' * RanksTwo;

  % Debug: display loop count and stopwatch
  sprintf('count(time) = %16.0f(%1.14f)', 1, toc)  
  
  % Loop through the simulations
  for DrawIndex = 2:DrawCount

    % Debug: start the stopwatch
    tic; 

    % Swap cards in the deck
    SwapIndex              = randi(SampleTwoCount, 2, 1);
    RanksTwo(SwapIndex, 1) = RanksTwo(SwapIndex([2; 1], 1), 1);

    % Compute the statistic
    CorrelationStatistic(DrawIndex, 1) = RanksOne' * RanksTwo;
    
    % Debug: display loop count and stopwatch
    sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%16.0f(%1.14f)', DrawIndex, toc)   
  end  
end
