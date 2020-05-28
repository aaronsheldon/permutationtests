% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Markov Chain Monte Carlo draws from a multinomial distribution
%              to find the expected maximum, minimum, and median enumerated 
%              voters in electoral divisions.
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
% parameter..: scalar(int) VoterCount     Number of voters to distribute among the electoral divisions
% parameter..: scalar(int) DivisionCount  Number of electoral divisions
% parameter..: scalar(int) DrawCount      Number of random draws
% return.....: array(int)  RangeStatistic Random samples of the minimum, median, and maximum of the electoral division enumeration
%
function RangeStatistic = MCMCEnumeration(VoterCount, DivisionCount, DrawCount)

  % Escape with bad data
  if ~ (isnumeric(VoterCount) && isnumeric(DivisionCount) && isnumeric(DrawCount))
    sprintf('%s', 'Data must be numeric.')
    RangeStatistic = [];
    return
  end

  % Rectify variables
  DrawCount     = sum(ceil(abs(DrawCount    (~ (isnan(DrawCount)     | isinf(DrawCount))))));
  VoterCount    = sum(ceil(abs(VoterCount   (~ (isnan(VoterCount)    | isinf(VoterCount))))));
  DivisionCount = sum(ceil(abs(DivisionCount(~ (isnan(DivisionCount) | isinf(DivisionCount))))));          

  % Escape with bad data
  if VoterCount < 1 || DivisionCount < 1 || DrawCount < 1
    sprintf('%s', 'Data out of bounds.')
    RangeStatistic = [];
    return
  end  
  
  % Debug: start the stopwatch
  tic;
  
  % Seed with a shuffled deck
  EnumerationCount = mnrnd(VoterCount, ones(1, DivisionCount) ./ DivisionCount);

  % Initialize the statistic
  RangeStatistic = zeros(DrawCount, 3);

  % Compute the statistic
  RangeStatistic(1, 1) = min   (EnumerationCount);
  RangeStatistic(1, 2) = median(EnumerationCount);
  RangeStatistic(1, 3) = max   (EnumerationCount);

  % Debug: display loop count and stopwatch
  sprintf('count(time) = %16.0f(%1.14f)', 1, toc)  
  
  % Loop through the simulations
  for DrawIndex = 2:DrawCount

    % Debug: start the stopwatch
    tic; 

    % Draw two cards from the deck
    IncrementIndex = randi(DivisionCount, 1);
    DecrementIndex = randi(DivisionCount, 1);

    % Binomial reassignment of enumerated voter counts
    if IncrementIndex ~= DecrementIndex
      TotalCount                          = EnumerationCount(1, DecrementIndex) + EnumerationCount(1, IncrementIndex);
      IncrementCount                      = binornd(TotalCount, 0.5);
      DecrementCount                      = TotalCount - IncrementCount;
      EnumerationCount(1, DecrementIndex) = DecrementCount;
      EnumerationCount(1, IncrementIndex) = IncrementCount;
    end

    % Compute the statistic
    RangeStatistic(DrawIndex, 1) = min   (EnumerationCount);
    RangeStatistic(DrawIndex, 2) = median(EnumerationCount);
    RangeStatistic(DrawIndex, 3) = max   (EnumerationCount);
    
    % Debug: display loop count and stopwatch
    sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%16.0f(%1.14f)', DrawIndex, toc)
  end  
end
