% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Markov Chain Monte Carlo draws from the multivariable hypergeometric
%              distibution, to uniformly sample from the possible outcomes
%              with fixed number of votes per party, and enumerated voters
%              per division. The first column is assumed to be the governing
%              party, the last column is assumed to be the number of enumerated
%              voters that did not cast a ballot. Returns Spearman's 
%              Correlation between the governing party and each of the other
%              parties, and the total enumerated voters.
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
% parameter..: array(int)  PartyTotals       List of total votes for each party
% parameter..: array(int)  TiesTwo           List of total votes in each electoral division
% parameter..: scalar(int) DrawCount         Number of random draws
% return.....: array(int)  ElectionStatistic Results of random elections
%
function ElectionStatistic = MCMCBallot(PartyTotals, DivisionTotals, DrawCount)
  
  % Escape with bad data
  if ~ (isnumeric(PartyTotals) && isnumeric(DivisionTotals) && isnumeric(DrawCount))
    sprintf('%s', 'Data must be numeric.')
    ElectionStatistic = [];
    return
  end

  % Rectify the data  
  DrawCount = sum(ceil(abs(DrawCount(~ (isnan(DrawCount) | isinf(DrawCount))))));
  
  % Coerce the party data
  PartyTotals = reshape                                                                       ...
  (                                                                                           ...
    ceil(abs((PartyTotals(~ (isnan(PartyTotals) | isinf(PartyTotals) | PartyTotals == 0))))), ...
    numel    (PartyTotals(~ (isnan(PartyTotals) | isinf(PartyTotals) | PartyTotals == 0))),   ...
    1                                                                                         ...
  );

  % Coerce the division data
  DivisionTotals = reshape                                                                              ...
  (                                                                                                     ...
    ceil(abs(DivisionTotals(~ (isnan(DivisionTotals) | isinf(DivisionTotals) | DivisionTotals == 0)))), ...
    numel   (DivisionTotals(~ (isnan(DivisionTotals) | isinf(DivisionTotals) | DivisionTotals == 0))),  ...
    1                                                                                                   ...
  ); 
  
  % Voter totals
  VoterPartyCount    = sum(PartyTotals);
  VoterDivisionCount = sum(DivisionTotals);
  
  % Escape with bad data
  if DrawCount < 1 || VoterPartyCount ~= VoterDivisionCount
    sprintf('%s', 'Data out of bounds.')
    ElectionStatistic = [];
    return
  end
    
  % Debug: start the stopwatch
  tic; 

  % Build the index of changes in the party allocation of voters
  PartyCount = numel(PartyTotals);
  PartyIndex = cumsum                 ...
  (                                   ...
    [                                 ...
      1;                              ...
      PartyTotals(1:(PartyCount - 1)) ...
    ]                                 ...
  );
  
  % Allocate voters to parties
  PartyVoters             = zeros(VoterPartyCount, 1);
  PartyVoters(PartyIndex) = 1;
  PartyVoters             = cumsum(PartyVoters);  

  % Build the index of changes in the division allocation of voters
  DivisionCount = numel(DivisionTotals); 
  DivisionIndex = cumsum                    ...
  (                                         ...
    [                                       ...
      1;                                    ...
      DivisionTotals(1:(DivisionCount - 1)) ...
    ]                                       ...
  );
  
  % Allocate voters to divisions
  DivsionVoters                = zeros(VoterDivisionCount, 1);
  DivsionVoters(DivisionIndex) = 1;
  DivisionVoters               = cumsum(DivsionVoters);
  
  % Seed with a shuffled deck
  DivisionVoters = DivisionVoters(randperm(VoterDivisionCount), 1);

  % Construct the outcome, rows are parties, columns are divisions
  ElectionOutcome = accumarray     ...
  (                                ...
    [PartyVoters, DivisionVoters], ...
    1,                             ...
    [PartyCount , DivisionCount ]  ...
  );

  % Initialize the loop
  ElectionStatistic = zeros(DrawCount, PartyCount, DivisionCount);

  % Compute the statistic
  ElectionStatistic(1, :, :) = ElectionOutcome;
  
  % Debug: display loop count and stopwatch
  sprintf('count(time) = %16.0f(%1.14f)', 1, toc)  
  
  % Loop through the simulations
  for DrawIndex = 2:DrawCount

    % Debug: start the stopwatch
    tic; 

    % Draw two cards from the divsions deck
    DivisionIndexOne = randi(DivisionCount    , 1);
    DivisionIndexTwo = randi(DivisionCount - 1, 1);
    
    % Make sure the cards were drawn without replacement
    if DivisionIndexOne == DivisionIndexTwo
      DivisionIndexTwo = DivisionCount;
    end
    
    % Draw two cards from the parties deck
    PartyIndexOne = randi(PartyCount    , 1);
    PartyIndexTwo = randi(PartyCount - 1, 1);
    
    % Make sure the cards were drawn without replacement
    if PartyIndexOne == PartyIndexTwo
      PartyIndexTwo = PartyCount;
    end

    % Compute the marginal totals
    DivisionTotalOne = ElectionOutcome(PartyIndexOne, DivisionIndexOne) + ElectionOutcome(PartyIndexTwo, DivisionIndexOne);
    DivisionTotalTwo = ElectionOutcome(PartyIndexOne, DivisionIndexTwo) + ElectionOutcome(PartyIndexTwo, DivisionIndexTwo);
    PartyTotalOne    = ElectionOutcome(PartyIndexOne, DivisionIndexOne) + ElectionOutcome(PartyIndexOne, DivisionIndexTwo);
    PartyTotalTwo    = ElectionOutcome(PartyIndexTwo, DivisionIndexOne) + ElectionOutcome(PartyIndexTwo, DivisionIndexTwo);
    
    % Generate a random hypergeometric draw for the marginals
    ElectionOutcome(PartyIndexOne, DivisionIndexOne) = hygernd(PartyTotalOne + PartyTotalTwo, PartyTotalOne, DivisionTotalOne);
    ElectionOutcome(PartyIndexOne, DivisionIndexTwo) = PartyTotalOne    - ElectionOutcome(PartyIndexOne, DivisionIndexOne);
    ElectionOutcome(PartyIndexTwo, DivisionIndexOne) = DivisionTotalOne - ElectionOutcome(PartyIndexOne, DivisionIndexOne);
    ElectionOutcome(PartyIndexTwo, DivisionIndexTwo) = DivisionTotalTwo - ElectionOutcome(PartyIndexOne, DivisionIndexTwo);
    
    % Compute the statistic
    ElectionStatistic(DrawIndex, :, :) = ElectionOutcome;
    
    % Debug: display loop count and stopwatch
    sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%16.0f(%1.14f)', DrawIndex, toc)
  end  
end
