% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Calculate the Kaplan Meier estimator of the cumulative distribution.
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
% parameter..: array(float) OutcomeEvents         Times of outcome events
% parameter..: array(float) CensorEvents          Times of censor events
% return.....: array(float) DistributionEstimator Time of events and cumulative frequency
%
function DistributionEstimator = KaplanMeierEstimator(OutcomeEvents, CensorEvents)
  
  % Escape with bad data
  if ~ (isnumeric(OutcomeEvents) && isnumeric(CensorEvents))
    sprintf('%s', 'Data must be numeric.')
    DistributionEstimator = [];
    return
  end
  
  % Coerce the outcome data
  OutcomeEvents = reshape                                                  ...
  (                                                                        ...
    abs  (OutcomeEvents(~ (isnan(OutcomeEvents) | isinf(OutcomeEvents)))), ...
    numel(OutcomeEvents(~ (isnan(OutcomeEvents) | isinf(OutcomeEvents)))), ...
    1                                                                      ...
  );

  % Coerce the censor data
  CensorEvents = reshape                                                ...
  (                                                                     ...
    abs  (CensorEvents(~ (isnan(CensorEvents) | isinf(CensorEvents)))), ...
    numel(CensorEvents(~ (isnan(CensorEvents) | isinf(CensorEvents)))), ...
    1                                                                   ...
  );

  % Total data points in the sample
  OutcomeCount = numel(OutcomeEvents);
  CensorCount  = numel(CensorEvents);
  CountEvents  = OutcomeCount + CensorCount;
  
  % Assign censor flags and sort by event time
  SampleEvents = sortrows                     ...
  (                                           ...
    [                                         ...
      [OutcomeEvents, ones(OutcomeCount, 1)]; ...
      [CensorEvents , zeros(CensorCount, 1)]  ...
    ],                                        ...
    1                                         ...
  );
  
  % Generate the Kaplan Meier estimate of the cumulative probability function
  SampleEvents(:, 2) = SampleEvents(:, 2) ./ (CountEvents:-1:1)';
  SampleEvents(:, 2) = cumsum                           ...
  (                                                     ...
    SampleEvents(:, 2) .*                               ...
    [                                                   ...
      1;                                                ...
      cumprod(1 - SampleEvents(1:(CountEvents - 1), 2)) ...
    ]                                                   ...
  );
  
  % Return the results
  DistributionEstimator = SampleEvents;
end
