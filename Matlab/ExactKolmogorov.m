% Author.....: Aaron Sheldon
% Date.......: April, 2010
% Description: Implemention of the determinant algorythm of K. Sarkadi,
%              to find the extact cumulative distribution of the 
%              Kolmogorov statistic, for large samples uses recursive 
%              gaussian elimination.
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
% parameter..: scalar(float) TestStatistic         Distance in the uniform norm between the empirical cumulative distribution and the hypothesized
% parameter..: scalar(int)   SampleSize            Number of samples in the empirical cumulative distribution
% return.....: scalar(float) CumulativeProbability Cumulative probability of the test statistic
%
function CumulativeProbability = ExactKolmogorov(TestStatistic, SampleSize)
  
  % Escape with bad data
  if ~ (isnumeric(TestStatistic) && isnumeric(SampleSize))
    sprintf('%s', 'Data must be numeric.')
    CumulativeProbability = [];
    return
  end

  % Coerece the sample size and test statistic
  SampleSize    = sum(ceil(abs(SampleSize   (~ (isnan(SampleSize)    | isinf(SampleSize))))));
  TestStatistic = sum     (abs(TestStatistic(~ (isnan(TestStatistic) | isinf(TestStatistic)))));

  % There are no samples
  if SampleSize < 1
    CumulativeProbability = 1;

  % The statistic is too small
  elseif TestStatistic <= (1 ./ (2 .* SampleSize))
    CumulativeProbability = 0;
    
  % The statistic is too large
  elseif TestStatistic >= 1  
    CumulativeProbability = 1;
  
  % Sample size is one, case is linear
  elseif SampleSize < 2
    CumulativeProbability = 2 .* TestStatistic - 1;
  
  % For small test statistics use the power formula
  elseif TestStatistic <= 1 ./ SampleSize
    CumulativeProbability = prod((TestStatistic - 1 ./ (2 .* SampleSize)) .* (1:SampleSize));
      
  % For small samples find the complete matrix determinant
  elseif SampleSize < 128

    % Loop intialization, preallocate the matrix
    InclusionExclusion = zeros(SampleSize, SampleSize);  
          
    % Loop through the indices of the inclusion-exclusion matrix on the first non-zero diagonal
    for row = 2:SampleSize,
      InclusionExclusion(row, row - 1) = row;
    end
    
    % Loop through the upper triangle
    for row = 1:SampleSize,        
      for column = row:SampleSize,
        
        % Upper bound
        upper = min                               ...
        (                                         ...
          1,                                      ...
          (row - 1) ./ SampleSize + TestStatistic ...
        );
      
        % Lower bound
        lower = max                            ...
        (                                      ...
          0,                                   ...
          column ./ SampleSize - TestStatistic ...
        );
      
        % Bounds are ordered properly
        if upper > lower
          InclusionExclusion(row, column) = row .* prod                            ...
          (                                                                        ...
            (upper - lower) .* ones(1, column - row + 1) ./ (1:(column - row + 1)) ...
          );
          
        % The bounds do not agree  
        else
          break
        end
      end
    end
    
    % Find the determinant of the inclusion-exclusion matrix
    CumulativeProbability=max   ...
    (                           ...
      0,                        ...
      min                       ...
      (                         ...
        1,                      ...
        det(InclusionExclusion) ...
      )                         ...
    );
    
  % For large samples use stepwise gaussian elimination
  else
    
    % Loop initialization
    RowStack = zeros(1, SampleSize);
      
    % Build the first row  
    for column = 1:SampleSize,
      
      % Upper bound
      upper = min     ...
      (               ...
        1,            ...
        TestStatistic ...
      );
      
      % Lower bound
      lower = max                            ...
      (                                      ...
        0,                                   ...
        column ./ SampleSize - TestStatistic ...
      );
      
      % Bounds are properly ordered
      if upper > lower
        RowStack(1, column) = prod                           ...
        (                                                    ...
          (upper - lower) .* ones(1, column) ./ (1:(column)) ...
        );
        
      % The bounds do not agree  
      else
        break
      end
    end
    
    % Assign the probability
    InclusionExclusion(1, 1) = RowStack(1, 1);
    
    % Loop through the rows, Gaussian elimanting from the previous row
    for row = 2:SampleSize,
      for column = row:SampleSize,
        
        % Upper Bound
        upper = min                               ...
        (                                         ...
          1,                                      ...
          (row - 1) ./ SampleSize + TestStatistic ...
        );
        
        % Lower bound
        lower = max                            ...
        (                                      ...
          0,                                   ...
          column ./ SampleSize - TestStatistic ...
        );
        
        % Subtraction allowed, calculate the entry
        if upper > lower
          RowStack(1, column) = row .*                                               ...
          (                                                                          ...
            prod                                                                     ...
            (                                                                        ...
              (upper - lower) .* ones(1, column - row + 1) ./ (1:(column - row + 1)) ...
            ) -                                                                      ...
            RowStack(1, column) ./ RowStack(1, row - 1)                              ...
          );
        
        % The bounds do not agree
        else
          break
        end
      end
      RowStack(1, row - 1) = 0;
      
      % Assign the eigen value
      InclusionExclusion(row, 1)= RowStack(1, row);
    end

    % Balanced addition, drop the imaginary part because the real part becomes the magnitude
    CumulativeProbability = max  ...
    (                            ...
      0,                         ...
      min                        ...
      (                          ...
        1,                       ...
        prod(InclusionExclusion) ...
      )                          ...
    );
  end
end
