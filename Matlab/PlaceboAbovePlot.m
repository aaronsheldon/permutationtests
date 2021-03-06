% Author.....: Aaron Sheldon
% Date.......: July, 2010
% Description: Plot Monte Carlo results of survival data for placebo above
%              on drug curve

% Intent to Treat and syncope only % % % % % % % % % % % % % % % % % % % %
subplot(2, 2, 1);
hold;

% Individual Monte Carlo points
scatter                          ...
(                                ...
  IntentSyncopeMonteCarlo(:, 1), ...
  IntentSyncopeMonteCarlo(:, 3), ...
  4,                             ...
  '+'                            ...
);

% Left bound extreme region
plot                              ...
(                                 ...
  [                               ...
    IntentSyncopeStatistic(1, 1); ...
    IntentSyncopeStatistic(1, 1)  ...
  ],                              ... 
  [                               ...
    IntentSyncopeStatistic(1, 3); ...
    1                             ...
  ],                              ...
  '--k'                           ...
);

% Lower bound extreme region
plot                              ...
(                                 ...
  [                               ...
    IntentSyncopeStatistic(1, 1); ...
    1                             ...
  ],                              ...
  [                               ...
    IntentSyncopeStatistic(1, 3); ...
    IntentSyncopeStatistic(1, 3)  ...
  ],                              ...
  '--k'                           ...
);

% As Treated and syncope only % % % % % % % % % % % % % % % % % % % % % % %
subplot(2, 2, 2);
hold;

% Individual Monte Carlo points
scatter                           ...
(                                 ...
  TreatedSyncopeMonteCarlo(:, 1), ...
  TreatedSyncopeMonteCarlo(:, 3), ...
  4,                              ...
  '+'                             ...
);

% Left bound extreme region
plot                               ...
(                                  ...
  [                                ...
    TreatedSyncopeStatistic(1, 1); ...
    TreatedSyncopeStatistic(1, 1)  ...
  ],                               ... 
  [                                ...
    TreatedSyncopeStatistic(1, 3); ...
    1                              ...
  ],                               ...
  '--k'                            ...
);

% Lower bound extreme region
plot                               ...
(                                  ...
  [                                ...
    TreatedSyncopeStatistic(1, 1); ...
    1                              ...
  ],                               ...
  [                                ...
    TreatedSyncopeStatistic(1, 3); ...
    TreatedSyncopeStatistic(1, 3)  ...
  ],                               ...
  '--k'                            ...
);

% Intent to Treat, syncope and treatment failure% % % % % % % % % % % % % %
subplot(2, 2, 3);
hold;

% Individual Monte Carlo points
scatter                                 ...
(                                       ...
  IntentSyncopeFailureMonteCarlo(:, 1), ...
  IntentSyncopeFailureMonteCarlo(:, 3), ...
  4,                                    ...
  '+'                                   ...
);

% Left bound extreme region
plot                                     ...
(                                        ...
  [                                      ...
    IntentSyncopeFailureStatistic(1, 1); ...
    IntentSyncopeFailureStatistic(1, 1)  ...
  ],                                     ... 
  [                                      ...
    IntentSyncopeFailureStatistic(1, 3); ...
    1                                    ...
  ],                                     ...
  '--k'                                  ...
);

% Lower bound extreme region
plot                                     ...
(                                        ...
  [                                      ...
    IntentSyncopeFailureStatistic(1, 1); ...
    1                                    ...
  ],                                     ...
  [                                      ...
    IntentSyncopeFailureStatistic(1, 3); ...
    IntentSyncopeFailureStatistic(1, 3)  ...
  ],                                     ...
  '--k'                                  ...
);

% As Treated, syncope and treatment failure % % % % % % % % % % % % % % % %
subplot(2, 2, 4);
hold;

% Individual Monte Carlo points
scatter                                  ...
(                                        ...
  TreatedSyncopeFailureMonteCarlo(:, 1), ...
  TreatedSyncopeFailureMonteCarlo(:, 3), ...
  4,                                     ...
  '+'                                    ...
);

% Left bound extreme region
plot                                      ...
(                                         ...
  [                                       ...
    TreatedSyncopeFailureStatistic(1, 1); ...
    TreatedSyncopeFailureStatistic(1, 1)  ...
  ],                                      ... 
  [                                       ...
    TreatedSyncopeFailureStatistic(1, 3); ...
    1                                     ...
  ],                                      ...
  '--k'                                   ...
);

% Lower bound extreme region
plot                                      ...
(                                         ...
  [                                       ...
    TreatedSyncopeFailureStatistic(1, 1); ...
    1                                     ...
  ],                                      ...
  [                                       ...
    TreatedSyncopeFailureStatistic(1, 3); ...
    TreatedSyncopeFailureStatistic(1, 3)  ...
  ],                                      ...
  '--k'                                   ...
);