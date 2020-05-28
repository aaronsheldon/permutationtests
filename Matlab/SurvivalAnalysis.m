% Author.....: Aaron Sheldon
% Date.......: July, 2010
% Description: Initial data analysis of survival data

% Intent to treat, syncope events only, placebo estimator % % % % % % % % %
IntentSyncopePlaceboEstimator = KaplanMeierEstimator                     ...
(                                                                        ...
  IntentSyncope(IntentSyncope(:, 3) == 0 & IntentSyncope(:, 4) == 1, 2), ...
  IntentSyncope(IntentSyncope(:, 3) == 1 & IntentSyncope(:, 4) == 1, 2)  ...
);

% Intent to treat, syncope events only, drug estimator
IntentSyncopeDrugEstimator = KaplanMeierEstimator                        ...
(                                                                        ...
  IntentSyncope(IntentSyncope(:, 3) == 0 & IntentSyncope(:, 4) == 0, 2), ...
  IntentSyncope(IntentSyncope(:, 3) == 1 & IntentSyncope(:, 4) == 0, 2)  ...
);

% Intent to treat, syncope events only, statistic
IntentSyncopeStatistic = KaplanMeierStatistic                            ...
(                                                                        ...
  IntentSyncope(IntentSyncope(:, 3) == 0 & IntentSyncope(:, 4) == 1, 2), ...
  IntentSyncope(IntentSyncope(:, 3) == 1 & IntentSyncope(:, 4) == 1, 2), ...
  IntentSyncope(IntentSyncope(:, 3) == 0 & IntentSyncope(:, 4) == 0, 2), ...
  IntentSyncope(IntentSyncope(:, 3) == 1 & IntentSyncope(:, 4) == 0, 2)  ...
);

% Intent to treat, syncope and failure events, placebo estimator % % % % %
IntentSyncopeFailurePlaceboEstimator = KaplanMeierEstimator                                   ...
(                                                                                             ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 0 & IntentSyncopeFailure(:, 4) == 1, 2), ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 1 & IntentSyncopeFailure(:, 4) == 1, 2)  ...
);

% Intent to treat, syncope and failure events, drug estimator
IntentSyncopeFailureDrugEstimator = KaplanMeierEstimator                                      ...
(                                                                                             ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 0 & IntentSyncopeFailure(:, 4) == 0, 2), ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 1 & IntentSyncopeFailure(:, 4) == 0, 2)  ...
);

% Intent to treat, syncope and failure events, statistic
IntentSyncopeFailureStatistic = KaplanMeierStatistic                                          ...
(                                                                                             ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 0 & IntentSyncopeFailure(:, 4) == 1, 2), ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 1 & IntentSyncopeFailure(:, 4) == 1, 2), ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 0 & IntentSyncopeFailure(:, 4) == 0, 2), ...
  IntentSyncopeFailure(IntentSyncopeFailure(:, 3) == 1 & IntentSyncopeFailure(:, 4) == 0, 2)  ...
);

% As treated, syncope events only, placebo estimator % % % % % % % % % % %
TreatedSyncopePlaceboEstimator = KaplanMeierEstimator                       ...
(                                                                           ...
  TreatedSyncope(TreatedSyncope(:, 3) == 0 & TreatedSyncope(:, 4) == 1, 2), ...
  TreatedSyncope(TreatedSyncope(:, 3) == 1 & TreatedSyncope(:, 4) == 1, 2)  ...
);

% As treated, syncope events only, drug estimator
TreatedSyncopeDrugEstimator = KaplanMeierEstimator                          ...
(                                                                           ...
  TreatedSyncope(TreatedSyncope(:, 3) == 0 & TreatedSyncope(:, 4) == 0, 2), ...
  TreatedSyncope(TreatedSyncope(:, 3) == 1 & TreatedSyncope(:, 4) == 0, 2)  ...
);

% As treated, syncope events only, statistic
TreatedSyncopeStatistic = KaplanMeierStatistic                              ...
(                                                                           ...
  TreatedSyncope(TreatedSyncope(:, 3) == 0 & TreatedSyncope(:, 4) == 1, 2), ...
  TreatedSyncope(TreatedSyncope(:, 3) == 1 & TreatedSyncope(:, 4) == 1, 2), ...
  TreatedSyncope(TreatedSyncope(:, 3) == 0 & TreatedSyncope(:, 4) == 0, 2), ...
  TreatedSyncope(TreatedSyncope(:, 3) == 1 & TreatedSyncope(:, 4) == 0, 2)  ...
);

% As treated, syncope and failure events, placebo estimator % % % % % % % %
TreatedSyncopeFailurePlaceboEstimator = KaplanMeierEstimator                                     ...
(                                                                                                ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 0 & TreatedSyncopeFailure(:, 4) == 1, 2), ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 1 & TreatedSyncopeFailure(:, 4) == 1, 2)  ...
);

% As treated, syncope and failure events, drug estimator
TreatedSyncopeFailureDrugEstimator = KaplanMeierEstimator                                        ...
(                                                                                                ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 0 & TreatedSyncopeFailure(:, 4) == 0, 2), ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 1 & TreatedSyncopeFailure(:, 4) == 0, 2)  ...
);

% As treated, syncope and failure events, statistic
TreatedSyncopeFailureStatistic = KaplanMeierStatistic                                            ...
(                                                                                                ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 0 & TreatedSyncopeFailure(:, 4) == 1, 2), ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 1 & TreatedSyncopeFailure(:, 4) == 1, 2), ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 0 & TreatedSyncopeFailure(:, 4) == 0, 2), ...
  TreatedSyncopeFailure(TreatedSyncopeFailure(:, 3) == 1 & TreatedSyncopeFailure(:, 4) == 0, 2)  ...
);