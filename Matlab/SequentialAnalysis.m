% Author.....: Aaron Sheldon
% Date.......: July, 2010
% Description: Initial data analysis of sequential survival data

% Intent to treat, syncope events only, placebo estimator % % % % % % % % %
IntentSyncopeSequencePlaceboEstimator = SequentialKaplanMeierEstimator ...
(                                                                      ...
  IntentSyncopeSequence(IntentSyncopeSequence(:, 4) == 1, 2)           ...
);

% Intent to treat, syncope events only, drug estimator
IntentSyncopeSequenceDrugEstimator = SequentialKaplanMeierEstimator ...
(                                                                   ...
  IntentSyncopeSequence(IntentSyncopeSequence(:, 4) == 0, 2)        ...
);

% Intent to treat, syncope events only, statistic
IntentSyncopeSequenceStatistic = SequentialKaplanMeierStatistic ...
(                                                               ...
  IntentSyncopeSequence(IntentSyncopeSequence(:, 4) == 1, 2),   ...
  IntentSyncopeSequence(IntentSyncopeSequence(:, 4) == 0, 2)    ...
);

% Intent to treat, syncope and failure events, placebo estimator % % % % %
IntentSyncopeFailureSequencePlaceboEstimator = SequentialKaplanMeierEstimator ...
(                                                                             ...
  IntentSyncopeFailureSequence(IntentSyncopeFailureSequence(:, 4) == 1, 2)    ...
);

% Intent to treat, syncope and failure events, drug estimator
IntentSyncopeFailureSequenceDrugEstimator = SequentialKaplanMeierEstimator ...
(                                                                          ...
  IntentSyncopeFailureSequence(IntentSyncopeFailureSequence(:, 4) == 0, 2) ...
);

% Intent to treat, syncope and failure events, statistic
IntentSyncopeFailureSequenceStatistic = SequentialKaplanMeierStatistic      ...
(                                                                           ...
  IntentSyncopeFailureSequence(IntentSyncopeFailureSequence(:, 4) == 1, 2), ...
  IntentSyncopeFailureSequence(IntentSyncopeFailureSequence(:, 4) == 0, 2)  ...
);

% As treated, syncope events only, placebo estimator % % % % % % % % % % %
TreatedSyncopeSequencePlaceboEstimator = SequentialKaplanMeierEstimator ...
(                                                                       ...
  TreatedSyncopeSequence(TreatedSyncopeSequence(:, 4) == 1, 2)          ...
);

% As treated, syncope events only, drug estimator
TreatedSyncopeSequenceDrugEstimator = SequentialKaplanMeierEstimator ...
(                                                                    ...
  TreatedSyncopeSequence(TreatedSyncopeSequence(:, 4) == 0, 2)       ...
);

% As treated, syncope events only, statistic
TreatedSyncopeSequenceStatistic = SequentialKaplanMeierStatistic ...
(                                                                ...
  TreatedSyncopeSequence(TreatedSyncopeSequence(:, 4) == 1, 2),  ...
  TreatedSyncopeSequence(TreatedSyncopeSequence(:, 4) == 0, 2)   ...
);

% As treated, syncope and failure events, placebo estimator % % % % % % % %
TreatedSyncopeFailureSequencePlaceboEstimator = SequentialKaplanMeierEstimator ...
(                                                                              ...
  TreatedSyncopeFailureSequence(TreatedSyncopeFailureSequence(:, 4) == 1, 2)   ...
);

% As treated, syncope and failure events, drug estimator
TreatedSyncopeFailureSequenceDrugEstimator = SequentialKaplanMeierEstimator  ...
(                                                                            ...
  TreatedSyncopeFailureSequence(TreatedSyncopeFailureSequence(:, 4) == 0, 2) ...
);

% As treated, syncope and failure events, statistic
TreatedSyncopeFailureSequenceStatistic = SequentialKaplanMeierStatistic       ...
(                                                                             ...
  TreatedSyncopeFailureSequence(TreatedSyncopeFailureSequence(:, 4) == 1, 2), ...
  TreatedSyncopeFailureSequence(TreatedSyncopeFailureSequence(:, 4) == 0, 2)  ...
);