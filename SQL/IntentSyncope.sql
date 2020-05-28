/*
 * Author.....: Aaron Sheldon
 * Date.......: June, 2010
 * Description: List first outcome or censor event for
 *              each patient, for intent to treat analysis,
 *              excluding treatment failures.
 */
ALTER PROCEDURE IntentSyncope AS

  /* Start the query */
  SELECT 

    /* Return the patient study number for auditing */
    P.PID AS StudyNumber,

    /* Time to first event */
    CASE

      /* There was no outcome, return the time to censor, cut off at a year */ 
      WHEN MIN(O.EventDate) IS NULL AND MAX(RC.EventDate) IS NULL THEN
        365

      /* There was no outcome, return the time to censor, cut of at one day */ 
      WHEN MIN(O.EventDate) IS NULL AND (MAX(RC.EventDate) - MIN(LC.EventDate)) < 1 THEN
        1

      /* There was no outcome, return the time to censor */ 
      WHEN MIN(O.EventDate) IS NULL THEN
        (MAX(RC.EventDate) - MIN(LC.EventDate))

      /* The outcome occured on the first day, lower bound with one day */
      WHEN (MIN(O.EventDate) - MIN(LC.EventDate)) < 1 THEN
        1

      /* There was an outcome, return the time to event */
      ELSE
        (MIN(O.EventDate) - MIN(LC.EventDate))      
    END AS EventInterval,

    /* Type of event */
    CASE

      /* There was no outcome, return the censor */ 
      WHEN MIN(O.EventDate) IS NULL THEN
        1

      /* There was an outcome */
      ELSE
        0
    END AS CensorFlag,

    /* 1=On placebo, 0=on drug*/
    MIN(CAST(P.Placebo AS tinyint)) AS PlaceboFlag,

    /* 1=female, 0=male */
    MIN(P.GenderId - 1) AS FemaleFlag,

    /* Right censor date */
    MIN(LC.EventDate) AS LeftCensorDate, 

    /* Outcome date */
    MIN(O.EventDate) AS OutcomeDate, 

    /* Left censor date */
    MIN(RC.EventDate) AS RightCensorDate

    /* Master patient record */
    FROM Patient AS P

    /* Left censor, proxy the first date of observation, excluding primary outcome */
    LEFT OUTER JOIN
    (
      
      /* Find the minimum of all observation dates not including outcome */
      SELECT G.Id AS Id, MIN(G.EventDate) AS EventDate
        FROM 
        (

          /* Study drug prescription date */
          SELECT D.PatientId AS Id, D.StartDate AS EventDate 
            FROM Dose AS D
          UNION ALL

          /* Drug compliance handout date */
          SELECT B.PatientId, B.DispenseDate
            FROM Bottle AS B
          UNION ALL

          /* Drug compliance return date */
          SELECT B.PatientId, B.ReturnDate
            FROM Bottle AS B
          UNION ALL

          /* First baseline interview */
          SELECT P.Id, P.FormOneDate
            FROM Patient AS P
          UNION ALL

          /* Second baseline interview */
          SELECT P.Id, P.FormTwoDate
            FROM Patient AS P
          UNION ALL

          /* Study observation date */
          SELECT F.PatientId, F.InterviewDate
            FROM FollowUp AS F
          UNION ALL 

          /* Study  withdrawal date */
          SELECT F.PatientId, F.WithdrawDate
            FROM FollowUp AS F
        ) AS G

      /* Only return a single record for each patient */
      GROUP BY G.Id
    ) AS LC

    /* Link the left censor data to the master patient data */
    ON P.Id = LC.Id

    /* Right censor, when present the first withdraw, else the last observation, truncated at one year */
    LEFT OUTER JOIN
    (

      /* Return either the first withdraw or the last observation including outcome */  
      SELECT P.Id AS Id, COALESCE(MIN(F.WithdrawDate), MAX(G.EventDate)) AS EventDate
        FROM Patient AS P

        /* Withdraw records */
        LEFT OUTER JOIN FollowUp AS F
          ON P.Id = F.PatientId

        /* All observations excluding withdraw */
        LEFT OUTER JOIN
        (

          /* Study drug prescription date */
          SELECT D.PatientId AS Id, D.StartDate AS EventDate 
            FROM Dose AS D
          UNION ALL

          /* Drug compliance handout date */
          SELECT B.PatientId, B.DispenseDate
            FROM Bottle AS B
          UNION ALL

          /* Drug compliance return date */
          SELECT B.PatientId, B.ReturnDate
            FROM Bottle AS B
          UNION ALL

          /* First baseline interview */
          SELECT P.Id, P.FormOneDate
            FROM Patient AS P
          UNION ALL

          /* Second baseline interview */
          SELECT P.Id, P.FormTwoDate
             FROM Patient AS P
          UNION ALL

          /* Study observation date */
          SELECT F.PatientId, F.InterviewDate
            FROM FollowUp AS F
          UNION ALL 

          /* Study  outcome date */
          SELECT S.PatientId, S.EpisodeDate
            FROM Syncope AS S
        ) AS G

        /* Link to the patient */
        ON P.Id = G.Id        

        /* Only return a single record for each patient */
        GROUP BY P.Id
    ) AS RC

    /* Link the left censor data to master patient data, dates must be within first year of right censor date */
    ON P.Id = RC.Id
      AND FLOOR(RC.EventDate - LC.EventDate) >= 0
      AND FLOOR(RC.EventDate - LC.EventDate) <= 365

    /* Primary outcome only, syncope */
    LEFT OUTER JOIN
    ( 

      /* Return all observed spells*/
      SELECT S.PatientId AS Id, S.EpisodeDate AS EventDate
        FROM Syncope AS S
    ) AS O

    /* Link the outcome to the master date, dates must fall between right and left censor dates */
    ON P.Id = O.Id
      AND FLOOR(O.EventDate - LC.EventDate) >= 0
      AND FLOOR(O.EventDate - LC.EventDate) <= 365
      AND (RC.EventDate IS NULL OR FLOOR(RC.EventDate - O.EventDate) >= 0)

    /* Only return a single record for each patient */
    GROUP BY P.PID
    ORDER BY P.PID

/* Exit */
RETURN