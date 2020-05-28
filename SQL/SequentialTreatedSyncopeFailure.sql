/**
 * Author.....: Aaron Sheldon
 * Date.......: June, 2010
 * Description: List every outcome or censor event for
 *              each patient, for as treated analysis,
 *              including treatment failures.
 */
ALTER PROCEDURE SequentialTreatedSyncopeFailure AS

  /* Start the query */
  SELECT 

    /* Return the patient study number for auditing */
    P.PID AS StudyNumber,

    /* Sequences of study events */
    CASE

      /* Put a break between patients */
      WHEN O.EventFlag = 0 THEN
        NULL

      /* First date is left censor date */
      WHEN O.EventFlag = 1 THEN
        LC.EventDate

      /* Sequence of outcome dates */
      WHEN O.EventFlag = 2 THEN
        O.EventDate

      /* Sequence of event failures */
      WHEN O.EventFlag = 2.1 THEN
        O.EventDate

      /* Right censor exceeded one year */
      WHEN RC.EventDate IS NULL THEN
        FLOOR(LC.EventDate + 365)

      /* Last event is right censor date */
      ELSE
        RC.EventDate      
    END AS EventDate,

    /* Type of event */
    O.EventFlag AS EventFlag,

    /* 1=On placebo, 0=on drug*/
    CAST(P.Placebo AS tinyint) AS PlaceboFlag,

    /* 1=female, 0=male */
    (P.GenderId - 1) AS FemaleFlag,

    /* Right censor date */
    LC.EventDate AS LeftCensorDate, 

    /* Outcome date */
    O.EventDate AS OutcomeDate, 

    /* Left censor date */
    RC.EventDate AS RightCensorDate

    /* Master patient record */
    FROM Patient AS P

    /* Left censor, first date of medication consultation */
    LEFT OUTER JOIN
    (
      
      /* Find the minimum of all medication consultations */
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
        ) AS G

      /* Only return a single record for each patient */
      GROUP BY G.Id
    ) AS LC

    /* Link the left censor data to the master patient data */
    ON P.Id = LC.Id

    /* Right censor, when present the first withdraw, else the last medication consultation, truncated at one year */
    LEFT OUTER JOIN
    (

      /* Return either the first withdraw or the last medication consultation */  
      SELECT P.Id AS Id, COALESCE(MIN(F.WithdrawDate), MAX(G.EventDate)) AS EventDate
        FROM Patient AS P

        /* Any withdraw record is a right censor event */
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

    /* Primary and secondary outcome, syncope and treatment failure */
    LEFT OUTER JOIN
    ( 

      /* Return all observed spells*/
      SELECT S.PatientId AS Id, S.EpisodeDate AS EventDate, 2 AS EventFlag
        FROM Syncope AS S
      UNION ALL

      /* Return first treatment failure events */
      SELECT F.PatientId, MIN(F.WithdrawDate), 2.1
        FROM FollowUp AS F
        WHERE F.WithdrawComments IS NOT NULL
          AND UPPER(F.WithdrawComments) NOT LIKE '%PACEMAKER%'
          AND UPPER(F.WithdrawComments) NOT LIKE '%EPILEPSY%'
          AND UPPER(F.WithdrawComments) NOT LIKE '%CANCER%'
          AND UPPER(F.WithdrawComments) NOT LIKE '%PREGNAN%'
        GROUP BY F.PatientId
      UNION ALL

      /* Return a record to break between patients */
      SELECT P.Id, NULL, 0
        FROM Patient AS P
      UNION ALL

      /* Return a record to store the left censor date */
      SELECT P.Id, NULL, 1
        FROM Patient AS P
      UNION ALL

      /* Return a record to store the right censor date */
      SELECT P.Id, NULL, 3
        FROM Patient AS P  
    ) AS O

    /* Link the outcome to the master date, dates must fall between right and left censor dates */
    ON P.Id = O.Id
      AND (O.EventDate IS NULL OR FLOOR(O.EventDate - LC.EventDate) >= 0)
      AND (O.EventDate IS NULL OR FLOOR(O.EventDate - LC.EventDate) <= 365)
      AND (O.EventDate IS NULL OR RC.EventDate IS NULL OR FLOOR(RC.EventDate - O.EventDate) >= 0)

    /* Order by patient then event sequence */
    ORDER BY P.PID, EventDate, O.EventFlag

/* Exit */
RETURN
