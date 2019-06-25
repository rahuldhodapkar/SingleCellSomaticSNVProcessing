

/* Get frequency of each "Top Consequence" and plot */
select "Top Consequence", count(*),
         repeat('■', (  100.0
                      * count(*)
                      / sum(count(*)) over()
                     )::integer
               ) as pct
    from annotation
group by "Top Consequence"
order by "Top Consequence";


/* Get frequency of each "Top Consequence" and plot */
select "DonorNum", count(*),
         repeat('■', (  100.0
                      * count(*)
                      / sum(count(*)) over()
                     )::integer
               ) as pct
    from variant
group by "DonorNum"
order by "DonorNum";



/* Get frequency of each "Top Consequence" and plot */
select "CellName", count(*),
         repeat('■', (  100.0
                      * count(*)
                      / sum(count(*)) over()
                     )::integer
               ) as pct
    from variant
group by "CellName"
order by "CellName";


/* Get number of cells from blood for MS and HC */
SELECT 
    "Group",
    SUM("n") as NumCells
FROM 
    (SELECT 
        "DonorNum", COUNT(DISTINCT "CellName") as n 
    FROM 
        variant 
    GROUP BY 
        "DonorNum", "TissueType" 
    HAVING 
        "TissueType"='B') as a
INNER JOIN 
    donor ON a."DonorNum"=donor."Sample"
GROUP BY
    "Group";


/* Extract all raw variant calls for a specific gene */

SELECT * FROM
  (SELECT * FROM
    (SELECT * FROM annotation
      WHERE
          "Gene"='GLYCTK') as a
    INNER JOIN hq_variant
    USING
      ("CHROM", "POS", "REF", "ALT")) as k
INNER JOIN
  donor
ON
  k."DonorNum" = donor."Sample";



/* Get number of cells with ERCC transcripts by donor and tissue type */
select SUBSTRING("Gene", 1, 4) as "Code", "DonorNum", "TissueType", "Group",
        count(*),
         repeat('■', (  100.0
                      * count(*)
                      / sum(count(*)) over()
                     )::integer
               ) as pct
    from expression_donor
group by SUBSTRING("Gene", 1, 4), "DonorNum", "TissueType", "Group"
HAVING
  SUBSTRING("Gene", 1, 4) = 'ERCC';

/* list all materialized views */
select schemaname as schema_name,
       matviewname as view_name
from pg_matviews
order by schema_name,
         view_name;

/* Get number of cells from blood for MS and HC, excluding 9 + 10 */
SELECT 
    "Group",
    SUM("n") as NumCells
FROM 
    (SELECT 
        "DonorNum", COUNT(DISTINCT "CellName") as n 
    FROM 
        (SELECT "DonorNum", "CellName", "TissueType"
          FROM variant WHERE "DonorNum" NOT IN (9,10)) as clean_vars
    GROUP BY 
        "DonorNum", "TissueType" 
    HAVING 
        "TissueType"='B') as a
INNER JOIN 
    donor ON a."DonorNum"=donor."Sample"
GROUP BY
    "Group";



/* Get number of cells with FOXP3 transcripts by donor and tissue type */
select "Gene", "DonorNum", "TissueType", "Group",
        count(*),
         repeat('■', (  100.0
                      * count(*)
                      / sum(count(*)) over()
                     )::integer
               ) as pct
    from expression_donor_unclean
group by "Gene", "DonorNum", "TissueType", "Group"
HAVING
  "Gene" = 'RPL28' AND "TissueType"='B'
ORDER BY
  "Group";

/* list all materialized views */
select schemaname as schema_name,
       matviewname as view_name
from pg_matviews
order by schema_name,
         view_name;




