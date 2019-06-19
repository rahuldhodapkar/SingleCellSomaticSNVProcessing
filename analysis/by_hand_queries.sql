

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