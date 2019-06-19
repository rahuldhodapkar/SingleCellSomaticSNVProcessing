CREATE TABLE variant (
	"DonorNum" DECIMAL NOT NULL, /*MOD*/
	"TissueType" VARCHAR NOT NULL, /*MOD*/
	"CellName" VARCHAR NOT NULL, /*MOD*/
	"ChromKey" DECIMAL, 
	"CHROM" VARCHAR, /*MOD*/ 
	"POS" DECIMAL, 
	"ID" BOOLEAN, 
	"REF" VARCHAR, 
	"ALT" VARCHAR, 
	"QUAL" DECIMAL, 
	"FILTER" BOOLEAN, 
	"AC" VARCHAR, /*MOD*/
	"AF" VARCHAR, /*MOD*/
	"AN" DECIMAL, 
	"BaseQRankSum" DECIMAL, 
	"ClippingRankSum" BOOLEAN, 
	"DP" DECIMAL, 
	"DS" BOOLEAN, 
	"ExcessHet" DECIMAL, 
	"FS" DECIMAL, 
	"HaplotypeScore" BOOLEAN, 
	"InbreedingCoeff" BOOLEAN, 
	"MLEAC" VARCHAR, /*MOD*/ 
	"MLEAF" VARCHAR, /*MOD*/ 
	"MQ" DECIMAL, 
	"MQRankSum" BOOLEAN, 
	"QD" DECIMAL, 
	"ReadPosRankSum" DECIMAL, 
	"SOR" DECIMAL, 
	"Indiv" VARCHAR, 
	"gt_AD" VARCHAR, /*MOD*/
	"gt_DP" DECIMAL, 
	"gt_GQ" DECIMAL, 
	"gt_GT" VARCHAR, 
	"gt_PL" VARCHAR, /*MOD*/
	"gt_GT_alleles" VARCHAR
);
