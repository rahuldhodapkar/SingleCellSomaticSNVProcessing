CREATE TABLE donor (
	"﻿SampleCode" VARCHAR NOT NULL, 
	"Sample" DECIMAL NOT NULL, 
	"Group" VARCHAR NOT NULL, 
	"LastDMTNotes" VARCHAR, 
	"AgeAtProcedure" DECIMAL NOT NULL, 
	"DateOfProcedure" DATE NOT NULL, 
	"DiagnosisYear" DECIMAL, 
	"YearsSinceDiagnosis" DECIMAL, 
	"cmHeight" DECIMAL NOT NULL, 
	"kgWeight" DECIMAL NOT NULL, 
	"BMI" DECIMAL NOT NULL, 
	"Gender" VARCHAR NOT NULL, 
	"Ethnicity" VARCHAR NOT NULL, 
	"scRNAseqDate" DATE NOT NULL, 
	"SequencedDate" DATE NOT NULL
);
