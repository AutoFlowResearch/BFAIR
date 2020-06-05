-- Table: IsaTab.ONTOLOGY_SOURCE_REFERENCE

-- DROP TABLE "IsaTab"."ONTOLOGY_SOURCE_REFERENCE";

CREATE TABLE "IsaTab"."ONTOLOGY_SOURCE_REFERENCE"
(
    "ID" oid,
    "Term_Source_Name" text COLLATE pg_catalog."default",
    "Term_Source_File" text COLLATE pg_catalog."default",
    "Term_Source_Version" text COLLATE pg_catalog."default",
    "Term_Source_Description" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."ONTOLOGY_SOURCE_REFERENCE"
    OWNER to postgres;

-- Table: IsaTab.INVESTIGATION

-- DROP TABLE "IsaTab"."INVESTIGATION";

CREATE TABLE "IsaTab"."INVESTIGATION"
(
    "ID" oid,
    "Investigation_Identifier" text COLLATE pg_catalog."default" NOT NULL,
    "Investigation_Title" text COLLATE pg_catalog."default",
    "Investigation_Description" text COLLATE pg_catalog."default",
    "Investigation_Submission_date" date,
    "Investigation_Public_Release_Date" date,
    "Investigation_Publications" text COLLATE pg_catalog."default",
    "Investigation_Contacts" text COLLATE pg_catalog."default",
    CONSTRAINT investigation_pkey PRIMARY KEY ("Investigation_Identifier")
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."INVESTIGATION"
    OWNER to postgres;

-- Table: IsaTab.INVESTIGATION_CONTACTS

-- DROP TABLE "IsaTab"."INVESTIGATION_CONTACTS";

CREATE TABLE "IsaTab"."INVESTIGATION_CONTACTS"
(
    "ID" oid,
    "Investigation_Person_Last_Name" text COLLATE pg_catalog."default",
    "Investigation_Person_First_Name" text COLLATE pg_catalog."default",
    "Investigation_Person_Mid_Initials" character varying(5) COLLATE pg_catalog."default",
    "Investigation_Person_Email" text COLLATE pg_catalog."default",
    "Investigation_Person_Phone" text COLLATE pg_catalog."default",
    "Investigation_Person_Fax" text COLLATE pg_catalog."default",
    "Investigation_Person_Address" text COLLATE pg_catalog."default",
    "Investigation_Person_Affiliation" text COLLATE pg_catalog."default",
    "Investigation_Person_Roles" text COLLATE pg_catalog."default",
    "Investigation_Person_Roles_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Investigation_Person_Roles_Term_Source_REF" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."INVESTIGATION_CONTACTS"
    OWNER to postgres;

-- Table: IsaTab.INVESTIGATION_PUBLICATIONS

-- DROP TABLE "IsaTab"."INVESTIGATION_PUBLICATIONS";

CREATE TABLE "IsaTab"."INVESTIGATION_PUBLICATIONS"
(
    "ID" oid,
    "Investigation_PubMed_ID" text COLLATE pg_catalog."default",
    "Investigation_Publication_DOI" text COLLATE pg_catalog."default",
    "Investigation_Publication_Author_List" text COLLATE pg_catalog."default",
    "Investigation_Publication_Title" text COLLATE pg_catalog."default",
    "Investigation_Publication_Status" text COLLATE pg_catalog."default",
    "Investigation_Publication_Status_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Investigation_Publication_Status_Term_Source_REF" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."INVESTIGATION_PUBLICATIONS"
    OWNER to postgres;

-- Table: IsaTab.STUDY

-- DROP TABLE "IsaTab"."STUDY";

CREATE TABLE "IsaTab"."STUDY"
(
    "ID" oid,
    "Study_Identifier" text COLLATE pg_catalog."default",
    "Study_Title" text COLLATE pg_catalog."default",
    "Study_Description" text COLLATE pg_catalog."default",
    "Study_Submission_Date" date,
    "Study_Public_Release_Date" date,
    "Study_File_Name" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."STUDY"
    OWNER to postgres;

-- Table: IsaTab.STUDY_ASSAY

-- DROP TABLE "IsaTab"."STUDY_ASSAY";

CREATE TABLE "IsaTab"."STUDY_ASSAY"
(
    "ID" oid,
    "Study_Assay_File_Name" text COLLATE pg_catalog."default",
    "Study_Assay_Measurement_Type" text COLLATE pg_catalog."default",
    "Study_Assay_Measurement_Type_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Assay_Measurement_Type_Term_Source_REF" text COLLATE pg_catalog."default",
    "Study_Assay_Technology_Type" text COLLATE pg_catalog."default",
    "Study_Assay_Technology_Type_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Assay_Technology_Type_Term_Source_REF" text COLLATE pg_catalog."default",
    "Study_Assay_Technology_Platform" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."STUDY_ASSAY"
    OWNER to postgres;

-- Table: IsaTab.STUDY_CONTACTS

-- DROP TABLE "IsaTab"."STUDY_CONTACTS";

CREATE TABLE "IsaTab"."STUDY_CONTACTS"
(
    "ID" oid,
    "Study_Person_Last_Name" text COLLATE pg_catalog."default",
    "Study_Person_First_Name" text COLLATE pg_catalog."default",
    "Study_Person_Mid_Initials" character varying(5) COLLATE pg_catalog."default",
    "Study_Person_Email" text COLLATE pg_catalog."default",
    "Study_Person_Phone" text COLLATE pg_catalog."default",
    "Study_Person_Fax" text COLLATE pg_catalog."default",
    "Study_Person_Address" text COLLATE pg_catalog."default",
    "Study_Person_Affiliation" text COLLATE pg_catalog."default",
    "Study_Person_Roles" text COLLATE pg_catalog."default",
    "Study_Person_Roles_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Person_Roles_Term_Source_REF" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."STUDY_CONTACTS"
    OWNER to postgres

-- Table: IsaTab.STUDY_DESIGN_DESCRIPTORS

-- DROP TABLE "IsaTab"."STUDY_DESIGN_DESCRIPTORS";

CREATE TABLE "IsaTab"."STUDY_DESIGN_DESCRIPTORS"
(
    "ID" oid,
    "Study_Design_Type" text COLLATE pg_catalog."default",
    "Study_Design_Type_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Design_Type_Term_Source_REF" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."STUDY_DESIGN_DESCRIPTORS"
    OWNER to postgres;

-- Table: IsaTab.STUDY_FACTORS

-- DROP TABLE "IsaTab"."STUDY_FACTORS";

CREATE TABLE "IsaTab"."STUDY_FACTORS"
(
    "ID" oid,
    "Study_Factor_Name" text COLLATE pg_catalog."default",
    "Study_Factor_Type" text COLLATE pg_catalog."default",
    "Study_Factor_Type_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Factor_Type_Term_Source_REF" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."STUDY_FACTORS"
    OWNER to postgres;

-- Table: IsaTab.STUDY_PROTOCOLS

-- DROP TABLE "IsaTab"."STUDY_PROTOCOLS";

CREATE TABLE "IsaTab"."STUDY_PROTOCOLS"
(
    "ID" oid,
    "Study_Protocol_Name" text COLLATE pg_catalog."default",
    "Study_Protocol_Type" text COLLATE pg_catalog."default",
    "Study_Protocol_Type_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Protocol_Type_Term_Source_REF" text COLLATE pg_catalog."default",
    "Study_Protocol_Description" text COLLATE pg_catalog."default",
    "Study_Protocol_URI" text COLLATE pg_catalog."default",
    "Study_Protocol_Version" text COLLATE pg_catalog."default",
    "Study_Protocol_Parameters_Name" text COLLATE pg_catalog."default",
    "Study_Protocol_Parameter_Name_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Protocol_Parameter_Name_Term_Source_REF" text COLLATE pg_catalog."default",
    "Study_Protocol_Components_Name" text COLLATE pg_catalog."default",
    "Study_Protocol_Components_Type" text COLLATE pg_catalog."default",
    "Study_Protocol_Components_Type_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Protocol_Components_Type_Term_Source_REF" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."STUDY_PROTOCOLS"
    OWNER to postgres;

-- Table: IsaTab.STUDY_PUBLICATIONS

-- DROP TABLE "IsaTab"."STUDY_PUBLICATIONS";

CREATE TABLE "IsaTab"."STUDY_PUBLICATIONS"
(
    "ID" oid,
    "Study_PubMed_ID" text COLLATE pg_catalog."default",
    "Study_Publication_DOI" text COLLATE pg_catalog."default",
    "Study_Publication_Author_List" text COLLATE pg_catalog."default",
    "Study_Publication_Title" text COLLATE pg_catalog."default",
    "Study_Publication_Status" text COLLATE pg_catalog."default",
    "Study_Publication_Status_Term_Accession_Number" text COLLATE pg_catalog."default",
    "Study_Publication_Status_Term_Source_REF" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."STUDY_PUBLICATIONS"
    OWNER to postgres;

-- Table: IsaTab.S_SAMPLES

-- DROP TABLE "IsaTab"."S_SAMPLES";

CREATE TABLE "IsaTab"."S_SAMPLES"
(
    "ID" oid,
    "Sample_Name" text COLLATE pg_catalog."default",
    "Characteristics" text[] COLLATE pg_catalog."default",
    "Material_Type" text COLLATE pg_catalog."default",
    "Protocol_REF" text COLLATE pg_catalog."default",
    "Term_Accession_Number" text COLLATE pg_catalog."default",
    "Term_Source_REF" text COLLATE pg_catalog."default",
    "Comment" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."S_SAMPLES"
    OWNER to postgres;

-- Table: IsaTab.S_SOURCES

-- DROP TABLE "IsaTab"."S_SOURCES";

CREATE TABLE "IsaTab"."S_SOURCES"
(
    "ID" oid,
    "Source_Name" text COLLATE pg_catalog."default",
    "Characteristics" text COLLATE pg_catalog."default",
    "Material_Type" text COLLATE pg_catalog."default",
    "Term_Source_REF" text COLLATE pg_catalog."default",
    "Term_Accession_Number" text COLLATE pg_catalog."default",
    "Provider" text COLLATE pg_catalog."default",
    "Description" text COLLATE pg_catalog."default",
    "Comment" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."S_SOURCES"
    OWNER to postgres;

-- Table: IsaTab.ASSAY

-- DROP TABLE "IsaTab"."ASSAY";

CREATE TABLE "IsaTab"."ASSAY"
(
    "ID" oid,
    "Assay_Identifier" text COLLATE pg_catalog."default",
    "Assay_Measurement_Type" text COLLATE pg_catalog."default",
    "Assay_Technology_Type" text COLLATE pg_catalog."default",
    "Assay_Technology_Platform" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."ASSAY"
    OWNER to postgres;

-- Table: IsaTab.A_TECHNOLOGY_MICROARRAY

-- DROP TABLE "IsaTab"."A_TECHNOLOGY_MICROARRAY";

CREATE TABLE "IsaTab"."A_TECHNOLOGY_MICROARRAY"
(
    "ID" oid,
    "Sample_Name" text COLLATE pg_catalog."default",
    "Extract_Name" text COLLATE pg_catalog."default",
    "Labeled_Extract_Name" text COLLATE pg_catalog."default",
    "Label" text COLLATE pg_catalog."default",
    "Hybridization_Assay_Name" text COLLATE pg_catalog."default",
    "Comment[ArrayExpress_Accession]" text COLLATE pg_catalog."default",
    "Comment[ArrayExpress_Raw_Data_URL]" text COLLATE pg_catalog."default",
    "Comment[ArrayExpress_Processed_Data_URL]" text COLLATE pg_catalog."default",
    "Array_Design_REF" text COLLATE pg_catalog."default",
    "Scan_Name" text COLLATE pg_catalog."default",
    "Array_Data_File" text COLLATE pg_catalog."default",
    "Data_Transformation_Name" text COLLATE pg_catalog."default",
    "Derived_Array_Data_File" text COLLATE pg_catalog."default",
    "Array_Data_Matrix_File" text COLLATE pg_catalog."default",
    "Derived_Array_Data_Matrix_File" text COLLATE pg_catalog."default",
    "Array_Design_File" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."A_TECHNOLOGY_MICROARRAY"
    OWNER to postgres;

-- Table: IsaTab.A_TECHNOLOGY_MS

-- DROP TABLE "IsaTab"."A_TECHNOLOGY_MS";

CREATE TABLE "IsaTab"."A_TECHNOLOGY_MS"
(
    "ID" oid,
    "Sample_Name" text COLLATE pg_catalog."default",
    "Extract_Name" text COLLATE pg_catalog."default",
    "Labeled_Extract_Name" text COLLATE pg_catalog."default",
    "Label" text COLLATE pg_catalog."default",
    "MS_Assay_Name" text COLLATE pg_catalog."default",
    "Raw_Spectral_Data_File" text COLLATE pg_catalog."default",
    "Normalization_Name" text COLLATE pg_catalog."default",
    "Protein_Assignment_File" text COLLATE pg_catalog."default",
    "Peptide_Assignment_File" text COLLATE pg_catalog."default",
    "Post_Translational_Modification_Assignment_File" text COLLATE pg_catalog."default",
    "Data_Transformation_Name" text COLLATE pg_catalog."default",
    "Derived_Spectral_Data_File" text COLLATE pg_catalog."default",
    "Factor_Value[limiting nutrient]" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."A_TECHNOLOGY_MS"
    OWNER to postgres;

-- Table: IsaTab.A_TECHNOLOGY_NMR

-- DROP TABLE "IsaTab"."A_TECHNOLOGY_NMR";

CREATE TABLE "IsaTab"."A_TECHNOLOGY_NMR"
(
    "ID" oid,
    "Sample_Name" text COLLATE pg_catalog."default",
    "Extract_Name" text COLLATE pg_catalog."default",
    "Labeled_Extract_Name" text COLLATE pg_catalog."default",
    "Label" text COLLATE pg_catalog."default",
    "NMR_Assay_Name" text COLLATE pg_catalog."default",
    "Free_Induction_Decay_Data_File" text COLLATE pg_catalog."default",
    "Acquisition_Parameter_Data_File" text COLLATE pg_catalog."default",
    "Derived_Spectral_Data_File" text COLLATE pg_catalog."default",
    "Normalization_Name" text COLLATE pg_catalog."default",
    "Data_Transformation_Name" text COLLATE pg_catalog."default",
    "Metabolite_Assignment_File" text COLLATE pg_catalog."default"
)

TABLESPACE pg_default;

ALTER TABLE "IsaTab"."A_TECHNOLOGY_NMR"
    OWNER to postgres;
    