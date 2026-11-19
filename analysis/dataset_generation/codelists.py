from ehrql import codelist_from_csv

ethnicity_codelist = codelist_from_csv(
    "codelists/opensafely-ethnicity-snomed-0removed.csv",
    column="code",
    category_column="Label_6",
)

shinglesvax_codelist = codelist_from_csv(
    "codelists/user-sw-gsk_shingles_vax_snomed.csv",
    column="code"
)

shinglesvax_prev_names = ['Shingles (Herpes Zoster) adjuvanted rcmb vacc powder and suspension for suspension inj 0.5ml vials',
                     'Shingles (Herpes Zoster) vaccine (live) powder and solvent for suspension for injection 0.65ml pfs',
                     'Shingrix',
                     'Zostavax']

shinglesvax_inci_names = ['Shingles (Herpes Zoster) adjuvanted rcmb vacc powder and suspension for suspension inj 0.5ml vials',
                     'Shingrix'
                     ]

#####Comorbidity
# db_snomed = codelist_from_csv(
#     "codelists/bristol-multimorbidity_diabetes.csv",
#     column="code"
# )

# db1_icd = codelist_from_csv(
#     "codelists/reducehf-type-1-diabetes-icd10.csv",
#     column="code"
# )

# db2_icd = codelist_from_csv(
#     "codelists/reducehf-type-2-diabetes-icd10.csv",
#     column="code"
# )

# ihd_snomed = codelist_from_csv(
#     "codelists/opensafely-ischaemic-heart-disease.csv",
#     column="code"
# )

# mi_snomed = codelist_from_csv(
#     "codelists/nhsd-primary-care-domain-refsets-mi_cod.csv",
#     column="code"
# )

# hf_snomed = codelist_from_csv(
#     "codelists/bristol-multimorbidity_heart-failure.csv",
#     column="code"
# )

# ihd_icd = codelist_from_csv(
#     "codelists/user-sw-gsk-comorbidity-ischaemic-heart-disease-icd10.csv",
#     column="code"
# )

# mi_icd = codelist_from_csv(
#     "codelists/reducehf-myocardial-infarction-icd10.csv",
#     column="code"
# )

# hf_icd = codelist_from_csv(
#     "codelists/reducehf-heart-failure-secondary-care.csv",
#     column="code"
# )

# stroke_snomed = codelist_from_csv(
#     "codelists/bristol-multimorbidity_stroketransient-ischemic-attack.csv",
#     column="code"
# )

# stroke_icd = codelist_from_csv(
#     "codelists/opensafely-stroke-secondary-care.csv",
#     column="icd"
# )

# depress_snomed = codelist_from_csv(
#     "codelists/nhsd-primary-care-domain-refsets-depr_cod.csv",
#     column="code"
# )

# depress_icd = codelist_from_csv(
#     "codelists/user-sw-gsk-comorbidity-depression-icd10.csv",
#     column="code"
# )

# parkinson_snomed = codelist_from_csv(
#     "codelists/bristol-burden-of-neurodegenerative-diseases-parkinsons-disease-snomed.csv",
#     column="code"
# )

# parkinson_icd = codelist_from_csv(
#     "codelists/bristol-burden-of-neurodegenerative-diseases-parkinsons-disease.csv",
#     column="code"
# )

# huntington_snomed = codelist_from_csv(
#     "codelists/bristol-burden-of-neurodegenerative-diseases-dementia-due-to-huntingtons-snomed.csv",
#     column="code"
# )

# huntington_icd = codelist_from_csv(
#     "codelists/bristol-burden-of-neurodegenerative-diseases-huntington-disease.csv",
#     column="code"
# )

# ms_snomed = codelist_from_csv(
#     "codelists/bristol-cambridge-multimorbidity-score-multiple-sclerosis.csv",
#     column="code"
# )

# ms_icd = codelist_from_csv(
#     "codelists/nhsd-multiple-sclerosis.csv",
#     column="code"
# )

# mci_snomed = codelist_from_csv(
#     "codelists/user-sw-gsk-comorbidity-mild-cognitive-depression-snomed.csv",
#     column="code"
# )

# mci_icd = codelist_from_csv(
#     "codelists/user-sw-gsk-comorbidity-mild-cognitive-depression-icd10.csv",
#     column="code"
# )

# hyperchol_snomed = codelist_from_csv(
#     "codelists/user-sw-gsk-comorbidity-hypercholesterolaemia-snomed.csv",
#     column="code"
# )

# hypertension_snomed = codelist_from_csv(
#     "codelists/nhsd-primary-care-domain-refsets-hyp_cod.csv",
#     column="code"
# )

# hypertension_icd = codelist_from_csv(
#     "codelists/reducehf-hypertension-icd10.csv",
#     column="code"
# )

#####IC status
lymphoma_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-leukaemia-lymphoma.csv",
    column="code"
)

lymphprolif_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-lymphoproliferative-disorders.csv",
    column="code"
)

hiv_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-hiv-aids.csv",
    column="code"
)

lymphocyte_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-functional-lymphocyte-disorder.csv",
    column="code"
)

cell_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-stem-cell-therapy.csv",
    column="code"
)

gvhd_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-gvhd.csv",
    column="code"
)

chemo_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-chemotherapy.csv",
    column="code"
)

radio_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-radiotherapy.csv",
    column="code"
)

organ_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-solid-organ-transplant.csv",
    column="code"
)

immune_med = codelist_from_csv(
    "codelists/user-sw-gsk-immunocompromised-medication.csv",
    column="code"
)


#####Outcome
dementia_prev_icd = codelist_from_csv(
    "codelists/user-sw-gsk-prev-dementia-icd.csv",
    column="code"
)

dementia_prev_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-prev-dementia-snomed.csv",
    column="code"
)

dementia_inci_icd = codelist_from_csv(
    "codelists/user-sw-gsk-outcome-dementia-icd.csv",
    column="code"
)

dementia_inci_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-outcome-dementia-snomed.csv",
    column="code"
)

dementia_med = codelist_from_csv(
    "codelists/user-sw-gsk-prev-dementia-dmd.csv",
    column="code"
)

hz_snomed = codelist_from_csv(
    "codelists/user-sw-gsk-outcome-herpes-zoster-infection-snomed.csv",
    column="code"
)

hz_icd = codelist_from_csv(
    "codelists/user-sw-gsk-outcome-herpes-zoster-infection-icd.csv",
    column="code"
)


