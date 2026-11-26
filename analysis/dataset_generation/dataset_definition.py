from ehrql import create_dataset, get_parameter, days, years, minimum_of, maximum_of, case, when
from ehrql.tables.tpp import (
    patients,
    ons_deaths,
    practice_registrations,
)
from codelists import *
from variable_helper_functions import (
    first_record_before,
    service_usage_before,
    reside_eng_on,
    get_latest_ethnicity,
    first_matching_tpp_between,
    first_matching_apc_between,
    first_matching_vax_between,
    prevalent_tpp,
    prevalent_apc,
    prevalent_med,
    prevalent_vax,
    check_ic
)
from datetime import date

#get parameter
end_date = get_parameter(name="end_date")

#create dataset
dataset = create_dataset()

# Configure dummy data
dataset.configure_dummy_data(population_size=5000)

# dob and age
dataset.var_date_birth = patients.date_of_birth - days(patients.date_of_birth.day - 15)

dataset.cov_int_age = (date(2023,8,31) - dataset.var_date_birth).years

dataset.var_date_index = case(
    when(dataset.cov_int_age.is_in([65,66])).then(date(2023,9,1)),
    when(dataset.cov_int_age.is_in([63])).then(dataset.var_date_birth + years(65)),
    when(dataset.cov_int_age.is_in([64])).then(dataset.var_date_birth + years(65)),
    otherwise = None
    )

# dod
dataset.var_date_death = minimum_of(patients.date_of_death, ons_deaths.date)

# sex
dataset.cov_cat_sex = patients.sex.map_values(
    {"female": "Female", "male": "Male", "intersex": "Intersex"},
    default="Unknown"
)

# ethnicity
dataset.cov_cat_ethnicity = get_latest_ethnicity(
    dataset.var_date_index, ethnicity_codelist, grouping=6
)

# address
patient_address = reside_eng_on(dataset.var_date_index)
dataset.cov_cat_imd = patient_address.imd_quintile.map_values(
    {"1 (most deprived)": "1", "2": "2", "3": "3", "4": "4", "5 (least deprived)": "5"},
    default="Unknown"
)
dataset.cov_cat_urban = case(
    when(patient_address.rural_urban_classification.is_in([1,2,3,4])).then('urban'),
    when(patient_address.rural_urban_classification.is_in([5,6,7,8])).then('rural'),
    otherwise="Unknown"
)
dataset.var_bin_residency = patient_address.exists_for_patient()
dataset.var_date_residencyend = patient_address.end_date

# gp registration
patient_reg = practice_registrations.for_patient_on(dataset.var_date_index)
dataset.cov_cat_region = case(
    when(patient_reg.practice_nuts1_region_name.is_in(["North East", "North West", "Yorkshire and The Humber", "East Midlands", "West Midlands",
        "East", "London", "South East", "South West"])).then(patient_reg.practice_nuts1_region_name),
    otherwise = "Unknown"
    )
dataset.var_bin_registration = patient_reg.exists_for_patient()
dataset.var_date_deregistered = patient_reg.end_date
reg_start = patient_reg.start_date

# lookback date
dataset.var_date_lookback = first_record_before(dataset.var_date_index, dataset.var_date_death)

# follow up end date
fu_end_date = minimum_of(dataset.var_date_death, end_date, dataset.var_date_deregistered, dataset.var_date_residencyend)

# prevalent dementia
tmp_prev_tpp = prevalent_tpp(dementia_prev_snomed, dataset.var_date_index+days(31), dataset.var_date_death)
tmp_prev_apc = prevalent_apc(dementia_prev_icd, dataset.var_date_index+days(31), dataset.var_date_death)
tmp_prev_med = prevalent_med(dementia_med, dataset.var_date_index+days(31), dataset.var_date_death)
dataset.pre_bin_dementia = maximum_of(tmp_prev_tpp, tmp_prev_apc, tmp_prev_med)==1

# outcome: incident dementia
tmp_inci_tpp = first_matching_tpp_between(dementia_inci_snomed, fu_end_date, dataset.var_date_death, start_date=dataset.var_date_index+days(31))
tmp_inci_apc = first_matching_apc_between(dementia_inci_icd, fu_end_date, dataset.var_date_death, start_date=dataset.var_date_index+days(31))
dataset.out_date_dementia = minimum_of(tmp_inci_tpp, tmp_inci_apc)

# outcome: incident HZ
tmp_inci_tpp = first_matching_tpp_between(hz_snomed, fu_end_date, dataset.var_date_death, start_date=dataset.var_date_index+days(31))
tmp_inci_apc = first_matching_apc_between(hz_icd, fu_end_date, dataset.var_date_death, start_date=dataset.var_date_index+days(31))
dataset.out_date_hz = minimum_of(tmp_inci_tpp, tmp_inci_apc)

# immunocompromised
dataset.pre_bin_ic = check_ic(dataset.var_date_index, dataset.var_date_death)

# RZV ZVL vax before index date
tmp_prev_tpp = prevalent_tpp(shinglesvax_codelist, dataset.var_date_index, dataset.var_date_death)
tmp_prev_vax = prevalent_vax(shinglesvax_prev_names, dataset.var_date_index, dataset.var_date_death)
dataset.pre_bin_vax_shingles = maximum_of(tmp_prev_tpp, tmp_prev_vax)==1

#previous times of hospitalisation
dataset.pre_rate_hospitalisation = service_usage_before(dataset.var_date_index, maximum_of(dataset.var_date_index-years(2), reg_start), dataset.var_date_death, service="secondary")

# previous gp visits
dataset.pre_rate_gpvisit = service_usage_before(dataset.var_date_index, maximum_of(dataset.var_date_index-years(2),reg_start), dataset.var_date_death, service="primary")

# 1st rzv vaccination after index
dataset.var_date_1vax_rzv = first_matching_vax_between(shinglesvax_inci_names, fu_end_date, dataset.var_date_death, start_date=dataset.var_date_index)

# 2rd rzv vax
dataset.var_date_2vax_rzv = first_matching_vax_between(shinglesvax_inci_names, fu_end_date, dataset.var_date_death, start_date=dataset.var_date_1vax_rzv+days(1))

dataset.define_population(patients.exists_for_patient())


