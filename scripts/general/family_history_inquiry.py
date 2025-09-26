import pandas
import os
from datetime import date
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# This query represents dataset "Cancer_WGS_Dataset" for domain "survey" and was generated for All of Us Controlled Tier Dataset v7
dataset_18935163_survey_sql = """
    SELECT
        answer.person_id,
        answer.survey,
        answer.question_concept_id,
        answer.question,
        answer.answer_concept_id,
        answer.answer  
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.ds_survey` answer   
    WHERE
        (
            question_concept_id IN (SELECT
                DISTINCT concept_id                         
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c                         
            JOIN
                (SELECT
                    CAST(cr.id as string) AS id                               
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr                               
                WHERE
                    concept_id IN (1740639)                               
                    AND domain_id = 'SURVEY') a 
                    ON (c.path like CONCAT('%', a.id, '.%'))                         
            WHERE
                domain_id = 'SURVEY'                         
                AND type = 'PPI'                         
                AND subtype = 'QUESTION')
        )  
        AND (
            answer.PERSON_ID IN (SELECT
                distinct person_id  
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
            WHERE
                cb_search_person.person_id IN (SELECT
                    person_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
                WHERE
                    has_whole_genome_variant = 1 ) 
                AND cb_search_person.person_id IN (SELECT
                    criteria.person_id 
                FROM
                    (SELECT
                        DISTINCT person_id, entry_date, concept_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                    WHERE
                        (concept_id IN(SELECT
                            DISTINCT c.concept_id 
                        FROM
                            `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                        JOIN
                            (SELECT
                                CAST(cr.id as string) AS id       
                            FROM
                                `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr       
                            WHERE
                                concept_id IN (443392)       
                                AND full_text LIKE '%_rank1]%'      ) a 
                                ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                                OR c.path LIKE CONCAT('%.', a.id) 
                                OR c.path LIKE CONCAT(a.id, '.%') 
                                OR c.path = a.id) 
                        WHERE
                            is_standard = 1 
                            AND is_selectable = 1) 
                        AND is_standard = 1 )) criteria ) 
                AND cb_search_person.person_id IN (SELECT
                    criteria.person_id 
                FROM
                    (SELECT
                        DISTINCT person_id, entry_date, concept_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                    WHERE
                        (concept_id IN (836768) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528469) 
                        OR  concept_id IN (836768) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528467) 
                        OR  concept_id IN (836768) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528470) 
                        OR  concept_id IN (836768) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530276) 
                        OR  concept_id IN (836768) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528471) 
                        OR  concept_id IN (836769) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530278) 
                        OR  concept_id IN (836769) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530273) 
                        OR  concept_id IN (836769) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528477) 
                        OR  concept_id IN (836769) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530277) 
                        OR  concept_id IN (836769) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528478) 
                        OR  concept_id IN (836770) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528484) 
                        OR  concept_id IN (836770) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528485) 
                        OR  concept_id IN (836770) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528482) 
                        OR  concept_id IN (836770) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530279) 
                        OR  concept_id IN (836770) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528483) 
                        OR  concept_id IN (836771) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528496) 
                        OR  concept_id IN (836771) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528494) 
                        OR  concept_id IN (836771) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528497) 
                        OR  concept_id IN (836771) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528493) 
                        OR  concept_id IN (836771) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528498) 
                        OR  concept_id IN (836772) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528504) 
                        OR  concept_id IN (836772) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528502) 
                        OR  concept_id IN (836772) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528505) 
                        OR  concept_id IN (836772) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528501) 
                        OR  concept_id IN (836772) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528506) 
                        OR  concept_id IN (836773) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528543) 
                        OR  concept_id IN (836773) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528544) 
                        OR  concept_id IN (836773) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528540) 
                        OR  concept_id IN (836831) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528569) 
                        OR  concept_id IN (836831) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528567) 
                        OR  concept_id IN (836831) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528570) 
                        OR  concept_id IN (836831) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528566) 
                        OR  concept_id IN (836831) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528571) 
                        OR  concept_id IN (836832) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528672) 
                        OR  concept_id IN (836832) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528671) 
                        OR  concept_id IN (836832) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530420) 
                        OR  concept_id IN (836832) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530419) 
                        OR  concept_id IN (836832) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528670) 
                        OR  concept_id IN (836833) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530423) 
                        OR  concept_id IN (836833) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530424) 
                        OR  concept_id IN (836833) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530421) 
                        OR  concept_id IN (836774) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530428) 
                        OR  concept_id IN (836774) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530427) 
                        OR  concept_id IN (836774) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530429) 
                        OR  concept_id IN (836774) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530426) 
                        OR  concept_id IN (836774) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530430) 
                        OR  concept_id IN (836775) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530433) 
                        OR  concept_id IN (836775) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530432) 
                        OR  concept_id IN (836775) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530434) 
                        OR  concept_id IN (836775) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530431) 
                        OR  concept_id IN (836775) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528697) 
                        OR  concept_id IN (836834) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528891) 
                        OR  concept_id IN (836834) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528889) 
                        OR  concept_id IN (836834) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528892) 
                        OR  concept_id IN (836834) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528888) 
                        OR  concept_id IN (836834) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43528893) 
                        OR  concept_id IN (836776) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529140) 
                        OR  concept_id IN (836776) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529138) 
                        OR  concept_id IN (836776) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529141) 
                        OR  concept_id IN (836776) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529137) 
                        OR  concept_id IN (836776) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529142) 
                        OR  concept_id IN (836777) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529188) 
                        OR  concept_id IN (836777) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529186) 
                        OR  concept_id IN (836777) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529189) 
                        OR  concept_id IN (836777) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529185) 
                        OR  concept_id IN (836777) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529190) 
                        OR  concept_id IN (836778) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529679) 
                        OR  concept_id IN (836778) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529680) 
                        OR  concept_id IN (836778) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529677) 
                        OR  concept_id IN (836779) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530572) 
                        OR  concept_id IN (836779) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530571) 
                        OR  concept_id IN (836779) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530573) 
                        OR  concept_id IN (836779) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530570) 
                        OR  concept_id IN (836779) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43530574) 
                        OR  concept_id IN (836780) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529735) 
                        OR  concept_id IN (836780) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529738) 
                        OR  concept_id IN (836780) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529739) 
                        OR  concept_id IN (836781) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529821) 
                        OR  concept_id IN (836781) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529819) 
                        OR  concept_id IN (836781) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529822) 
                        OR  concept_id IN (836781) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529818) 
                        OR  concept_id IN (836781) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529823) 
                        OR  concept_id IN (836782) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529871) 
                        OR  concept_id IN (836782) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529869) 
                        OR  concept_id IN (836782) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529872) 
                        OR  concept_id IN (836782) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529868) 
                        OR  concept_id IN (836782) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529873) 
                        OR  concept_id IN (836783) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529914) 
                        OR  concept_id IN (836783) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529912) 
                        OR  concept_id IN (836783) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529915) 
                        OR  concept_id IN (836783) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529911) 
                        OR  concept_id IN (836783) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529916) 
                        OR  concept_id IN (836835) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529629) 
                        OR  concept_id IN (836835) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529627) 
                        OR  concept_id IN (836835) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529630) 
                        OR  concept_id IN (836835) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529626) 
                        OR  concept_id IN (836835) 
                        AND is_standard = 0  
                        AND  value_source_concept_id IN (43529631))) criteria ) )
        )"""

dataset_18935163_survey_df = pandas.read_gbq(
    dataset_18935163_survey_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")

substring = "Including yourself, who in your family has had"

# Filter rows containing both the substring and 'cancer'
df_filtered = df[
    df['question'].str.contains(substring, na=False) &
    df['question'].str.contains("cancer", na=False, case=False) &
    ~df['question'].str.contains("other cancer(s)", na=False, case=False) &
    ~df['question'].str.contains("non-cancer", na=False, case=False)
]

# Remove rows where answer is 'PMI: Skip'
df_filtered_clean = df_filtered[df_filtered['answer'] != 'PMI: Skip'].copy()

df_filtered_clean = df_filtered_clean[~df_filtered_clean['answer'].str.contains("grandparent|grandchild", case=False, na=False)]
df_filtered_clean.to_csv("family_cancer_filtered.csv", index=False)




# Load family history filtered file
family = pd.read_csv("family_cancer_filtered.csv")

# Remove self-reported cancers
family = family[~family['answer'].str.contains("Self", na=False)]

# Map family history questions to cancer type
question_to_cancer = {
    'bladder cancer': 'Bladder',
    'blood or soft tissue cancer': 'Blood_Soft_Tissue',
    'bone cancer': 'Bone',
    'brain cancer': 'Brain',
    'breast cancer': 'Breast',
    'cervical cancer': 'Cervix',
    'esophageal cancer': 'Esophagus',
    'eye cancer': 'Eye',
    'kidney cancer': 'Kidney',
    'lung cancer': 'Lung',
    'ovarian cancer': 'Ovary',
    'pancreatic cancer': 'Pancreas',
    'prostate cancer': 'Prostate',
    'skin cancer': 'Skin',
    'stomach cancer': 'Stomach',
    'thyroid cancer': 'Thyroid',
    'colon cancer/rectal cancer': 'Colorectal',
    'endocrine cancer': 'Endocrine',
    'endometrial cancer': 'Uterus',
    'head and neck cancer': 'Head_Neck'
}

family['cancer_type'] = family['question'].map(
    lambda q: next((v for k, v in question_to_cancer.items() if k in q.lower()), None)
)

# Drop rows where mapping failed
family = family.dropna(subset=['cancer_type'])

# Collapse by sample_id → collect unique family cancers
family_collapsed = (
    family.groupby('person_id')['cancer_type']
    .apply(lambda x: ";".join(sorted(set(x))))
    .reset_index()
    .rename(columns={'person_id': 'sample_id', 'cancer_type': 'cancers_in_FDRs'})
)

# Write with a header block
out_file = "dfci-ufc.aou.family_cancers.tsv"
with open(out_file, "w") as f:
    f.write(f"Created: # {date.today().isoformat()}\n")
    f.write("# File listing cancers reported in first-degree relatives (FDRs) for each sample_id\n")
    f.write("# This file contains everyone in the AoU cohort, and needs to be trimmed down. \n")
    family_collapsed.to_csv(f, sep="\t", index=False)


# --- filter cohort cancers first ---
cohort_counts = cohort_matrix.sum()
cohort_keep = cohort_counts[cohort_counts >= 5].index  # only cancers with ≥5 cases
cohort_keep = cohort_keep.difference(["Environmental_miscellaneous","Control","Urinary_system","Gynecologic_system", "Parathyroid","Basal_cell_carcinoma","Squamous_cell_carcinoma","Gastrointestinal","Skin","Other"])

cohort_matrix_filt = cohort_matrix[cohort_keep]

# --- build co-occurrence counts ---
cooccurrence_counts = pd.DataFrame(
    0, index=family_matrix.columns, columns=cohort_matrix_filt.columns
)

for fam_cancer in family_matrix.columns:
    person_has_fam_cancer = (family_matrix[fam_cancer] > 0).reindex(
        cohort_matrix_filt.index, fill_value=False
    )
    for cohort_cancer in cohort_matrix_filt.columns:
        # count number of people with both
        count = cohort_matrix_filt.loc[person_has_fam_cancer, cohort_cancer].sum()
        cooccurrence_counts.loc[fam_cancer, cohort_cancer] = count

# --- save to file ---
cooccurrence_counts.to_csv("family_vs_cohort_cancer_counts.tsv", sep="\t")

# --- heatmap visualization ---
plt.figure(figsize=(16, 10))
sns.heatmap(
    cooccurrence_counts,
    cmap="Reds",
    annot=True,
    fmt="d",
    cbar_kws={'label': 'Number of individuals'}
)
plt.title("Family Cancer vs Patient Cancer (Counts)")
plt.xlabel("Patient (Cohort) Cancer")
plt.ylabel("Family Cancer History")
plt.tight_layout()
plt.savefig("Family_Patient_Count.png")
plt.show()
