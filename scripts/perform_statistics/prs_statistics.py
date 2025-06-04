import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

def load_data(metadata_path, pgs_path):
    meta = pd.read_csv(metadata_path, sep='\t')
    pgs = pd.read_csv(pgs_path, sep='\t')
    merged = pd.merge(meta, pgs, left_on = 'original_id', right_on='sample', how='inner')
    return merged

def run_ttest(df, case_col='case_control', pgs_col='PGS', case_label='case', control_label='control'):
    cases = df.loc[df[case_col] == case_label, pgs_col]
    controls = df.loc[df[case_col] == control_label, pgs_col]
    stat, pval = ttest_ind(cases, controls, equal_var=False)  # Welch's t-test
    return stat, pval

def run_logistic_regression(df, outcome_col='case_control', pgs_col='PGS',
                            covariates=['sex_label', 'PC1', 'PC2', 'PC3', 'PC4', 'age'],
                            case_label='case', control_label='control'):
    # Encode outcome as binary 0/1
    df = df.copy()
    df['y'] = df[outcome_col].apply(lambda x: 1 if x == case_label else 0)

    # Z-score covariates except sex_label
    covars_to_scale = [c for c in covariates if c != 'sex_label']
    scaler = StandardScaler()
    df[covars_to_scale] = scaler.fit_transform(df[covars_to_scale])

    # Encode sex_label as binary (assuming male/female)
    df['sex_label_bin'] = df['sex_label'].apply(lambda x: 1 if x.lower() == 'male' else 0)

    # Build design matrix with intercept
    X = df[[pgs_col] + covars_to_scale + ['sex_label_bin']]
    X = sm.add_constant(X)
    y = df['y']

    model = sm.Logit(y, X).fit(disp=0)

    # Collect results
    results = pd.DataFrame({
        'beta': model.params,
        'p_value': model.pvalues
    })
    return results

def write_results(ttest_stat, ttest_pval, logreg_results, outfile):
    with open(outfile, 'w') as f:
        f.write(f"Simple t-test comparing PGS between cases and controls\n")
        f.write(f"T-statistic: {ttest_stat:.4f}\n")
        f.write(f"P-value: {ttest_pval:.4e}\n\n")

        f.write("Logistic regression results: outcome = case/control\n")
        f.write("Predictors: PGS, sex_label (binary male=1), PC1-4, age (covariates z-scored except sex_label)\n\n")
        f.write(logreg_results.to_string(float_format='%.4e'))
        f.write('\n')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='PGS association tests')
    parser.add_argument('--metadata', required=True, help='TSV file with sample metadata (sample, case_control, sex_label, PC1-PC4, age)')
    parser.add_argument('--pgs', required=True, help='TSV file with sample and PGS columns')
    parser.add_argument('--out', required=True, help='Output text file for results')
    args = parser.parse_args()

    df = load_data(args.metadata, args.pgs)

    t_stat, t_pval = run_ttest(df)
    logreg_res = run_logistic_regression(df)

    write_results(t_stat, t_pval, logreg_res, args.out)
