"""Convenience loaders for the infra_cptac assets. Usage:
    import sys; sys.path.insert(0, '/workspace/scripts/infra_cptac'); import cptac_load as cl
    ph  = cl.matrix('LUAD', 'phospho')                 # site_key x sample (float32); also 'protein','rna','cnv_log2','cnv_gistic','mut_binary', 'input_tc', ...
    sam = cl.samples()                                  # pan-cancer sample table
    z   = cl.activity('tc', 'corev3_z')                 # kinase x sample (all cohorts); variants tc|paired|tc_padj|paired_padj
    t   = cl.tumour_ids('LUAD')                         # tumour sample ids with phospho
"""
import pandas as pd

ROOT = '/workspace/results/infra_cptac'
COHORTS = ['BRCA', 'CCRCC', 'COAD', 'GBM', 'HNSCC', 'LSCC', 'LUAD', 'OV', 'PDAC', 'UCEC']


def matrix(cohort, name):
    return pd.read_parquet(f'{ROOT}/cohorts/{cohort}/{name}.parquet')


def table(cohort, name):            # clinical | annotations | mutations | samples | subtypes_cbioportal
    return pd.read_parquet(f'{ROOT}/cohorts/{cohort}/{name}.parquet')


def samples():
    return pd.read_parquet(f'{ROOT}/pancan_samples.parquet')


def clinical():
    return pd.read_parquet(f'{ROOT}/pancan_clinical.parquet')


def activity(variant='tc', what='corev3_z'):   # corev3_z | corev3_zcore | corev3_n | meanz_z | meanz_n
    return pd.read_parquet(f'{ROOT}/activity/{variant}/{what}.parquet')


def tumour_ids(cohort=None):
    s = samples()
    s = s[(s.sample_type == 'Tumor') & s.has_phospho]
    return list(s[s.cohort == cohort].sample_id) if cohort else list(s.sample_id)


def site_map():
    return pd.read_parquet(f'{ROOT}/mapping/site_map.parquet')
