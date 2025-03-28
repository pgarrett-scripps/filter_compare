import math
from statistics import mean
from io import StringIO

import pandas as pd
import streamlit as st
import serenipy.dtaselectfilter as filter
import plotly.express as px

st.set_page_config(
    page_title="FilterCompare",
    layout="wide"
)

st.title('DTASelect-filter compare! :bar_chart:')

files = st.file_uploader(label='DTASelect-filter files', type='.txt', accept_multiple_files=True)


def strip_modifications(peptide_sequence: str) -> str:
    return ''.join([c for c in peptide_sequence if c in 'ACDEFGHIKLMNPQRSTVWY'])


if st.button(label='Run :runner:'):

    data = []
    for file in files:
        name = file.name.split('.txt')[0]
        if len(file.name.split('_DTASelect-filter.txt')) > 1:
            name = file.name.split('_DTASelect-filter.txt')[0]
        elif len(file.name.split('DTASelect-filter.txt')) > 1:
            name = file.name.split('DTASelect-filter.txt')[0]
        _, _, results, end_lines = filter.from_dta_select_filter(StringIO(file.getvalue().decode('utf-8')))
        target_results = [result for result in results if
                          any(['Reverse_' not in protein_line.locus_name for protein_line in result.protein_lines])]
        target_protein_locuses = len(
            [protein_line.locus_name for result in target_results for protein_line in result.protein_lines])
        target_protein_groups = len([result.protein_lines[0].locus_name for result in target_results])
        target_stripped_peptides = len(
            {strip_modifications(peptide_line.sequence[2:-2]) for result in target_results for peptide_line in
             result.peptide_lines})
        target_charged_stripped_peptides = len(
            {(peptide_line.charge, strip_modifications(peptide_line.sequence[2:-2])) for result in target_results for
             peptide_line in result.peptide_lines})
        target_peptides = len(
            {(peptide_line.charge, peptide_line.sequence[2:-2]) for result in target_results for peptide_line in
             result.peptide_lines})
        target_spectra = sum([result.protein_lines[0].spectrum_count for result in target_results])
        target_protein_coverage = mean(
            [protein_line.sequence_coverage for result in target_results for protein_line in result.protein_lines])

        decoy_results = [result for result in results if
                         all(['Reverse_' in protein_line.locus_name for protein_line in result.protein_lines])]
        decoy_protein_locuses = len(
            [protein_line.locus_name for result in decoy_results for protein_line in result.protein_lines])
        decoy_protein_groups = len([result.protein_lines[0].locus_name for result in decoy_results])
        decoy_stripped_peptides = len(
            {strip_modifications(peptide_line.sequence[2:-2]) for result in decoy_results for peptide_line in
             result.peptide_lines})
        decoy_charged_stripped_peptides = len(
            {(peptide_line.charge, strip_modifications(peptide_line.sequence[2:-2])) for result in decoy_results for
             peptide_line in result.peptide_lines})
        decoy_peptides = len(
            {(peptide_line.charge, peptide_line.sequence[2:-2]) for result in decoy_results for peptide_line in
             result.peptide_lines})
        decoy_spectra = sum([result.protein_lines[0].spectrum_count for result in decoy_results])
        decoy_protein_coverage = mean(
            [protein_line.sequence_coverage for result in decoy_results for protein_line in result.protein_lines])

        protein_fdr = round(decoy_protein_groups/(decoy_protein_groups + target_protein_groups)*100,4)
        peptide_fdr = round(decoy_peptides / (decoy_peptides + target_peptides) * 100, 4)
        spectra_fdr = round(decoy_spectra / (decoy_spectra + target_spectra) * 100, 4)

        d = {'file': name,
             'target_protein_locuses': target_protein_locuses,
             'target_protein_groups': target_protein_groups,
             'target_stripped_peptides': target_stripped_peptides,
             'target_charged_stripped_peptides': target_charged_stripped_peptides,
             'target_peptides': target_peptides,
             'target_spectra': target_spectra,
             'target_protein_coverage': target_protein_coverage,

             'decoy_protein_locuses': decoy_protein_locuses,
             'decoy_protein_groups': decoy_protein_groups,
             'decoy_stripped_peptides': decoy_stripped_peptides,
             'decoy_charged_stripped_peptides': decoy_charged_stripped_peptides,
             'decoy_peptides': decoy_peptides,
             'decoy_spectra': decoy_spectra,
             'decoy_protein_coverage': decoy_protein_coverage,

             'protein_fdr': protein_fdr,
             'peptide_fdr': peptide_fdr,
             'spectra_fdr': spectra_fdr,

             'fdr_adj_protein_groups': math.floor(target_protein_groups - target_protein_groups*protein_fdr/100),
             'fdr_adj_peptides': math.floor(target_peptides - target_peptides * peptide_fdr / 100),
             'fdr_adj_spectra': math.floor(target_spectra - target_spectra * spectra_fdr / 100),

             }
        data.append(d)

    df = pd.DataFrame(data)
    st.dataframe(df)
    st.download_button(label='Download Table', data=df.to_csv(), file_name='filter_compare.csv')

    fig = px.bar(df, x="file", y=['target_protein_locuses', 'target_protein_groups', 'target_stripped_peptides',
                                  'target_charged_stripped_peptides', 'target_peptides', 'target_spectra'],
                 barmode="group",
                 text_auto=True, title='Target Protein & Peptide & Spectra Stats')
    st.plotly_chart(fig)

    fig = px.bar(df, x="file", y=['decoy_protein_locuses', 'decoy_protein_groups', 'decoy_stripped_peptides',
                                  'decoy_charged_stripped_peptides', 'decoy_peptides', 'decoy_spectra'],
                 barmode="group",
                 text_auto=True, title='Decoy Protein & Peptide & Spectra Stats')
    st.plotly_chart(fig)

    fig = px.bar(df, x="file", y=['target_protein_coverage', 'decoy_protein_coverage'], barmode="group", text_auto=True,
                 title='Target vs Decoy Protein Coverage')
    st.plotly_chart(fig)

    fig = px.bar(df, x="file", y=['protein_fdr', 'peptide_fdr', 'spectra_fdr'], barmode="group", text_auto=True,
                 title='Protein & Peptide & Spectra FDR')
    st.plotly_chart(fig)

    fig = px.bar(df, x="file", y=['fdr_adj_protein_groups', 'fdr_adj_peptides', 'fdr_adj_spectra'],
                 barmode="group",
                 text_auto=True, title='fdr adjusted Protein & Peptide & Spectra counts')
    st.plotly_chart(fig)


