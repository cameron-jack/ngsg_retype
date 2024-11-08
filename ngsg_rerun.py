#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@created: 31 Oct 2024
@author: Cameron Jack, ANU Bioinformatics Consultancy,
        JCSMR, Australian National University

Streamlit application to read a genotyping report and output a custom manifest
for rerunning failed assays.

        Run with: streamlit --run ngsgeno.py

"""
from csv import field_size_limit
from msilib.schema import File
from xmlrpc.client import boolean
import jsonpickle
import os
import sys
from pathlib import Path  
from subprocess import check_output, CalledProcessError, STDOUT
import subprocess
import pandas as pd
import chardet  # For automatic encoding detection
import streamlit as st
import extra_streamlit_components as stx
from time import sleep

###
# Use Eslam Ibrahim's code below to check for non-ASCII characters
###

######################################################################################################
# This script is used to validate the format of FASTA files.
# It is designed to be used in NGS genotyping pipelines.
# The script utilizes Streamlit to provide a web interface for the end-user.
# To run the script, use the command: #"python -m streamlit run fasta_checker.py"
# This script was developed and tested by Eslam Ibrahim (Eslam.ibrahim@anu.edu.au) on OCT. 10, 2024.
######################################################################################################

def check_non_ascii(file_content):
    """ A function to check for non-ASCII characters """
    issue_num = 0
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line_num, line in enumerate(lines):
        for idx, char in enumerate(line):
            if ord(char) > 127:  # ASCII characters have values from 0 to 127
                issue_num += 1
                issues.append({
                    'Issue Number': issue_num,
                    'Line Number': line_num,
                    'Issue': f"Invalid character in sequence header: '{char}' at position {idx + 1}"
                })
        line_num += 1  # Move to the next line number
    return issues  
#=====================================
# Function to check that sequences contain only A, T, C, G, N
# Spaces and tabs will not be reported here, only checked in check_gaps
def check_valid_sequence(file_content):
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line in lines:
        if line.startswith(">") or not line.strip():  # Ignore headers and blank lines
            line_num += 1
            continue
        
        sequence = line.strip().replace(" ", "").replace("\t", "")  # Clean the sequence

        for idx, char in enumerate(sequence):
            if char not in "AaTtCcGgNn":  # Check for valid characters
                global issue_num
                issue_num += 1
                issues.append({
                    'Issue Number': issue_num,
                    'Line Number': line_num,
                    'Issue': f"Invalid character in sequence: '{char}' at position {idx + 1}"
                })
        
        line_num += 1  # Move to the next line number

    return issues
#=====================================
# Function to check for gaps (spaces or tabs) within sequences
def check_gaps(file_content):
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line in lines:
        if line.startswith(">") or not line.strip():  # Ignore headers and blank lines
            line_num += 1
            continue

        if " " in line or "\t" in line:
            global issue_num
            issue_num += 1
            issues.append({
                'Issue Number': issue_num,
                'Line Number': line_num,
                'Issue': "Gap (space or tab) in sequence"
            })
        
        line_num += 1  # Move to the next line number

    return issues


def check_blank_lines(file_content):
    """
    FASTA files should not contain blank lines
    """
    issues = []
    lines = file_content.splitlines()
    line_num = 1  # Start with the first line

    for line in lines:
        if not line.strip():  # Blank line
            global issue_num
            issue_num += 1
            issues.append({
                'Issue Number': issue_num,
                'Line Number': line_num,
                'Issue': "Blank line found"
            })
        line_num += 1  # Move to the next line number
    return issues


def check_fasta_file(file_content):
    """
    Issue feature checks for FASTA specific file features
    """
    issues = []
    issues.extend(check_non_ascii(file_content))
    issues.extend(check_valid_sequence(file_content))
    issues.extend(check_gaps(file_content))
    issues.extend(check_blank_lines(file_content))
    return issues   


def parse_failed_genotyping_results(file_content):
    """
    Genotyping returns a tab-delimited report file with the results of the genotyping assessment for each assay.
    The header looks like this:
        barcode	code_assays	plate	wellLocation	sex	alleleSymbol	alleleKey	assayKey\
        passFail	seqName1	seqName2	seqName3	seqName4	seqName5	seqName6	\
        seqName7	efficiency	alleleRatio	alleleRatioAdjusted	genotype	args	reason
    Presumably we are working with passFail?
    
    We need code_assays (sample barcode plus assay), plate, wellLocation, passFail
    """
    failed_assays = {}  # (barcode,plate,wellLocation,sex) = [failed_assays]
    for line in file_content:
        if line.startswith('barcode'):
            continue  # header
        #print(line)
        cols = [c.strip() for c in line.split(',')]
        print(f'{cols=}')
        if len(cols) < 21:
            continue
        code_assays = cols[1]
        barcode, assay = code_assays.split(';')
        plate = cols[2]
        wellLocation = cols[3]
        sex = cols[4]
        alleleSymbol = cols[5].strip()
        alleleKey = cols[6].strip()
        assayKey = cols[7].strip()
        passFail = cols[8]
        genotype = cols[19]
        reason = cols[21]
        ident = (barcode, plate, wellLocation, sex)
        if '?' in genotype:
            if ident not in failed_assays:
                failed_assays[ident] = []
            failed_assays[ident].append((assay,alleleSymbol))
    return failed_assays


def generate_rerun_manifest(failed_assays, manifest_name='../Downloads/rerun_failures.csv'):
    """
    We need to output a CSV (comma delimited) with the following headers:
        Sample no	plateBarcode	well	sampleBarcode	Assay	Assay	Assay	Assay	\
        Assay	Assay	Assay	clientName	sampleName 	alleleSymbol
        
    Sample no. is automatically incremented from 1 for each sample.
    
    inputs: failed_assays {(barcode,plate,well,sex)=[(assay,alleleSymbol)]}
    """
    with open(manifest_name, 'wt') as f:
        # header
        print(','.join(['Sample no','plateBarcode','well','sampleBarcode','Assay','Assay',
                'Assay','Assay','Assay','Assay','Assay','clientName','sampleName','alleleSymbol']), file=f)
        
        for i, ident in enumerate(failed_assays):
            barcode, plate, wellLocation, sex = ident
            all_assays = []
            all_alleleSymbols = []
            for assay, alleleSymbol in failed_assays[ident]:
                all_assays.append(assay)
                all_alleleSymbols.append(alleleSymbol)
            row_array = [str(i+1),plate, wellLocation, barcode] + all_assays + ['']*(7-len(all_assays)) + ['rerun','',';'.join(all_alleleSymbols)]
            print(','.join(row_array), file=f)
    return manifest_name

#=====================================
# Streamlit UI

def main():
    """
    The NGSgeno "Xplorer" application. Allows full control of all sections of the pipeline,
    and displays all aspects of the experiment state at any time.
    """    
    st.set_page_config(
        page_title="NGSG: Rerun failed assays",
        page_icon="ngsg_icon.png",
        layout="wide"
    )
    st.title("NGS Genotyping: rerun failed assays")
    st.write("Upload a genotyping report file and a custom "+\
            "manifest will be generated to rerun failed assays")

    uploaded_file = st.file_uploader("Choose a file")
    if uploaded_file is not None:
        # Read the raw file
        rawfile = uploaded_file.read()
    
        # Detect encoding and Decode the file content
        result = chardet.detect(rawfile)
        charenc = result['encoding']
        file_content = rawfile.decode(charenc)
    
        # Check for issues in the file content
        issues = check_non_ascii(file_content)
    
        if issues:
            # Create DataFrame and sort by Line Number if necessary
            df = pd.DataFrame(issues)
            # df = df.sort_values(by="Line Number")  # Uncomment if you want to sort by line number
        
            # Count the number of issues found
            issue_num = len(issues)
        
            # Display the number of issues found and the table with custom size
            st.write(f"There are {issue_num} issues found in the file:")
            st.dataframe(df, height=400, width=1000)  # Increase table size
        else:
            file_content_lines = [line.strip() for line in file_content.split('\n')]
            print(f'{file_content_lines=}')
            st.write("The file is valid. Proceeding to build custom manifest for assay reruns")
            failed_assays = parse_failed_genotyping_results(file_content_lines)
            st.write('The following sample/assay combinations appear to have failed:')
            scrollable = st.container(height=400, border=True)
            for ident in failed_assays:
                barcode, plate, wellLocation, sex = ident
                assays = '; '.join([a.strip() for a,ak in failed_assays[ident] if a.strip() != ''])
                allele_keys = '; '.join([ak.strip() for a,ak in failed_assays[ident] if ak.strip() != ''])
                print(allele_keys)
                row = ' '.join(['Barcode:', barcode, '  plate:', plate, '  well:', wellLocation, '  assays:', assays, '  allele_keys:', allele_keys])
                scrollable.text(row)
            manifest_name = generate_rerun_manifest(failed_assays)
            st.write('')
            st.write(f'Rerun manifest written to {manifest_name}')

if __name__ == '__main__':
    main()