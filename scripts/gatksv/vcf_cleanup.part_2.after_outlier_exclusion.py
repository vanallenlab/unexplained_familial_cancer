#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 Ryan L. Collins, Noah Fields, and the Van Allen Laboratory at Dana-Farber Cancer Institute
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Second of two passes for final cleanup of UFC WGS SV callset
Based on a script from Gillani & Collins et al., Science, 2025
Part 2: after outlier sample exclusion
"""


new_infos = ['##INFO=<ID=OLD_ID,Number=1,Type=String,Description="Original ' + \
             'GATK-SV variant ID before polishing">',
             '##INFO=<ID=FAILED_COHORT_COMPARISONS,Number=.,Type=String,Description=' + \
             '"Pairs of cohorts with significantly different frequencies">']
new_filts = ['##FILTER=<ID=MANUAL_FAIL,Description="This variant failed ' +
             'post hoc manual review and should not be trusted.">']
infos_to_rm = ['NCR_TMP']



import argparse
import csv
import numpy as np
import pandas as pd
import pybedtools as pbt
import pysam
from scipy.stats import fisher_exact
from sys import stdin, stdout


def is_multiallelic(record):
    """
    Reports whether a pysam.VariantRecord is an mCNV
    """

    return len(record.alleles) > 2 \
           or 'MULTIALLELIC' in record.filter.keys() \
           or record.info.get('SVTYPE') in 'CNV MCNV'.split()


def is_depth_only(record):
    """
    Checks whether a record is a read depth-only CNV
    """

    if record.info['ALGORITHMS'] == ('depth', ) \
    and 'PE' not in record.info['EVIDENCE'] \
    and 'SR' not in record.info['EVIDENCE']:
        return True
    else:
        return False


def update_af(record):
    """
    Update AC, AN, and AF for a single record
    """

    ac, an = 0, 0
    af = None
    for sdat in record.samples.values():
        GT = [a for a in sdat['GT'] if a is not None]
        an += len(GT)
        ac += len([a for a in GT if a > 0])
    if an > 0:
        af = ac / an
    record.info['AC'] = ac
    record.info['AN'] = an
    record.info['AF'] = af

    return record


def is_artifact_deletion(record):
    """
    Checks whether a record is a deletion in the artifact zone
    """

    svtype = record.info['SVTYPE']
    svlen = record.info['SVLEN']
    if svtype == "DEL" \
    and svlen > 400 \
    and svlen < 1000 \
    and record.info.get('AC', (0, ))[0] / record.info.get('AN', 1) < 0.05:
        return True
    else:
        return False


def calc_ncr(record, exclude_samples=[]):
    """
    Calculate no-call rate for a record's genotypes
    """

    total = 0
    nocalls = 0
    for sid, sinfo in record.samples.items():
        
        if sid in exclude_samples:
            continue
        
        total += 1

        if all([a is None for a in sinfo.get('GT', (None, ))]):
            nocalls += 1

    if total > 0:
        return nocalls / total
    else:
        return None


def recalibrate_qual(record):
    """
    Recalibrate record QUAL (quality)
    Defined as median GQ among all non-reference GTs
    """

    gqs = []
    for sinfo in record.samples.values():
        a = [a for a in sinfo.get('GT', (None, None)) if a is not None]
        if np.nansum(a) > 0:
            gqs.append(sinfo.get('GQ', None))

    if len(gqs) > 0:
        return np.nanmedian(gqs)
    else:
        return 0


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf')    
    parser.add_argument('vcf_out', help='output .vcf')
    parser.add_argument('-p', '--ped-file', help='PLINK-style .ped file. Used ' +
                        'for handling no-call rates on allosomes. Required.',
                        required=True, type=str)
    parser.add_argument('--fail-variants', help='list of variant IDs to be ' +
                        'marked as having failed manual review')
    parser.add_argument('--version-number', help='callset version number for ' +
                        'tagging variant IDs.', type=str)
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)
    samples = [s for s in invcf.header.samples]

    mf_infos = [k for k in invcf.header.info.keys() \
                if k.startswith('MALE_') or k.startswith('FEMALE_')]

    # Reformat header
    header = invcf.header
    for key in mf_infos:
        header.info.remove_header(key)
    for info in new_infos:
        header.add_line(info)
    for filt in new_filts:
        header.add_line(filt)

    # Make list of male/female samples for handling sex chromosome NCRs
    male_ids = set()
    female_ids = set()
    with open(args.ped_file) as ped_in:
        for fid, sid, dad, mom, sex, pheno in csv.reader(ped_in, delimiter='\t'):
            if sid not in samples:
                continue
            if int(sex) == 1:
                male_ids.add(sid)
            if int(sex) == 2:
                female_ids.add(sid)

    # Load list of variant IDs to manually fail, if optioned
    if args.fail_variants is not None:
        with open(args.fail_variants) as fin:
            fail_vids = [l.rstrip() for l in fin.readlines()]
    else:
        fail_vids = []

    # Open connection to output vcf
    if args.vcf_out in '- stdout /dev/stdout':
        outvcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=header)

    # Iterate over records in invcf, clean up, and write to outvcf
    svtype_counter = {}
    for record in invcf.fetch():

        # Get reused record info
        svtype = record.info.get('SVTYPE')
        svlen = record.info.get('SVLEN', 0)

        # Clear unnecessary INFOs
        for key in infos_to_rm:
            if key in record.info.keys():
                record.info.pop(key)

        # Check if record should be marked as manual fail
        if record.id in fail_vids:
            record.filter.add('MANUAL_FAIL')

        # Clear all MALE/FEMALE AF annotations
        for key in mf_infos:
            if key in record.info.keys():
                record.info.pop(key)

        # Update AC/AN/AF
        if not is_multiallelic(record):
            record = update_af(record)

        # Skip empty records
        if record.info.get('AC', (1, ))[0] == 0 \
        and not is_multiallelic(record):
            continue

        # Recompute NCR
        if not is_multiallelic(record):
            if record.chrom == 'chrX':
                record.info['NCR'] = calc_ncr(record, exclude_samples=male_ids)
            elif record.chrom == 'chrY':
                record.info['NCR'] = calc_ncr(record, exclude_samples=female_ids)
            else:
                record.info['NCR'] = calc_ncr(record)

        # Clear old high NCR FILTER and reannotate based on updated NCR
        original_filters = [k for k in record.filter.keys()]
        record.filter.clear()
        for k in original_filters:
            if k not in 'HIGH_NCR HIGH_PCRMINUS_NOCALL_RATE'.split():
                record.filter.add(k)
        if is_artifact_deletion(record):
            if record.info.get('NCR', 0) > 1/250:
                record.filter.add('HIGH_NCR')
        elif record.info.get('NCR', 0) >= 0.04:
            record.filter.add('HIGH_NCR')

        # Recalibrate QUAL score
        if is_multiallelic(record):
            record.qual = 99
        else:
            record.qual = recalibrate_qual(record)

        # Rename record
        if record.chrom not in svtype_counter.keys():
            svtype_counter[record.chrom] = {}
        if svtype not in svtype_counter[record.chrom].keys():
            svtype_counter[record.chrom][svtype] = 0
        svtype_counter[record.chrom][svtype] += 1
        if args.version_number is None:
            id_prefix = 'dfci-ufc'
        else:
            id_prefix = 'dfci-ufc.v' + str(args.version_number)
        new_id = '_'.join([id_prefix, svtype, record.chrom,
                           str(svtype_counter[record.chrom][svtype])])
        record.info['OLD_ID'] = record.id
        record.id = new_id

        # Write to outvcf
        outvcf.write(record)

    # Close connection to output file to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

