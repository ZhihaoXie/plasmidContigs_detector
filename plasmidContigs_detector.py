#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# FileName:  plasmidContigs_detector.py
# Author:    Zhihao Xie  \(#^o^)/
# Date:      2017/12/8 11:13
# Version:   v1.0.0
# CopyRight: Copyright ©Zhihao Xie, All rights reserved.

import re, sys, os
import argparse
import shutil
from Bio import SeqIO

script_path = os.path.split(os.path.realpath(sys.argv[0]))[0]

def get_params():
    usage = "%(prog)s [options]"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('-f', dest='contigs', help="contigs fasta")
    parser.add_argument('-o', dest='outdir', default="plasmid_detectOut", help="the output dir")
    parser.add_argument('-x', dest='prefix', help="the output file prefix")
    parser.add_argument('-p', dest='pfam', default=os.path.join(script_path,"plasmid_pfam.id"), help="the pfam of plasmid copy family")
    parser.add_argument('-c', dest="cbar", default="/sdd/Biosoft/cBar.1.2/cBar.pl", help="the cBar path, you can get it from http://csbl.bmb.uga.edu/~ffzhou/cBar/")
    parser.add_argument('-a', dest="abricate", default="/sdd/Biosoft/abricate/bin/abricate", help="the abricate path, you can get it from https://github.com/tseemann/abricate")
    parser.add_argument('-i', dest="ipr", default="/sdd/userLogin/xiezh/biosoft/interproscan-5.16-55.0/interproscan.sh",
                        help="Interproscan path, get it from https://www.ebi.ac.uk/interpro/")
    options = parser.parse_args()
    if not options.contigs or not options.prefix:
        parser.print_help()
        sys.exit(1)

    return options

def intersection(file1, file2):
    """
    Args:     
    Returns:
    """
    set1 = set()
    with open(file1) as fh:
        for line in fh:
            if re.search(r"^#|^\s*$", line):
                continue
            if re.search('Plasmid', line):
                array = line.strip().split("\t")
                set1.add(array[0])
    set2 = set()
    with open(file2) as fh:
        for line in fh:
            if re.search(r"^#|^\s*$", line):
                continue
            array = line.strip().split("\t")
            set2.add(array[1])
    if len(set1) > 0 and len(set2) > 0:
        interList = list(set1 & set2)
    else:
        if len(set1) > 0 and len(set2) == 0:
            interList = list(set1)
        elif len(set2) > 0 and len(set1) == 0:
            interList = list(set2)
        else:
            interList = []
    return interList

def choosePlasmid(pfam, iprout, gff, contigs, outfasta1, outfasta2):

    pfam_hash = {}
    with open(pfam) as fh:
        for line in fh:
            if re.search(r"^#|^\s*$|^Pfam_ID", line):
                continue
            array = line.strip().split("\t")
            pfam_hash[array[0]] = 1

    gff_hash = {}
    with open(gff) as fh:
        for line in fh:
            if re.search("CDS", line):
                fields = line.strip().split("\t")
                cds_id = re.findall(r"locus_tag=([^;]+);?", fields[-1])[0]
                gff_hash[cds_id] = fields[0]

    plasmid_set = set()
    with open(iprout) as fh:
        for line in fh:
            if re.search(r"^#|^\s*$", line):
                continue
            array = line.strip().split("\t")
            if array[2] in pfam_hash:
                plasmid_set.add(gff_hash[array[0]])

    outfastafh1 = open(outfasta1, 'w')
    outfastafh2 = open(outfasta2, 'w')
    for rec in SeqIO.parse(contigs, 'fasta'):
        if rec.id in plasmid_set:
            SeqIO.write(rec, outfastafh1, 'fasta')
        else:
            SeqIO.write(rec, outfastafh2, 'fasta')
    outfastafh1.close()
    outfastafh2.close()

def chack_exist(file_path):
    if not os.path.exists(file_path):
        sys.stderr.write("Error: the %s don't exist.\n" % (file_path))
        sys.exit(1)

def main():

    options = get_params()
    working_dir = os.path.abspath(options.outdir)
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    else:
        sys.stderr.write("Warning: the %s had exist.\n" % (working_dir))
        sys.stderr.flush()
    out_prefix = options.prefix
    contigs = os.path.abspath(options.contigs)
    plasmid_pfamF = os.path.abspath(options.pfam)
    cBar_pl = os.path.abspath(options.cbar)
    abricate = os.path.abspath(options.abricate)
    ipr_sh = os.path.abspath(options.ipr)
    chack_exist(contigs)
    chack_exist(plasmid_pfamF)
    chack_exist(cBar_pl)
    chack_exist(abricate)
    chack_exist(ipr_sh)

    # run cBar
    cbar_output = os.path.join(working_dir, out_prefix+".cBar.txt")
    cmd_sh = """{cbar} {fasta} {out}""".format(cbar=cBar_pl, fasta=contigs, out=cbar_output)
    os.system(cmd_sh)
    # run plasmidFinder
    plasmidfinder_out = os.path.join(working_dir, out_prefix+".plasmidfinder.txt")
    cmd_sh = """{ab} --db plasmidfinder {fasta} > {out} 2>/dev/null""".format(ab=abricate, fasta=contigs, out=plasmidfinder_out)
    os.system(cmd_sh)
    if not os.path.exists(cbar_output) and not os.path.exists(plasmidfinder_out):
        sys.stderr.write("Error: the results about cBar, plasmidFinder don't exist.\n")
        sys.exit(1)

    interList = intersection(cbar_output, plasmidfinder_out)
    if len(interList) == 0:
        sys.stderr.write("Info: %s don't have plasmid.\n" % (os.path.basename(contigs)))
        sys.stderr.flush()

    # extract contigs
    tmpcontigs = os.path.join(working_dir, out_prefix+".tmpContigs.fasta")
    contigs_out = open(tmpcontigs, 'w')
    for rec in SeqIO.parse(contigs, 'fasta'):
        if rec.id in interList:
            SeqIO.write(rec, contigs_out, 'fasta')
    contigs_out.close()

    # predict
    tmp_dir = os.path.join(working_dir, out_prefix+"_predict")
    os.system("python3 %s --fasta %s -d %s --name %s --tr_table 11" % (script_path+os.sep+"src"+os.sep+"gene_predict4pro.py", tmpcontigs, tmp_dir, out_prefix))
    protein_faa = tmp_dir+os.sep+out_prefix+".protein.faa"
    protein_gff = tmp_dir+os.sep+out_prefix+".gene.gff"
    if not os.path.exists(protein_faa):
        sys.stderr.write("Error: the protein of %s don't exists.\n" % protein_faa)
        sys.stderr.flush()

    # run ipr
    cmd_str = r"""#!/bin/bash
cd {dir}
{ipr} -i {fasta} -o {out} -f tsv -appl Pfam -goterms -iprlookup -pa
echo -e "Protein_Accession\tSequence_Length\tPfam_ID\tDescription\tStart_location\tStop_location\tEvalue" > {base}.pfam.xls && cut -f 1,3,5-9 {out} >> {base}.pfam.xls
""".format(dir=tmp_dir,fasta=protein_faa,out=out_prefix+".gene_iprscan.tsv",base=out_prefix, ipr=ipr_sh)
    cmd_sh = open(os.path.join(working_dir,out_prefix+".ipr.sh"), 'w')
    cmd_sh.write(cmd_str)
    cmd_sh.close()
    os.system("sh %s" % (os.path.join(working_dir,out_prefix+".ipr.sh")))
    pfam_out = os.path.join(tmp_dir, out_prefix+".pfam.xls")

    # choose isplasmid contigs
    plasmid_contigs = os.path.join(working_dir, out_prefix+".plasmid.fasta")
    chromosome_contigs = os.path.join(working_dir, out_prefix+".chromosome.fasta")
    choosePlasmid(plasmid_pfamF, pfam_out, protein_gff, contigs, plasmid_contigs, chromosome_contigs)
    if os.path.exists(plasmid_contigs) and os.path.getsize(plasmid_contigs) > 0:
        sys.stderr.write("Info: %s has plasmid.\n" % out_prefix)
        sys.stderr.flush()
    else:
        sys.stderr.write("Info: %s don't have plasmid.\n" % out_prefix)
        sys.stderr.flush()

    #if os.path.getsize(plasmid_contigs) > 0:
    # remove temp files
    if os.path.exists(plasmid_contigs) and os.path.exists(chromosome_contigs) and os.path.getsize(chromosome_contigs)>0:
        os.remove(tmpcontigs)
        os.remove(os.path.join(working_dir,out_prefix+".ipr.sh"))
        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    main()
