#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, re
import argparse
from collections import deque
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


ScriptPath = os.path.split(os.path.realpath(sys.argv[0]))[0]
lib_path = os.path.join(ScriptPath, "..", "lib")
sys.path.append(lib_path)
import FileUtilx

## Option params ##
def ParseArgs():
    usage = "python3 %(prog)s [options] --fasta assembly.fasta -d gene_prodict_outdir"
    parser = argparse.ArgumentParser(usage=usage, description="%(prog)s v1.0.0")
    parser.add_argument("--fasta", dest="fasta", action="store", help="the assembly fasta. required.")
    parser.add_argument("-d", "--dir", dest="gene_prodict_outdir", action="store", help="out directory of gene_prodict, a new directory. required.")
    parser.add_argument("--name", dest="name", help="the analysis number as sequence id(prefix). required.")
    parser.add_argument("-f", dest="format", action="store", default="gff", help="Select gene_prodict output format (gbk, gff, or sco). default is gff")
    parser.add_argument("--max_olap", dest="max_olap", action="store", type=int,default=50, help="Set maximum overlap length to <n>.Overlaps this short or shorter are ignored. default 50")
    parser.add_argument("--threshold", dest="threshold", action="store", type=int, default=30, help="Set threshold score for calling as gene to n.  If the in-frame score >= <n>, then the region is given a number and considered a potential gene. default 30")
    parser.add_argument("--gene_len", dest="gene_len", action="store", type=int,default=110, help="Set minimum gene length to <n>, default is 110")
    parser.add_argument("--start_codons" ,dest="start_codons", action="store", default= "atg,gtg,ttg", help="Use comma-separated list of codons as start codons, default is 'atg,gtg,ttg'")
    parser.add_argument("--tr_table", dest="tr_table", action="store", type=int,default=11, help="Specify a translation table to use, default is 11")
    parser.add_argument("--conf", dest="conf", default=os.path.join(ScriptPath, "..", "pipeline_env.conf"),help="config file of this script.")
    parser.add_argument("--continue", dest="continues", action="store_true", help="continue run assemble")

    options = parser.parse_args()
    if not options.gene_prodict_outdir or not options.fasta or not options.name:
        parser.print_help()
        sys.exit(1)
    return options

## gene_prodict ##
def prodictProdigal(fasta, outdir):
    global running_status
    fasta = os.path.abspath(fasta)
    outdir = os.path.abspath(outdir)
    if os.path.exists(outdir):
        sys.stderr.write("Warning: the %s directory is exists.\n" % outdir)
    else:
        os.makedirs(outdir)

    gene_prodict_outd = os.path.join(outdir, options.name+'.gene.'+options.format)
    cmd_sh = "{0} -i {1} -o {2} -f {3} -g {4}\n".format(pipe_confs["prodigal"], fasta, gene_prodict_outd, options.format, options.tr_table)
    shFile = open(os.path.join(outdir, "run_prodigal.sh"), 'w')
    shFile.write(cmd_sh)
    shFile.close()
    status = os.system("sh %s 2> %s" % (os.path.join(outdir, "run_prodigal.sh"), os.path.join(outdir, "run_prodigal.log")))
    if status == 0:
        sys.stderr.write("gene prodict finished.\n")
        running_status = True
        return gene_prodict_outd
    else:
        sys.stderr.write("gene prodict failed.\n")
        running_status = False
        sys.exit(1)


def prodictGlimmer(fasta, outdir):
    global running_status
    fasta = os.path.abspath(fasta)
    outdir = os.path.abspath(outdir)

    if os.path.exists(outdir):
        sys.stderr.write("Warning: the %s directory is exists.\n" % outdir)
    else:
        os.makedirs(outdir)

    cmd_sh = """cd {0} &&\\
{1} {2} prodict1 &&\\
{3} -o {4} -g {5} -t {6} -A {7} -l -z {8} -b {9} {2} {10} prodict2
""".format(outdir, g3_iterated_csh, fasta, glimmer3, options.max_olap, options.gene_len, options.threshold, options.start_codons, 1, "prodict1.motif", "prodict1.icm")
    shFile = open(os.path.join(outdir, "run_glimmer3.sh"), 'w')
    shFile.write(cmd_sh)
    shFile.close()
    status = os.system("sh %s 2> %s" % (os.path.join(outdir, "run_glimmer3.sh"),os.path.join(outdir, "run_glimmer3.log")))
    prodict2_predict = os.path.join(outdir, "prodict2.predict")
    if status == 0:
        sys.stderr.write("gene prodict finished.\n")
        running_status = True
        return prodict2_predict
    else:
        sys.stderr.write("gene prodict failed.\n")
        running_status = False
        sys.exit(1)

## get the result of predict ##
def Prodigal_Parse(fastaf, gfff, outPrefix, codon=11):
    fastaf = os.path.abspath(fasta)
    gfff = os.path.abspath(gfff)

    # read fasta sequences
    seq_Dicts = {}
    for seq_record in SeqIO.parse(fastaf, 'fasta'):
            seq_id = seq_record.id
            seq_Dicts[seq_id] = seq_record

    # read the gff file of prodigal result
    number = 1
    nuc_out = open(outPrefix+".gene.fnn", 'w')
    protein_out = open(outPrefix+".protein.faa", 'w')
    gff_out = open(outPrefix+".gene.gff", 'w')
    gff_out.write("##gff-version  3\n")
    with open(gfff, 'r')as fh:
        for line in fh:
            if re.search(r"^#|^\s*$", line):
                continue
            if re.search(r"CDS", line):
                fields = line.strip().split("\t")
                scaf_id = fields[0]
                # filter short genes
                if (int(fields[4]) - int(fields[3]) + 1) < options.gene_len:
                    continue
                start_site = int(fields[3]) - 1
                end_site = int(fields[4])
                substr = seq_Dicts[scaf_id].seq[start_site:end_site]
                subrecord = SeqRecord(substr, id=scaf_id+"_orf"+str(number), description='')
                if fields[6] == '-':
                    rc_seq = subrecord.seq.reverse_complement()
                    subrecord = SeqRecord(rc_seq, id=subrecord.id, description='')
                #输出核酸序列
                SeqIO.write(subrecord, nuc_out, 'fasta')
                #protein
                protein = subrecord.seq.translate(table=codon)
                protein_record = SeqRecord(Seq(str(protein).rstrip("*"), IUPAC.protein), id=subrecord.id, description='')
                SeqIO.write(protein_record, protein_out, 'fasta')
                # gff out
                gff_out.write("\t".join(fields[:2]+['gene']+fields[3:8])+"\tlocus_tag="+subrecord.id+"\n")
                gff_out.write("\t".join(fields[:8])+"\tlocus_tag="+subrecord.id+";translation="+str(protein).rstrip("*")+"\n")
                number += 1

    nuc_out.close()
    protein_out.close()
    gff_out.close()


# main #
if __name__ == "__main__":
    # get params
    options = ParseArgs()
    fasta = os.path.abspath(options.fasta)

    # config file and return a dict of configs
    pipe_confs = FileUtilx.read_Param(os.path.abspath(options.conf))
    # software pathway
    glimmer3 = os.path.join(pipe_confs["glimmer3"], "bin", "glimmer3")
    g3_iterated_csh = os.path.join(pipe_confs["glimmer3"], "scripts", "g3-iterated.csh")

    # output directory
    working_out_dir = os.path.abspath(options.gene_prodict_outdir)
    if not os.path.exists(working_out_dir):
        os.makedirs(working_out_dir)
    else:
        sys.stderr.write("Warning: The %s directory exists.\n" % working_out_dir)
        sys.stderr.flush()

    # the log file of genome assemble
    if options.continues:
        gene_predict_Log = open(os.path.join(working_out_dir, "gene_predict.log"), 'a')
    else:
        gene_predict_Log = open(os.path.join(working_out_dir, "gene_predict.log"), 'w')

    # record script running status
    running_status = True

    # GC opinion 
    GC_file = os.path.join(working_out_dir,"GC_file.txt")
    status = os.system("perl %s %s > %s" % (os.path.join(ScriptPath, "GC_count.pl"), fasta, GC_file))
    tail_1line = deque(open(GC_file), 1)
    k = float(tail_1line[0].strip())
    """
    with open(GC_file, 'r') as fh:
        for line in fh:
            if re.search(r"^[0-9]", line):
                k = float(line.strip())
    """

    # gene predict
    tmp_out_dir = working_out_dir
    #tmp_out_dir = os.path.join(working_out_dir, "gene_predict_out")

    if k > 70:
        gene_predict_Log.write("gene predict toolkit: prodigal.\n")
        gene_predict_Log.write("start gene predict...\n")
        test_peodict_file = prodictProdigal(fasta, tmp_out_dir) # gene predict result
        out_prefix = os.path.join(tmp_out_dir, options.name)
        status = Prodigal_Parse(fasta, test_peodict_file, out_prefix, options.tr_table) # the result of predict
        gene_file = out_prefix + ".gene.fnn"
        protein_file = out_prefix + ".protein.faa"
        if running_status:
            gene_predict_Log.write("finish gene predict...\n")
        else:
            gene_predict_Log.write("the gene predict failed...\n")
            sys.exit(1)
    else:
        gene_predict_Log.write("gene predict toolkit: glimmer.\n")
        gene_predict_Log.write("start gene predict...\n")
        out_prefix = os.path.join(tmp_out_dir, options.name)
        test_peodict_file = prodictGlimmer(fasta, tmp_out_dir) # gene predict result
        status = os.system("perl %s %s %s %s %s" % (ScriptPath+os.sep+"glimmer_predict_parse.pl", test_peodict_file, fasta, out_prefix, options.tr_table))
        gene_file = out_prefix + ".gene.fnn"
        protein_file = out_prefix + ".protein.faa"
        if status == 0:
            running_status = True
        else:
            running_status = False
        if running_status:
            gene_predict_Log.write("finish gene predict...\n")
        else:
            gene_predict_Log.write("the gene predict failed...\n")
            sys.exit(1)

    # summary and analyze
    if running_status:
        status = os.system("perl %s %s %s %s" % (ScriptPath+os.sep+"summary_gene.pl", fasta, gene_file, os.path.join(working_out_dir, options.name+".gene_summary.xls")))
        status = os.system("python3 %s %s %s" % (ScriptPath+os.sep+"gene_histPlot_v1.1.py", gene_file, os.path.join(working_out_dir, "gene")))
        if status == 0:
            gene_predict_Log.write("gene predict finished.\n")
            gene_predict_Log.write("All completed.\n")
        else:
            gene_predict_Log.write("gene predict failed.\n")

