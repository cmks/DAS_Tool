#!/usr/bin/env python
import argparse
import subprocess
import sys
import os
import pandas as pd

from Bio import SeqIO


BLAST6_HEADER = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]

USEARCH_HEADER = ["query", "target", "id", "alnlen", "ql", "tl", "mism", 
                  "opens", "qlo", "qhi", "tlo", "thi", "evalue", "bits"]

HEADERS = {
    "diamond": BLAST6_HEADER,
    "blast": BLAST6_HEADER,
    "usearch": USEARCH_HEADER
}


def make_db(search_engine, proteins, db_name):
    cmd = ""
    if search_engine == "diamond":
        cmd = [
            "diamond", "makedb",
            "--in", proteins,
            "-d", db_name
        ]
    elif search_engine == "blast":
        cmd = [
            "makeblastdb",
            "-in", proteins,
            "-dbtype", "prot",
            "-out", db_name
        ]
    elif search_engine == "usearch":
        cmd = [
            "usearch",
            "-makeudb_ublast",
            proteins,
            "-output", db_name
        ]
    else:
        print("don't support %s" % search_engine)
        sys.exit(-1)

    cmd_str = " ".join(cmd)
    print(cmd_str)
    subprocess.call(cmd_str, shell=True, stdout=sys.stdout, stderr=sys.stderr)


def run_blast(search_engine, query, database, blast_out, threads, first=True):
    cmd = ""
    if search_engine == "diamond":
        cmd = [
            "diamond", "blastp",
            "--query", query,
            "--db", database,
            "--outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
            "--out", blast_out,
            "--threads", threads
        ]
        if first:
            cmd += ["--evalue", "0.01", "--max-target-seqs", "0"]
        else:
            cmd += ["--evalue", "0.00001", "--max-target-seqs", "1"]

    elif search_engine == "blast":
        cmd = [
            "blastp",
            "-db", database,
            "-query", query,
            "-outfmt", "'6 qseqid sseqid pident length qlen slen evalue bitscore'",
            "-out", blast_out,
            "-num_threads", threads
        ]
        if first:
            cmd += ["-evalue", "0.01"]
        else:
            cmd += ["-evalue", "0.00001", "-max_target_seqs", "1"]

    elif search_engine == "usearch":
        cmd = [
            "usearch",
            "-ublast", query,
            "-db", database,
            "-threads", threads,
            "-userout", blast_out,
            "-userfields", "query+target+id+alnlen+ql+tl+mism+opens+qlo+qhi+tlo+thi+evalue+bits"
        ]
        if first:
            cmd += ["-evalue", "0.01"]
        else:
            cmd += ["-evalue", "0.00001", "-maxhits", "1"]

    else:
        print("don't support %s" % search_engine)
        sys.exit(-1)

    cmd_str = " ".join(cmd)
    print(cmd_str)
    subprocess.call(cmd_str, shell=True, stdout=sys.stdout, stderr=sys.stderr)


def extract_seq(faa, blast_out, faa_out):
    gene_id_list = []
    with open(blast_out, 'r') as ih:
        for line in ih:
            gene_id_list.append(line.strip().split("\t")[1])

    with open(faa_out, 'w') as oh:
        for seq in SeqIO.parse(faa, "fasta"):
            if seq.id in gene_id_list:
                SeqIO.write(seq, oh, "fasta")


def find_scg(search_engine, faa, scg_faa, db, all_db, scg_lookup,
             output_dir, threads):
    blast_out = os.path.join(
        output_dir, os.path.basename(faa) + ".findSCG.b6")
    run_blast(search_engine, scg_faa, db, blast_out, threads, True)

    # extract proteins
    scg_faa_ = os.path.join(
        output_dir, os.path.basename(faa) + ".scg.candiates.faa")
    extract_seq(faa, blast_out, scg_faa_)

    # verify SCGs by blasting against all proteins of all genomes
    blast_out = os.path.join(output_dir, os.path.basename(faa) + ".all.b6")
    run_blast(search_engine, scg_faa_, all_db, blast_out, threads, False)

    # scg annotation
    scg_lookup_df = pd.read_csv(scg_lookup, sep='\t', names=[
                                "scg_id", "scg_annotation"])\
                      .astype({"scg_id": str})

    # parse blast out
    blast_out_df = pd.read_csv(blast_out, sep='\t', names = HEADERS[search_engine])\
                     .rename(columns={
                         "qseqid": "gene_id",
                         "sseqid": "scg_id",
                         "query": "gene_id",
                         "target": "scg_id",
                         "alnlen": "length",
                         "tl": "slen"})\
                     .astype({"scg_id": str,
                              "length": int,
                              "slen": int})

    scg_df = blast_out_df[blast_out_df.scg_id.isin(scg_lookup_df.scg_id)]\
        .query('length > slen * 0.5')\
        .merge(scg_lookup_df)

    if scg_df.empty:
        return None
    else:
        scg = os.path.join(output_dir, os.path.basename(faa) + ".scg")
        scg_df.loc[:, ["gene_id", "scg_annotation"]].to_csv(
            scg, sep='\t', index=False, header=None)
        return scg


def set_db_name(search_engine, faa, output_dir):
    db_suffix = {
        "diamond": ".dmnd",
        "blast": ".blastdb",
        "usearch": ".udb"
    }
    return os.path.join(output_dir, os.path.basename(faa) + db_suffix[search_engine])


def main():
    parser = argparse.ArgumentParser(description="find single copy gene")
    parser.add_argument("--search-engine", choices=["diamond", "blast", "usearch"],
                        dest="search_engine", help="search engine")
    parser.add_argument("--proteins", dest="faa", help="gene")
    parser.add_argument("--all-proteins",
                        dest="all_faa", help="all gene")
    parser.add_argument("--scg-proteins",
                        dest="scg_faa", help="single copy gene")
    parser.add_argument("--lookup", dest="scg_lookup",
                        help="single copy gene annotation")
    parser.add_argument("--threads", dest="threads", help="threads", default=7)
    parser.add_argument("--output-dir", dest="output_dir",
                        help="output directory")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # build db for faa and all_faa
    db = set_db_name(args.search_engine, args.faa, args.output_dir)
    all_db = set_db_name(args.search_engine, args.all_faa, args.output_dir)
    make_db(args.search_engine, args.faa, db)
    make_db(args.search_engine, args.all_faa, all_db)

    # find single copy gene
    find_scg(args.search_engine, args.faa, args.scg_faa, db, all_db,
             args.scg_lookup, args.output_dir, args.threads)


if __name__ == '__main__':
    main()
