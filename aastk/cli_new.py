#!/usr/bin/env python3
import argparse
from contextlib import contextmanager

@contextmanager
def subparser(parser, name, description):
    yield parser.add_parser(name, conflict_handler='resolve', help=description)

@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)

def __block(group, required=False):
    group.add_argument('-b', '--block', type=int, required=required,
                       help='Choose diamond blastp sequence block size in billions of letters. (Default: 6)')

def __chunk(group, required=False):
    group.add_argument('-c', '--chunk', type=int, required=required,
                       help='Choose number of chunks for diamond blastp index processing. (Default: 2)')

def __db(group, required=False):
    group.add_argument('-d', '--db', type=str, required=required,
                       help='Specify path to DIAMOND database')

def __extracted(group, required=False):
    group.add_argument('-e', '--extracted', type=str, required=required,
                       help='Path to FASTA file containing extracted matching sequences')

def __key_column(group, required=False):
    group.add_argument('-k', '--key_column', type=int, default=0, required=required,
                       help='Column index in the BLAST tab file to pull unique IDs from (default is 0)')

def __matrix(group, required=False):
    group.add_argument('-m', '--matrix', choices=['BLOSUM45', 'BLOSUM62'], required=required,
                       help='Choose BLOSUM substitution matrix (BLOSUM 45 or BLOSUM 62)')

def __max_scores(group, required=False):
    group.add_argument('-m', '--max_scores', type=str, required=required,
                       help='Path to file containing max self scores')

def __output(group, required=False):
    group.add_argument('-o', '--output', type=str, required=required,
                       help='Desired output location (default: current working directory)')

def __protein_name(group, required=False):
    group.add_argument('-p', '--protein_name', type=str, default=None, required=required,
                       help='Name of protein of interest')

def __query(group, required=False):
    group.add_argument('-q', '--query', type=str, default=None, required=required,
                       help='Path to query FASTA')

def __seed(group, required=False):
    group.add_argument('-s', '--seed', type=str, default=None, required=required,
                       help='Path to FASTA files containing seed sequences for database creation')

def __sensitivity(group, required=False):
    group.add_argument('--sensitivity', choices=['fast', 'sensitive', 'mid-sensitive', 'very-sensitive', 'ultra-sensitive', 'faster'], required=required,
                       help='Set the sensitivity level for the DIAMOND search: fast, sensitive, mid-sensitive, very-sensitive, '
                            'ultra-sensitive, or faster (default: fast)')

def __tabular(group, required=False):
    group.add_argument('-t', '--tabular', type=str, default=None, required=required,
                       help='Path to tabular BLAST/DIAMOND output file')

def __threads(group, required=False):
    group.add_argument('-n', '--threads', type=int, default=1, required=required,
                       help='Number of threads to be used (default: 1)')

def get_main_parser():
    main_parser = argparse.ArgumentParser(
        prog='aastk', add_help=False, conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    with subparser(sub_parsers, 'build', 'Build DIAMOND database from seed sequence(s)') as parser:
        with arg_group(parser, 'Required arguments') as grp:
            __protein_name(grp, required=True)
            __seed(grp, required=True)
            __db(grp, required=True)
        with arg_group(parser, 'Optional') as grp:
            __threads(grp)

    with subparser(sub_parsers, 'search', 'Search DIAMOND reference database for homologous sequences') as parser:
        with arg_group(parser, 'Required arguments') as grp:
            __db(grp, required=True),
            __query(grp, required=True),
            __protein_name(grp, required=True),
        with arg_group(parser, 'Optional') as grp:
            __output(grp)
            __threads(grp)
            __block(grp)
            __chunk(grp)
            __sensitivity(grp)

    with subparser(sub_parsers, 'extract', 'Extract reads that have DIAMOND hits against custom database') as parser:
        with arg_group(parser, 'Required arguments') as grp:
            __tabular(grp, required=True),
            __query(grp, required=True),
            __output(grp, required=True)
        with arg_group(parser, 'Optional') as grp:
            __key_column(grp)

    with subparser(sub_parsers, 'calculate', 'Calculate max scores for extracted sequences using BLOSUM matrix') as parser:
        with arg_group(parser, 'Required arguments') as grp:
            __extracted(grp, required=True)
            __matrix(grp, required=True)
        with arg_group(parser, 'Optional') as grp:
            __output(grp)

    with subparser(sub_parsers, 'bsr', 'Compute BSR (Blast Score Ratio) using a BLAST tab file and max scores from a TSV.') as parser:
        with arg_group(parser, 'Required arguments') as grp:
            __tabular(grp, required=True)
            __max_scores(grp, required=True)
        with arg_group(parser, 'Optional') as grp:
            __output(grp)
            __key_column(grp)

    with subparser(sub_parsers, 'pasr', 'PASR: protein alignment score ratio') as parser:
        with arg_group(parser, 'Required arguments') as grp:
            __matrix(grp, required=True)
            __protein_name(grp, required=True)
            __query(grp, required=True)
            __seed(grp, required=True)
        with arg_group(parser, 'Optional') as grp:
            __db(grp)
            __output(grp)
            __threads(grp)
            __key_column(grp)
            __block(grp)
            __chunk(grp)
            __sensitivity(grp)

    return main_parser
