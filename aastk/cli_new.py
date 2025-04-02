#!/usr/bin/env python3
import argparse
import tempfile
from contextlib import contextmanager

@contextmanager
def subparser(parser, name, description):
    yield parser.add_parser(name, conflict_handler='resolve', help=description)

@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)

def __matrix(group, required=False):
    group.add_argument('-m', '--matrix', choices=['BLOSUM45', 'BLOSUM62'],
                       help='Choose BLOSUM substitution matrix (BLOSUM 45 or BLOSUM 62)')

def __protein_name(group, required=False):
    group.add_argument('-p', '--protein_name', type=str, default=None, required=required,
                       help='Name of protein of interest')

def __query_fasta(group, required=False):
    group.add_argument('-q', '--query', type=str, default=None, required=required,
                       help='Path to query FASTA')

def __seed_fasta(group, required=False):
    group.add_argument('-s', '--seeds', type=str, default=None, required=required,
                       help='Path to FASTA files containing seed sequences for database creation')

def __threads(group, required=False):
    group.add_argument('-t', '--threads', type=int, default=1, required=required,
                       help='Number of threads to be used (default: 1)')






def get_main_parser():
    main_parser = argparse.ArgumentParser(
        prog='aastk', add_help=False, conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    with subparser(sub_parsers, 'pasr', 'PASR: protein alignment score ratio') as parser:
        with arg_group(parser, 'Required arguments') as grp:
            __matrix(grp, required=True)
            __protein_name(grp, required=True)
            __query_fasta(grp, required=True)
            __seed_fasta(grp, required=True)
        with arg_group(parser, 'Optional') as grp:
            __threads(grp, required=False)

    return main_parser
