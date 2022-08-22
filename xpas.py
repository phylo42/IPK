#!/usr/bin/env python3

"""
XPAS wrapper script.
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
import click
import subprocess
from pathlib import Path
import shutil
import json


@click.group()
def xpas():
    """
    xpas

    N. Romashchenko, B. Linard, F. Pardi, E. Rivals
    """
    pass


# See: https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model
NUCL_MODELS = ['JC', 'K80', 'F81', 'HKY', 'TN93ef',
               'TN93', 'K81', 'K81uf', 'TPM2', 'TPM2uf', 'TPM3', 'TPM3uf',
               'TIM1', 'TIM1uf', 'TIM2', 'TIM2uf', 'TIM3', 'TIM3uf', 'TVMef',
               'TVM', 'SYM', 'GTR']
AMINO_MODELS = ['Blosum62', 'cpREV', 'Dayhoff', 'DCMut', 'DEN', 'FLU', 'HIVb', 'HIVw', 'JTT',
                'JTT-DCMut', 'LG', 'mtART', 'mtMAM', 'mtREV', 'mtZOA', 'PMB', 'rtREV', 'stmtREV',
                'VT', 'WAG', 'LG4M', 'LG4X', 'PROTGTR']
ALL_MODELS = NUCL_MODELS + AMINO_MODELS


KMER_FILTERS = ["no-filter", "mif0", "mif1", "random"]


def validate_filter(ctx, param, value):
    value = value.lower()
    if value not in KMER_FILTERS:
        valid_values = ', '.join(v for v in KMER_FILTERS)
        raise click.BadParameter('Filter must be one of: ' + valid_values)
    return value


def validate_model(ctx, param, value):
    if ('ar_config' in ctx.params) or (value and value in ALL_MODELS):
        return value

    raise click.BadParameter(f'Please define a valid evolutionary model either via --model or in a config file '
                             f'via --ar-config. Valid values: {ALL_MODELS}')

@xpas.command()
@click.option('-b', '--ar',
              type=click.Path(exists=True),
              required=False,
              help="Path for the ancestral reconstruction software used. Currently,"
                   "PhyML and RAxML-ng are supported. If not given, XPAS will try to"
                   "find RAxML-ng using environmental variables.")
@click.option('-r', '--refalign',
              type=click.Path(exists=True),
              required=True,
              help="""Reference alignment in fasta format.
                  It must be the multiple alignment used to build the reference tree.""")
@click.option('-t', '--reftree',
              type=click.Path(exists=True),
              required=True,
              help=" Reference tree in the newick format")
@click.option('-s', '--states',
              type=click.Choice(['nucl', 'amino']),
              default='nucl', show_default=True,
              required=True,
              help="States used in analysis.")
@click.option('-v', '--verbosity',
              type=int,
              default=0, show_default=True,
              help="Verbosity level: -1=none ; 0=default ; 1=high")
@click.option('-w', '--workdir',
              required=True,
              type=click.Path(dir_okay=True, file_okay=False),
              help="Working directory for temp files.")
@click.option('--write-reduction',
              type=click.Path(file_okay=True, dir_okay=False),
              help="Write reduced alignment to file.")
@click.option('-a', '--alpha',
              type=float,
              default=1.0, show_default=True,
              help="Gammma shape parameter used in ancestral reconstruction.")
@click.option('-c', '--categories',
              type=int,
              default=4, show_default=True,
              help="Number of categories used in ancestral reconstruction.")
@click.option('-k', '--k',
             type=int,
             default=8, show_default=True,
             help="k-mer length used at DB build.")
@click.option('-m', '--model', type=click.UNPROCESSED, callback=validate_model, required=False,
             help="Model used in AR, one of the following:\n"
                  f"nucl: {', '.join(x for x in NUCL_MODELS)}\n"
                  f"amino: {', '.join(x for x in AMINO_MODELS)}")
@click.option('--convert-uo',
              is_flag=True,
              help="U, O amino acids are converted to C, L.")
@click.option('--no-reduction',
              is_flag=True,
              help="""Do not operate alignment reduction. This will 
                  keep all sites of input reference alignment and 
                  may produce erroneous ancestral k-mers.""")
@click.option('--reduction-ratio',
              type=float,
              default=0.99, show_default=True,
              help="""Ratio for alignment reduction, e.g. sites 
                holding >99% gaps are ignored.""")
@click.option('--omega',
              type=float,
              default=1.5, show_default=True,
              help="""Modifier levelling the threshold used during
                  phylo-kmer filtering, T=(omega/#states)^k""")
@click.option('--filter',
              callback=validate_filter,
              default="no-filter", show_default=True)
@click.option('-f', type=float, default=1.0)
@click.option('-u', '--mu',
              type=float,
              default=0.5, show_default=True,
              help="""K-mer filter threshold""")
@click.option('--use-unrooted',
              is_flag=True,
              help="""Confirms you accept to use an unrooted reference
                  tree (option -t). The trifurcation described by the
                  newick file will be considered as root. Be aware that
                  meaningless roots may impact accuracy.""")
@click.option('--merge-branches',
              is_flag=True,
              default=False,
              help="""Builds a databases of phylo k-mers, merging phylo
                  k-mers of different branches. Thus, for every k-mer
                  the only one maximum score among all the branches
                  will be saved.""")
@click.option('--ar-dir',
             type=click.Path(exists=True, dir_okay=True, file_okay=False),
             help="""Skip ancestral sequence reconstruction, and 
                  uses outputs from the specified directory.""")
@click.option('--ar-only',
             is_flag=True,
             default=False, show_default=True,
             help="Dev option. Run only ancestral reconstruction and tree extension. No database will be built")
@click.option('--ar-config',
              required=False,
              type=click.Path(exists=True),
              help="A .json-formatted config file for ancestral reconstruction parameters. "
                   "See xpas.readthedocs.org for help")
@click.option('--keep-positions',
              is_flag=True,
              default=False,
              help="""Keeps phylo k-mers positions in the alignment. Makes databases larger in size.""")
@click.option('--uncompressed',
              is_flag=True,
              default=False,
              help="""Stores databases uncompressed.""")
@click.option('--threads',
             type=int,
             default=4, show_default=True,
             help="Number of threads used.")
def build(ar,
          refalign, reftree, states,
          verbosity,
          workdir, write_reduction, #dbfilename,
          alpha, categories, #ghosts,
          k, model, convert_uo, #gap_jump_thresh,
          no_reduction, reduction_ratio, omega,
          filter, f, mu, use_unrooted, merge_branches,
          ar_dir, ar_only, ar_config,
          keep_positions, uncompressed,
          threads):
    """
    Builds a database of phylo-k-mers.

    Minimum usage:

    \tpython xpas.py build -s [nucl|amino] -b ARbinary -w workdir -r alignment.fasta -t tree.newick

    """
    build_database(ar,
                   refalign, reftree, states,
                   verbosity,
                   workdir, write_reduction, #dbfilename,
                   alpha, categories, #ghosts,
                   k, model, convert_uo, #gap_jump_thresh,
                   no_reduction, reduction_ratio, omega,
                   filter, f, mu, use_unrooted, merge_branches,
                   ar_dir, ar_only, ar_config,
                   keep_positions, uncompressed,
                   threads)


def find_raxmlng():
    """Try to find RAxML-ng"""
    path = shutil.which("raxml-ng")
    if not path:
        raise RuntimeError("RAxML-ng not found. Please check it exists in your PATH or provide a full filename")
    return path


def parse_config(ar_config: str) -> str:
    """
    Parse .json-formatted config for ancestral reconstruction parameters
    """
    with open(ar_config, 'r') as f:
        content = json.load(f)
        if "arguments" not in content:
            raise RuntimeError(f"Error parsing {ar_config}: 'arguments' not found")
        arguments = content["arguments"]
        return " ".join(f"--{param} {value}" for param, value in arguments.items())


def build_database(ar,
                   refalign, reftree, states,
                   verbosity,
                   workdir, write_reduction, #dbfilename,
                   alpha, categories, #ghosts,
                   k, model, convert_uo, #gap_jump_thresh,
                   no_reduction, reduction_ratio, omega,
                   filter, f, mu, use_unrooted, merge_branches,
                   ar_dir, ar_only, ar_config,
                   keep_positions, uncompressed,
                   threads):

    if not ar:
        ar = find_raxmlng()

    # create working directory
    Path(workdir).mkdir(parents=True, exist_ok=True)
    current_dir = os.path.dirname(os.path.realpath(__file__))

    ar_parameters = parse_config(ar_config) if ar_config else ""

    # run rappas2
    if states == 'nucl':
        if keep_positions:
            raise RuntimeError("--keep-positions is not supported for DNA.")
        else:
            bin = f"{current_dir}/bin/xpas/xpas-dna"
    else:
        if keep_positions:
            raise NotImplementedError()
            bin = f"{current_dir}/bin/xpas/xpas-aa-pos"
        else:
            bin = f"{current_dir}/bin/xpas/xpas-aa"


    command = [
        bin,
        "--ar-binary", str(ar),
        "--refalign", str(refalign),
        "-t", str(reftree),
        "-w", str(workdir),
        "-k", str(k),
        #"--model", model,
        "--alpha", str(alpha),
        "--categories", str(categories),
        "--reduction-ratio", str(reduction_ratio),
        "-o", str(omega),
        "--" + filter.lower(),
        "-u", str(mu),
        "-j", str(threads)
    ]

    if model:
        command.append("--model")
        command.append(model)
    if ar_only:
        command.append("--ar-only")
    if ar_dir:
        command.append("--ar-dir")
        command.append(ar_dir)
    if ar_parameters:
        command.append("--ar-parameters")
        command.append(f'"{ar_parameters}"')
    if no_reduction:
        command.append("--no-reduction")
    if merge_branches:
        command.append("--merge-branches")
    if use_unrooted:
        command.append("--use-unrooted")
    if uncompressed:
        command.append("--uncompressed")

    # remove the temporary folder just in case
    hashmaps_dir = f"{workdir}/hashmaps"
    #subprocess.call(["rm", "-rf", hashmaps_dir])

    command_str = " ".join(s for s in command)
    print(command_str)
    p = subprocess.run(command_str, shell=True, check=True)

    # clean after
    subprocess.call(["rm", "-rf", hashmaps_dir])

    if p.returncode != 0:
        raise RuntimeError(f"XPAS returned error: {return_code}")


if __name__ == "__main__":
    xpas()
