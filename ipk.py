#!/usr/bin/env python3

"""
IPK wrapper script.
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
import click
import subprocess
from pathlib import Path
import shutil
import json


# See: https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model
NUCL_MODELS = ['JC', 'K80', 'F81', 'HKY', 'TN93ef',
               'TN93', 'K81', 'K81uf', 'TPM2', 'TPM2uf', 'TPM3', 'TPM3uf',
               'TIM1', 'TIM1uf', 'TIM2', 'TIM2uf', 'TIM3', 'TIM3uf', 'TVMef',
               'TVM', 'SYM', 'GTR']
AMINO_MODELS = ['Blosum62', 'cpREV', 'Dayhoff', 'DCMut', 'DEN', 'FLU', 'HIVb', 'HIVw', 'JTT',
                'JTT-DCMut', 'LG', 'mtART', 'mtMAM', 'mtREV', 'mtZOA', 'PMB', 'rtREV', 'stmtREV',
                'VT', 'WAG', 'LG4M', 'LG4X', 'PROTGTR']
ALL_MODELS = NUCL_MODELS + AMINO_MODELS


KMER_FILTERS = ["mif0", "random"]

GHOST_STRATEGIES = ["inner-only", "outer-only", "both"]


@click.group()
def ipk():
    """
    IPK

    N. Romashchenko, B. Linard, F. Pardi, E. Rivals
    """
    pass


def validate_filter(ctx, param, value):
    value = value.lower()
    if value not in KMER_FILTERS:
        valid_values = ', '.join(v for v in KMER_FILTERS)
        raise click.BadParameter('Filter must be one of: ' + valid_values)
    return value

def validate_ghosts(ctx, param, value):
    value = value.lower()
    if value not in GHOST_STRATEGIES:
        valid_values = ', '.join(v for v in GHOST_STRATEGIES)
        raise click.BadParameter('Strategy must be one of: ' + valid_values)
    return value


def validate_model(ctx, param, value):
    if ('ar_config' in ctx.params) or (value and value in ALL_MODELS):
        return value

    raise click.BadParameter(f'Please define a valid evolutionary model either via --model or in a config file '
                             f'via --ar-config. Valid values: {ALL_MODELS}')


@ipk.command()
@click.option('-b', '--ar',
              type=click.Path(exists=True),
              required=False,
              help="""Path for the ancestral reconstruction software to use to infer ancestral states. 
                      Currently, :program:`PhyML` and :program:`RAxML-ng` are supported.
                      If no value is provided, IPK searches for ``raxml-ng`` in your PATH.""")
@click.option('-r', '--refalign',
              type=click.Path(exists=True),
              required=True,
              help="""Reference multiple sequence alignment in FASTA format. 
              Reference sequences should be in 1-to-1 correspondance with the tips 
              of the tree given by :option:`--reftree`.""")
@click.option('-t', '--reftree',
              type=click.Path(exists=True),
              required=True,
              help="""Reference phylogenetic tree in Newick format. 
              The tips should be in 1-to-1 correspondance with the sequences given by :option:`--refalign`. 
              If the tree is unrooted, the flag :option:`--use-unrooted` should be used.""")
@click.option('-s', '--states',
              type=click.Choice(['nucl', 'amino']),
              default='nucl', show_default=True,
              required=True,
              help="States used in the analysis, either `nucl` for DNA or `amino` for proteins.")
@click.option('-v', '--verbosity',
              type=int,
              default=0, show_default=True,
              help="Verbosity level: -1=none ; 0=default ; 1=high")
@click.option('-w', '--workdir',
              required=True,
              type=click.Path(dir_okay=True, file_okay=False),
              help="Working directory for temporary files.")
@click.option('--write-reduction',
              type=click.Path(file_okay=True, dir_okay=False),
              help="Flag indicating to write reduced alignment to file.")
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
             help="k-mer length used. Must be a value in [2, 31].")
@click.option('-m', '--model', type=click.UNPROCESSED, callback=validate_model, required=False,
             help="Phylogenetic model used for ancestral state reconstruction.\n\n"
                  "Must be one of the following:\n\n"
                  f"nucl: {', '.join(x for x in NUCL_MODELS)}\n\n"
                  f"amino: {', '.join(x for x in AMINO_MODELS)}")
@click.option('--convert-uo',
              is_flag=True,
              help="Convert U, O amino acids to C, L.")
@click.option('--no-reduction',
              is_flag=True,
              help="""Disables alignment reduction. This will 
                  keep all sites of the input reference alignment and 
                  may produce erroneous ancestral k-mers.""")
@click.option('--reduction-ratio',
              type=float,
              default=0.99, show_default=True,
              help="""Ratio for alignment reduction. 
              Sites holding a higher percentage of gaps than this value will be removed.""")
@click.option('--omega',
              type=float,
              default=1.5, show_default=True,
              help="""Score threshold modifier. Determines the 
              minimal phylo-k-mer score considered according to the following formula:
              (omega / #states)^k""")
@click.option('--filter',
              callback=validate_filter,
              default="mif0", show_default=True,
              help="""Filtering function used to filter phylo-k-mers. Currently supported values:
              mif0, random.""")
@click.option('-u', '--mu',
              type=float,
              default=1.0, show_default=True,
              help="""K-mer filtering threshold. Determines the fraction of most informative k-mers 
              that will be saved in the resulting database.""")
@click.option('--ghosts',
              callback=validate_ghosts,
              default="both", show_default=True,
              help="""The strategy to process ghost nodes. Supported values:
              inner-only, outer-only, both""")
@click.option('--use-unrooted',
              is_flag=True,
              help="""Confirms you accept to use an unrooted reference
                  tree provided with :option:`-t`. The trifurcation described by the
                  newick file will be considered as the root. WARNING: meaningless rooting
                  can influence the accuracy of phylo-k-mer-based applications.""")
@click.option('--merge-branches',
              is_flag=True,
              default=False,
              help="""Creates a database by merging phylo-k-mers of different branches
              together. Thus, for every k-mer, the only one maximum score 
              among all the branches will be stored.""")
@click.option('--ar-dir',
             type=click.Path(exists=True, dir_okay=True, file_okay=False),
             help="""Skip ancestral sequence reconstruction
             and use outputs from the specified directory.""")
@click.option('--ar-only',
             is_flag=True,
             default=False, show_default=True,
             help="""If set, IPK stops after tree extension and ancestral reconstruction. 
             No database will be created.""")
@click.option('--ar-config',
              required=False,
              type=click.Path(exists=True),
              help="""A .json-formatted config file for ancestral reconstruction parameters.
                   See phylo-k-mers.readthedocs.org for help""")
@click.option('--keep-positions',
              is_flag=True,
              default=False,
              help="""Keeps alignment positions where phylo-k-mers obtain
              reported scores. Makes output databases larger in size.""")
@click.option('--uncompressed',
              is_flag=True,
              default=False,
              help="""Disables database compression. 
              Saves time during database load and save, but requires more disk space.""")
@click.option('--threads',
             type=int,
             default=4, 
             show_default=True,
             help="Number of threads used to compute phylo-k-mers.")
def build(ar,
          refalign, reftree, states,
          verbosity,
          workdir, write_reduction, #dbfilename,
          alpha, categories,
          k, model, convert_uo, #gap_jump_thresh,
          no_reduction, reduction_ratio, omega,
          filter, mu, ghosts, use_unrooted, merge_branches,
          ar_dir, ar_only, ar_config,
          keep_positions, uncompressed,
          threads):
    """
    Computes a database of phylo-k-mers.
    """
    build_database(ar,
                   refalign, reftree, states,
                   verbosity,
                   workdir, write_reduction, #dbfilename,
                   alpha, categories,
                   k, model, convert_uo, #gap_jump_thresh,
                   no_reduction, reduction_ratio, omega,
                   filter, mu, ghosts, use_unrooted, merge_branches,
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
                   alpha, categories,
                   k, model, convert_uo, #gap_jump_thresh,
                   no_reduction, reduction_ratio, omega,
                   filter, mu, ghosts, 
                   use_unrooted, merge_branches,
                   ar_dir, ar_only, ar_config,
                   keep_positions, uncompressed,
                   threads):

    if not ar:
        ar = find_raxmlng()

    # create working directory
    Path(workdir).mkdir(parents=True, exist_ok=True)
    current_dir = os.path.dirname(os.path.realpath(__file__))

    ar_parameters = parse_config(ar_config) if ar_config else ""

    # If IPK is installed, look for the binary in the installed location,
    # otherwise it is run from sources
    ipk_bin_dir = f"{current_dir}" if os.path.exists(f"{current_dir}/ipk-dna") else f"{current_dir}/bin/ipk"
    
    # run IPK
    if states == 'nucl':
        if keep_positions:
            raise RuntimeError("--keep-positions is not supported for DNA.")
        else:
            bin = f"{ipk_bin_dir}/ipk-dna"
    else:
        if keep_positions:
            bin = f"{ipk_bin_dir}/ipk-aa-pos"
        else:
            bin = f"{ipk_bin_dir}/ipk-aa"


    command = [
        bin,
        "--ar-binary", str(ar),
        "--refalign", str(refalign),
        "-t", str(reftree),
        "-w", str(workdir),
        "-k", str(k),
        "--model", model,
        "--alpha", str(alpha),
        "--categories", str(categories),
        "--reduction-ratio", str(reduction_ratio),
        "--omega", str(omega),
        "--" + filter.lower(),
        "--" + ghosts.lower(),
        "-u", str(mu),
        "-j", str(threads)
    ]

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
    print("Running", command_str)
    print()
    try:
        subprocess.run(command_str, shell=True, check=True)

        # clean after
        subprocess.call(["rm", "-rf", hashmaps_dir])

    except subprocess.CalledProcessError as error:
        print(error.output)


if __name__ == "__main__":
    ipk()
