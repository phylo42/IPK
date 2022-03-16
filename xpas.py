#!/usr/bin/env python3
"""
xpas wrapper script.
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
from os import path
from pathlib import Path
import click
import subprocess
from pathlib import Path
from enum import Enum


@click.group()
def xpas():
    """
    xpas

    N. Romashchenko, B. Linard, F. Pardi, E. Rivals
    """
    pass


NUCL_MODELS = ['JC69', 'HKY85', 'K80', 'F81', 'TN93', 'GTR']
AMINO_MODELS = ['LG', 'WAG', 'JTT', 'Dayhoff', 'DCMut', 'CpREV', 'mMtREV', 'MtMam', 'MtArt']
ALL_MODELS = NUCL_MODELS + AMINO_MODELS


KMER_FILTERS = ["no-filter", "mif0", "mif1", "random"]


class ScoreModel(Enum):
    MAX = 1
    EXISTS = 2

SCORE_MODELS = [e.name.lower() for e in ScoreModel]


def validate_filter(ctx, param, value):
    value = value.lower()
    if value not in KMER_FILTERS:
        valid_values = ', '.join(v for v in KMER_FILTERS)
        raise click.BadParameter('Filter must be one of: ' + valid_values)
    return value


def validate_score_model(ctx, param, value):
    value = value.lower()
    if value not in SCORE_MODELS:
        valid_values = ', '.join(v for v in SCORE_MODELS)
        raise click.BadParameter('Filter must be one of: ' + valid_values)
    return value


@xpas.command()
@click.option('-b', '--arbinary',
              type=click.Path(exists=True),
              required=True,
              help="Binary file for marginal AR, currently 'phyml' and "
                   "'baseml' (from PAML) are supported.")
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
              help=" Write reduced alignment to file.")
#@click.option('--dbfilename',
#              type=str,
#              default="db.rps", show_default=True, 
#              help="Output database filename.")
@click.option('-a', '--alpha',
              type=float,
              default=1.0, show_default=True,
              help="Gammma shape parameter used in ancestral reconstruction.")
@click.option('-c', '--categories',
              type=int,
              default=4, show_default=True,
              help="Number of categories used in ancestral reconstruction.")
#@click.option('-g', '--ghosts',
#              type=int,
#              default=1, show_default=True,
#              help="Number of ghost nodes injected per branch.")
@click.option('-k', '--k',
             type=int,
             default=8, show_default=True,
             help="k-mer length used at DB build.")
@click.option('-m', '--model',
             type=click.Choice(ALL_MODELS),
             required=True,
             help="Model used in AR, one of the following:\n"
                  f"nucl: {', '.join(x for x in NUCL_MODELS)}\n"
                  f"amino: {', '.join(x for x in AMINO_MODELS)}")
@click.option('--arparameters',
             type=str,
             help="""Parameters passed to the software used for
                  anc. seq. reconstuct. Overrides -a,-c,-m options.
                  Value must be quoted by ' or ". Do not set options
                  -i,-u,--ancestral (managed by RAPPAS).""")
@click.option('--convert-uo',
              is_flag=True,
              help="U, O amino acids are converted to C, L.")
#@click.option('--gap-jump-thresh',
#              type=float,
#              deafult=0.3, show_default=True,
#              help="Gap ratio above which gap jumps are activated.")
@click.option('--no-reduction',
              is_flag=True,
              help="""Do not operate alignment reduction. This will 
                  keep all sites of input reference alignment and 
                  may produce erroneous ancestral k-mers.""")
@click.option('--ratio-reduction',
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
              type=click.Choice(KMER_FILTERS),
              default="no-filter", show_default=True)
@click.option('--score-model',
              type=click.Choice(SCORE_MODELS),
              default="max", show_default=True)
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
@click.option('--ardir',
             type=click.Path(exists=True, dir_okay=True, file_okay=False),
             help="""Skip ancestral sequence reconstruction, and 
                  uses outputs from the specified directory.""")
@click.option('--aronly',
             is_flag=True,
             default=False, show_default=True,
             help="Dev option. Run only ancestral reconstruction and tree extension. No database will be built")
@click.option('--keep-positions',
              is_flag=True,
              default=False,
              help="""Keeps phylo k-mers positions in the alignment. Makes databases larger in size.""")
@click.option('--threads',
             type=int,
             default=4, show_default=True,
             help="Number of threads used.")
def build(arbinary, #database,
          refalign, reftree, states, verbosity,
          workdir, write_reduction, #dbfilename,
          alpha, categories, #ghosts,
          k, model, arparameters, convert_uo, #gap_jump_thresh,
          no_reduction, ratio_reduction, omega,
          filter, score_model, f, mu, use_unrooted, merge_branches,
          ardir, keep_positions,
          threads, aronly):
    """
    Builds a database of phylo k-mers.

    Minimum usage:

    \tpython xpas.py build -s [nucl|amino] -b ARbinary -w workdir -r alignment.fasta -t tree.newick

    """
    # create working directory
    Path(workdir).mkdir(parents=True, exist_ok=True)
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # auxilary files produced by RAPPAS
    extended_tree = f"{workdir}/extended_trees/extended_tree_withBL.tree"

    # run rappas2
    if states == 'nucl':
        if keep_positions:
            raise RuntimeError("--keep-positions is not supported for DNA.")
        else:
            rappas_bin = f"{current_dir}/bin/build/xpas-build-dna"
    else:
        if keep_positions:
            rappas_bin = f"{current_dir}/bin/build/xpas-build-aa-pos"
        else:
            rappas_bin = f"{current_dir}/bin/build/xpas-build-aa"

    command = [
        rappas_bin,
        "--ar-binary", arbinary,
        "--refalign", str(refalign),
        "-t", str(reftree),
        "-w", str(workdir),
        "-k", str(k),
        "--model", model,
        "--reduction-ratio", str(ratio_reduction),
        "-o", str(omega),
        "--" + filter.lower(),
        "--" + score_model.lower(),
        "-u", str(mu),
        "-j", str(threads)
    ]

    if aronly:
        command.append("--ar-only")
    if ardir:
        command.append("--ar-dir")
        command.append(ardir)
    if no_reduction:
        command.append("--no-reduction")
    if merge_branches:
        command.append("--merge-branches")
    if use_unrooted:
        command.append("--use-unrooted")

    # remove the temporary folder just in case
    hashmaps_dir = f"{workdir}/hashmaps"
    subprocess.call(["rm", "-rf", hashmaps_dir])

    print(" ".join(s for s in command))
    subprocess.call(command)

    # clean after
    subprocess.call(["rm", "-rf", hashmaps_dir])


if __name__ == "__main__":
    xpas()
