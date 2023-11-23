# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Standalone dedicated to taxonomic content (kraken)"""
import argparse
import os
import sys

import colorlog
import rich_click as click
from sequana import KrakenDownload, KrakenPipeline, KrakenSequential
from sequana import sequana_config_path as cfg
from sequana import sequana_config_path as scfg
from sequana.modules_report.kraken import KrakenModule
from sequana.taxonomy import Taxonomy
from sequana.utils import config
from sequana import version as sequana_version


logger = colorlog.getLogger(__name__)


click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.SHOW_ARGUMENTS = True

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

click.rich_click.OPTION_GROUPS = {
    "sequana_taxonomy": [
        {
            "name": "Required",
            "options": ["--input-file1", "--input-file2", "--databases"],
        },
        {
            "name": "Detection Algorithm Options",
            "options": [
                "--confidence",
                "--thread",
            ],
        },
        {
            "name": "Download utilities",
            "options": ["--download", "--update-taxonomy"],
        },
        {
            "name": "Output files",
            "options": ["--show-html", "--no-multiqc", "--output-directory", "--keep-temp-files", 
                    "--unclassified-out", "--classified-out"],
        },
        {
            "name": "Behaviour",
            "options": ["--version", "--level", "--help"],
        },
    ],
}



@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("-1", "--file1", "file1", type=click.Path(), help="""Please use --input-file1 """)
@click.option("-2", "--file2", "file2", type=click.Path(), help="""Please use --input-file2 """)
@click.option("-1", "--input-file1", "file1", type=click.Path(), required=True, help="""R1 fastq file (zipped) """)
@click.option("-2", "--input-file2", "file2", type=click.Path(), help="""R2 fastq file (zipped) """)
@click.option(
    "--databases",
    "databases",
    type=click.Path(),
    nargs="+",
    required=True,
    help="Path to a valid Kraken database(s). If you do not have any, use--download option. You may use several, in which case, an iterative taxonomy is performed as explained in online sequana  documentation (see HierarchicalKRaken"
    "",
)
@click.option(
    "--output-directory",
    "directory",
    type=click.STRING,
    help="""name of the output directory""",
    default="taxonomy",
)
@click.option(
    "--keep-temp-files",
    "keep_temp_files",
    default=False,
    is_flag=True,
    help="keep temporary files (hierarchical case with several databases",
)
@click.option("--thread", "thread", type=click.INT, help="number of threads to use (default 4)", default=4)
@click.option(
    "--show-html",
    "html",
    is_flag=True,
    help="Results are stored in report/ directory and results are not shown by default",
)
@click.option(
    "--download",
    "download",
    default=None,
    type=click.Choice(["toydb"]),
    help="A toydb example to be downloaded.",
)
@click.option("--unclassified-out", help="save unclassified sequences to filename")
@click.option("--classified-out", help="save classified sequences to filename")
@click.option("--confidence", type=click.FLOAT, default=0, show_default=True, help="confidence (kraken2 DB only)")
@click.option(
    "--update-taxonomy",
    is_flag=True,
    show_default=True,
    help="Update the local NCBI taxonomy database to the last version",
)
@click.option("--level", "level", default="INFO", type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]))
@click.version_option(sequana_version)
def main(**kwargs):
    """Welcome to SEQUANA - Taxonomy standalone

    ----

    This standalone tool takes as input one or two (paired-end) FastQ files. 
    Utilizing Kraken, Krona, and some codecs from Sequana, it aims to identify 
    the taxon/organism that matches each read found in the input files.

    Kraken requires an input database. We provide three databases, but you 
    can also use Sequana to help you build a customized one.

    - **toydb**: This database contains only a few measles viruses. Its size
       is only 32Mb and should be used for testing and examples only.
    - **minikraken**: Provided by Kraken's authors, it is about 4Gb and 
      contains viruses and bacteria only.
    - **Sequana-built database**: This 8Gb database is stored on the Synapse 
      website, and you will need an account on synapse.org. It contains 
      approximately 22,000 species, including viruses, bacteria, Homo 
      sapiens, fungi, and more.

    Each database can be downloaded using the following command:

        sequana_taxonomy --download toydb

    After downloading, you can use the tool with a command similar to the following:

        sequana_taxonomy --file1 R1.fastq --file2 R2.fastq
           --databases /home/user/.config/sequana/kraken_toydb
           --show-html --thread 4

    ----

    AUTHORS: Thomas Cokelaer
    Documentation: http://sequana.readthedocs.io
    Issues: http://github.com/sequana/sequana
    """
    from easydev import AttrDict

    options = AttrDict(**kwargs)

    logger.setLevel(options.level)

    if options.update_taxonomy is True:
        tax = Taxonomy()
        logger.info("Will overwrite the local database taxonomy.dat in {}".format(cfg))
        tax.download_taxonomic_file(overwrite=True)
        sys.exit(0)

    if options.download:
        kd = KrakenDownload()
        kd.download(options.download)
        sys.exit()

    fastq = []
    if options.file1:
        if os.path.exists(options.file1) is False:
            logger.error(f"{options.file1} not found")
            sys.exit(1)
        fastq.append(options.file1)
    if options.file2:
        if os.path.exists(options.file2) is False:
            logger.error(f"{options.file1} not found")
            sys.exit(1)
        fastq.append(options.file2)

    if options.databases is None:
        logger.critical("You must provide a database")
        sys.exit(1)

    databases = []
    for database in options.databases:
        guess = os.sep.join([scfg, "kraken2_dbs", database])
        if os.path.exists(guess):  # in Sequana path
            databases.append(guess)
        elif os.path.exists(database):  # local database
            databases.append(database)
        else:
            msg = (
                f"Invalid database name {database}. Neither found locally "
                + f"or in the sequana path {scfg}/kraken2_dbs; Use the --download option"
            )
            logger.error(msg)
            raise IOError

    output_directory = options.directory + os.sep + "kraken"
    os.makedirs(output_directory, exist_ok=True)

    # if there is only one database, use the pipeline else KrakenHierarchical
    _pathto = lambda x: f"{options.directory}/kraken/{x}" if x else x
    if len(databases) == 1:
        logger.info("Using 1 database")
        k = KrakenPipeline(
            fastq,
            databases[0],
            threads=options.thread,
            output_directory=output_directory,
            confidence=options.confidence,
        )

        summary = k.run(
            output_filename_classified=_pathto(options.classified_out),
            output_filename_unclassified=_pathto(options.unclassified_out),
        )
    else:
        logger.info("Using %s databases" % len(databases))
        k = KrakenSequential(
            fastq,
            databases,
            threads=options.thread,
            output_directory=output_directory + os.sep,
            force=True,
            keep_temp_files=options.keep_temp_files,
            output_filename_unclassified=_pathto(options.unclassified_out),
            confidence=options.confidence,
        )
        summary = k.run(output_prefix="kraken")

    with open(output_directory + "/summary.json", "w") as fh:
        import json

        json.dump(summary, fh, indent=4)

    # This statements sets the directory where HTML will be saved
    config.output_dir = options.directory

    # output_directory first argument: the directory where to find the data
    # output_filename is relative to the config.output_dir defined above
    kk = KrakenModule(output_directory, output_filename="summary.html")

    logger.info("Open ./%s/summary.html" % options.directory)
    logger.info("or ./%s/kraken/kraken.html" % options.directory)

    if options.html is True:
        ss.onweb()


if __name__ == "__main__":
    main(sys.argv)
