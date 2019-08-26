"""
Wrapper to Julia's setup_databases function to allow setuptools to install it.
"""

import argparse
import subprocess
import pkg_resources

_JULIA_SCRIPT = pkg_resources.resource_filename('phylosofs',
                                                'src/reconstruct_pir.jl')


def main():
    parser = argparse.ArgumentParser(
        description="""
    This script set up the needed databases for PhyloSofS in the `--output`
    path. If the `--pdb`, `--uniclust` or `--pdb70` argument is used, the
    script is going to create a symbolic link to the indicated folder
    instead of downloading the database.

    It creates a `databases` folder in `--output` containing three folders:
    `pdb`, `uniclust` and `pdb70`.
    """,
        epilog="""
    It has been developed at LCQB (Laboratory of Computational and
    Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=".",
        help="path where the database folder is going to be created.")
    parser.add_argument(
        "--pdb",
        type=str,
        default="",
        help=
        "path to an already existing local folder containing the entire pdb "
        "in mmCIF format.")
    parser.add_argument(
        "--uniclust",
        type=str,
        default="",
        help="path to an already existing local folder containing the uniclust "
        "database from the HH-suite databases.")
    parser.add_argument("--uniclust_version",
                        type=str,
                        default="2018_08",
                        help="Uniclust30 version to be downloaded: YYYY_MM")
    parser.add_argument(
        "--pdb70",
        type=str,
        default="",
        help="path to an already existing local folder containing the "
        "pdb70_from_mmcif database from the HH-suite databases.")
    parser.add_argument("--julia",
                        "-j",
                        type=str,
                        default="julia",
                        help="path to the Julia executable.")

    args = parser.parse_args()

    subprocess.call([
        args.julia, _JULIA_SCRIPT, "--output", args.output, "--pdb", args.pdb,
        "--uniclust", args.uniclust, "--uniclust_version",
        args.uniclust_version, "--pdb70", args.pdb70
    ])


if __name__ == '__main__':
    main()
