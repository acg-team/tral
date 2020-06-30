"""TRAL configuration tool
"""


import argparse
import logging
import shutil
import tarfile
import os
from . import paths


def install(args):
    "Install tral config from package"
    tralhome = args.path

    if os.path.exists(tralhome):
        logging.warn("The TRAL configuration directory {} already exists.".format(tralhome))
    else:
        pkg_config = os.path.join(paths.PACKAGE_DIRECTORY, "tral_configuration")
        shutil.copytree(pkg_config, tralhome)
        logging.info("Installing tral configuration to {}".format(tralhome))


def download(args):
    "Install tral configs, plus download the full data files"
    # Minimum file size; used to detect download errors
    MIN_FILE_SIZE = 0

    install(args)

    # Download
    expanded = os.path.join(args.path, "data", "pvalue")
    if os.path.isdir(expanded):
        logging.warn("TRAL data already found at {}".format(expanded))
        return

    url = paths.CONFIG_DATA_URL + "pvalue.tar.gz"
    compressed = os.path.join(args.path, "data", "pvalue.tar.gz")
    if not os.path.isfile(compressed) or os.path.getsize(compressed) <= MIN_FILE_SIZE:
        logging.info("Downloading pvalue files (2.7 GB). This may take some time, depending on your connection.")
        paths.download(url, compressed)

    logging.info("Expanding pvalue files")
    with tarfile.open(compressed, "r:gz") as tar:
        tar.extractall(os.path.join(args.path, "data"))
        assert os.path.isdir(expanded)

    # Clean up
    os.remove(compressed)


def main(args=None):
    parser = argparse.ArgumentParser(description='TRAL configuration tool', prog="tral")
    parser.add_argument("-v", "--verbose", help="Verbose logging",
                        dest="verbose", default=logging.INFO,
                        action="store_const", const=logging.DEBUG)
    parser.add_argument("-q", "--quiet", help="Quiet logging",
                        dest="verbose", action="store_const",
                        const=logging.ERROR)
    subparsers = parser.add_subparsers(title="subcommands", description="valid subcommands")

    install_parser = subparsers.add_parser("install", help="Install Basic TRAL configs")
    install_parser.add_argument("path", type=str,
                                nargs='?',
                                default=paths.CONFIG_DIR,
                                help="TRAL configuration directory path")
    install_parser.set_defaults(func=install)

    install_parser = subparsers.add_parser("download", help="Download complete TRAL data")
    install_parser.add_argument("path", type=str,
                                nargs='?',
                                default=paths.CONFIG_DIR,
                                help="TRAL configuration directory path")
    install_parser.set_defaults(func=download)

    parsed_args = parser.parse_args(args)

    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=parsed_args.verbose)

    if hasattr(parsed_args, 'func'):
        parsed_args.func(parsed_args)
    else:
        parser.error("No command given")


if __name__ == "__main__":
    main()
