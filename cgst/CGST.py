#!/usr/bin/env python3
"""
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/Nilad/CGST.git

This script is the main step of CGST.

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
import os
import argparse
import logging
import textwrap

from cgst.helpFormatter import MyParser, MyHelpFormatter

from cgst.build_database import build_database
from cgst.version import __version__
from cgst.analysis import main_analysis
from cgst.log import section_header, explanation, log
from cgst.utils import mentalist_path_and_version, available_species_mentalist
from cgst.detection import main_detection


def check_for_required_tools():
    """
    The function that check the dependencies
    """
    section_header('Checking requirements')
    explanation('CGST requires MentaLiST to run, so it checks for this tool now.')
    mentalist_path, mentalist_version, mentalist_status = mentalist_path_and_version('mentalist')
    if mentalist_status == 'good':
        log(f'MentaLiST found: {mentalist_path} (v{mentalist_version})')
    elif mentalist_status == 'not found':
        os.sys.exit('Error: MentaLiST not found - make sure it is in your PATH before running '
                    'CGST')
    elif mentalist_status == 'bad':
        os.sys.exit('Error: unable to determine MentaLiST version')
    log()


def status_genome(database_path):
    """
    The function that display the genome installed
    :param database_path: The CGST database directory path
    """
    for file in os.listdir(database_path):
        log(f"Type of database : {file}\n")
        file_path = os.path.join(database_path, file)
        if os.path.isdir(file_path):
            log("Species:")
            for file2 in os.listdir(file_path):
                log(f" - {file2}")
        log()


def pre_main(args):
    """
    The pre-main function that receives arguments
    :param args: the arguments
    """
    force = output = species = ""
    threads = 0

    if args.force:
        force = os.path.abspath(args.force)
    if args.which == "check":
        check_for_required_tools()
    elif args.which == "detection" or args.which == "analysis":
        database = work_dir = ""
        if args.threads:
            threads = args.threads
        if args.database:
            database_path_abs = os.path.abspath(args.database[0])
            if not os.path.exists(database_path_abs):
                log(f"The Database folder {args.database[0]} do not exist\n")
                os.sys.exit()
            else:
                database = database_path_abs
        if args.work_dir:
            work_dir = os.path.abspath(args.work_dir[0])
        if args.species:
            species = args.species
        else:
            log("The Species is not provided\n")
            os.sys.exit()
        if args.which == "analysis":
            detection_dir = ""
            if args.detectionDir:
                if not os.path.exists(args.detectionDir[0]):
                    log(f"The detection directory {args.detectionDir[0]} do not exist\n")
                    os.sys.exit()
                else:
                    detection_dir = os.path.abspath(args.detectionDir[0])
            main_analysis(detection_dir, database, work_dir, species, force, threads)
        if args.which == "detection":
            r1 = r2 = ""
            only_mentalist = False
            if args.output:
                output = args.output[0]
            else:
                log("A name for the process is required\n")
                os.sys.exit()

            if args.only_mentalist:
                only_mentalist = True
            if args.R1:
                if not os.path.exists(args.R1[0]):
                    log(f"The Reads file {args.R1[0]} do not exist\n")
                    os.sys.exit()
                else:
                    r1 = os.path.abspath(args.R1[0])
            if args.R2:
                if not os.path.exists(args.R2[0]):
                    log(f"The Reads file {args.R2[0]} do not exist\n")
                    os.sys.exit()
                else:
                    r2 = os.path.abspath(args.R2[0])
            main_detection(r1, r2, database, work_dir, output, species, force, only_mentalist, threads)
    elif args.which == "build-database":
        cgst_database_path = ""
        if args.database:
            cgst_database_path = os.path.abspath(args.database[0])
        else:
            log("The Database is not provided\n")
            os.sys.exit()
        if args.species:
            species = args.species
        else:
            log("The Species is not provided\n")
            os.sys.exit()
        build_database(cgst_database_path, species, threads)
    elif args.which == "available-species":
        available_species_mentalist("")
    elif args.which == "status-genome":
        cgst_database_path = ""
        if args.database:
            cgst_database_path = os.path.abspath(args.database[0])
        status_genome(cgst_database_path)
    else:
        logging.error("The command is not provided. Please follow usage.\n")
        os.sys.exit()


def run():
    parser = MyParser(
        prog='CGST',
        usage="CGST [commands][options] ",
        description=f'CGST: Core-Genome Sequence Typing - Version {__version__}',

        add_help=False, formatter_class=MyHelpFormatter,
        epilog=textwrap.dedent('''\
            Tool available on GitHub https://github.com/CNRResistanceAntibiotic/CGST
                --Created by Aurélien Birer--
                abirer36@gmail.com
                --CNR "Résistance aux Antibiotiques" - Clermont-Ferrand--
                2020
         ''')
    )

    subparsers = parser.add_subparsers(help='Description', metavar="<commands>")

    # ALL Commands
    check_parser = subparsers.add_parser('check', help='check Dependencies')
    available_species_parser = subparsers.add_parser('available-species', help='Available Species in cgmlst.org')
    build_database_parser = subparsers.add_parser('build-database', help='Build Database by species')
    status_genome_parser = subparsers.add_parser('status-genome', help='Resume of genomes available')
    detection_parser = subparsers.add_parser('detection', help='Detection Sequence Typing')
    analysis_parser = subparsers.add_parser('analysis', help='Analysis strains that had a cgMLST detection')

    # Check Args
    check_parser.set_defaults(which='check')

    # Check Analysis
    analysis_parser.set_defaults(which='analysis')
    analysis_parser = analysis_parser.add_argument_group('Positional arguments')
    analysis_parser.add_argument('-dd', '--detect_dir', dest="detectionDir", default='', nargs=1,
                                 help='Directory with all detection files "final" provided by the detection',
                                 required=True)
    analysis_parser.add_argument('-db', '--db_path', dest="database", default='', nargs=1,
                                 help='CGST database Path', required=True)
    analysis_parser.add_argument('-w', '--work_dir_path', dest="work_dir", default='', nargs=1,
                                 help='Working Directory Path', required=True)
    analysis_parser.add_argument('-sp', '--species', dest="species", default='',
                                 help='Quoted scientific name of species wanted like "Escherichia coli"',
                                 required=True)
    analysis_parser.add_argument('-t', '--threads', dest='threads', default='4',
                                 help="Threads define (Default 4)")

    # Available Species
    available_species_parser.set_defaults(which='available-species')

    # Build Database
    build_database_parser.set_defaults(which='build-database')
    build_database_parser = build_database_parser.add_argument_group('Positional arguments')
    build_database_parser.add_argument('-db', '--db_path', dest="database", default='', nargs=1,
                                       help='CGST database Path', required=True)
    build_database_parser.add_argument('-sp', '--species', dest="species", default='',
                                       help='Quoted scientific name of species wanted like "Escherichia coli"',
                                       required=True)
    build_database_parser.add_argument('-t', '--threads', dest='threads', default='4',
                                       help="Threads define (Default 4)")

    # Detection Args
    detection_parser.set_defaults(which='detection')
    detection_parser = detection_parser.add_argument_group('Positional arguments')
    detection_parser.add_argument('-1', dest="R1", default='', nargs=1, help='R1 fastq file', required=True)
    detection_parser.add_argument('-2', dest="R2", default='', nargs=1, help='R2 fastq file', required=True)
    detection_parser.add_argument('-db', '--db_path', dest="database", default='', nargs=1, help='CGST database Path',
                                  required=True)
    detection_parser.add_argument('-w', '--work_dir_path', dest="work_dir", default='', nargs=1,
                                  help='Working Directory Path', required=True)
    detection_parser.add_argument('-o', '--output', dest="output", default='', nargs=1,
                                  help='Output Name of file. If not provided, take name of the R1 fastq file',
                                  required=True)
    detection_parser.add_argument('-sp', '--species', dest="species", default='',
                                  help='Quoted scientific name of species wanted like "Escherichia coli"',
                                  required=True)
    detection_parser.add_argument('-mentalist', '--only-mentalist', dest="only_mentalist", default=False,
                                  help='Just run Mentalist in detection', required=False)
    detection_parser.add_argument('-t', '--threads', dest='threads', default='4', help="Threads define (Default 4)")

    # Check Status
    status_genome_parser.set_defaults(which='status-genome')

    status_genome_parser.add_argument('-db', '--db_path', dest="database", default='', nargs=1,
                                      help='CGST database Path', required=True)

    # Global args
    parser.add_argument('-f', '--force', dest="force", default='False', nargs='?', help="Overwrite output directory")
    parser.add_argument('-V', '--version', action='version', version=f'version-{__version__}',
                        help="Prints version number")
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                        help='Show this help message and exit')

    args = parser.parse_args()
    if len(os.sys.argv) == 1:
        parser.print_help()
        os.sys.exit(1)
    pre_main(args)


if __name__ == '__main__':
    run()
