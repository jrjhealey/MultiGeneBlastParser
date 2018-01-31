"""
Extract and parse data from the text output of
MultiGeneBlast. This script will produce a table of
starts and ends, that can be used with the genbank slicer
script to subset operons for each resulting hit.
"""

import argparse
import io
import os
import re
import sys
import traceback

import pandas as pd

pd.set_option('expand_frame_repr', False)

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "MGBparser"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"



def get_args():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser(
            description='Parse and extract data from the output of MultiGeneBlast.')
        parser.add_argument(
            '-v',
            '--verbose',
            action='store_true',
            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument(
            '-r',
            '--references',
            action='store_true',
            help='Display relevant references. This option cannot be used with any others.')
        parser.add_argument(
            '-q',
            '--query',
            action='store_true',
            help='Output the query sequence information (and write a file).')
        parser.add_argument(
            '-s',
            '--sighits',
            action='store_true',
            help='Output the table of significant hits (and write a file).')
        parser.add_argument(
            '-b',
            '--blastfile',
            action='store_true',
            help='Output the table of BLAST hits for each match (and write a file).')
        parser.add_argument(
            '-o',
            '--outfile',
            action='store',
            help='Name stem for all output files (extensions are set internally).')
        parser.add_argument(
            '-m',
            '--max_result',
            type=int,
            default=50,
            action='store',
            help='Return results for just the top n details sections (Def: 50).')
        parser.add_argument(
            '--clusterfile',
            action='store',
            default='./clusterblast_output.txt',
            help='The text file of hits output by MGB. By default this is called \'clusterblast_output.txt\'.')

    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()


def display_refs():
    """Display relevant references"""
    print('''
    - MultiGeneBlast has been previously published here:
        Medema MH, Takano E, Breitling R. 
        "Detecting Sequence Homology at the Gene Cluster Level with MultiGeneBlast."
         Molecular Biology and Evolution. 2013;30(5):1218-1223. doi:10.1093/molbev/mst025.
    ''')
    sys.exit(0)


def create_header_class(header_section):
    """
    Clean up and parse the header section in to query-able formats
    """

    class Subject(object):
        """
        This class captures the attributes of the query sequence present at the top of the MGB output file.
        This is a separate class, since this information occurs once within a file, and would be inefficient to
        re-parse as part of the MGB_hit class for every new hit.

        Attributes:
            Filename: The input file provided to MGB originally (line 1 of output)
            Table: The tab separated list of features
            Columns: The column names from the section passed to
        """

        def __init__(self, filename, table, columns):
            self.filename = filename
            self.table = table
            self.columns = columns

    # Extract filename from line 1 of the file (this is a bit brute force, might break since it disregards
    # the string at the start of the line).
    Subject.filename = os.path.basename(header_section[0])

    # Define the column headers for the section since the file's are too verbose and ambiguous
    Subject.columns = ["Locus", "Start", "Stop", "Strand", "Annotation", "Comment"]

    # Store the table of loci and associated data (tab separated, removing last blank column.
    # NB > Not sure if the right hand column is always blank. Check with additional files.

    # Use StringIO object to imitate a file, which means that we can use read_table and have the dtypes
    # assigned automatically (necessary for functions like min() to work correctly on integers)
    Subject.table = pd.read_table(io.StringIO(u'\n'.join([row.rstrip('\t') for row in header_section[2:]])),
                                  names=Subject.columns)

    return Subject


def create_sighits_class(sighits_section):
    """
    Clean up and parse the significant hits section into a query-able class.
    """
    class SigHit(object):
        """
        Class to hold the table of significant hits.

        Attributes:
            Table: Table of significant hits in rank order
        """

        def __init__(self, Table):
            self.Table = Table

    # Define the column headers for the section since the file's are too verbose and ambiguous
    SigHit.Columns = ["Rank", "ID", "Description"]

    # Store the table of loci and associated data (tab separated, removing last blank column.

    # Use StringIO object to imitate a file, which means that we can use read_table and have the dtypes
    # assigned automatically (necessary for functions like min() to work correctly on integers)
    SigHit.Table = pd.read_table(io.StringIO(u'\n'.join([row.rstrip('.') for row in sighits_section])),
                                 sep='\.|\t', engine='python', names=SigHit.Columns)

    return SigHit


def create_hit_class(hit_sublist):
    """
    Clean up and parse the Detailed hits section in to a class. This class is the meat of this script.
    """
    class Hit(object):
        """
        Store each MGB hit as a class object so as to group all the attributes.

        Attributes:
            hit_no: The rank number of the hit returned from MultiGeneBlast.
            hit_id: The name assigned to the hit rank.
            source: The extended description of the hit/its sequence origin.
            protein_no: The number of proteins with hits in the detected cluster.
            MGB_score: The weighted MGB score used to rank the hits with synteny etc.
            cubit_score: The cumulative bit-score of all the BLAST hits within the cluster.
            location_table: The table of gene locations as a pandas dataframe
            location_columns: The names of the columns in the pandas dataframe (set internally)
            blast_table: The table of BLAST gene hits, as a pandas dataframe
            blast_columns: The names of the columns in the pandas dataframe (set internally)
            operon_start: Beginning base index for the match gene cluster
            operon_end: End base index for the matched gene cluster
        """

        def __init__(self, hit_no, hit_id, source, protein_no, MGB_score, cubit_score,
                           location_table, location_columns, blast_table, blast_columns,
                           operon_start, operon_end, operon_length,
                           dominant_strand, start_locus, end_locus):
            """Initialise a MGB hit object"""

            self.hit_no = hit_no
            self.hit_id = hit_id
            self.source = source
            self.protein_no = protein_no
            self.MGB_score = MGB_score
            self.cubit_score = cubit_score
            self.location_table = location_table
            self.location_columns = location_columns
            self.blast_table = blast_table
            self.blast_columns = blast_columns
            self.operon_length = operon_start
            self.operon_length = operon_end
            self.operon_length = operon_length
            self.dominant_strand = dominant_strand
            self.start_locus = start_locus
            self.end_locus = end_locus

    # Collect hit 'metadata'
    Hit.hit_no, Hit.hit_id = [a.lstrip(' ') for a in hit_sublist[0].split('.')]
    Hit.source = hit_sublist[1].replace('Source: ', '').rstrip('.')
    Hit.protein_no = re.match('.*?([0-9]+)$', hit_sublist[2]).group(1)
    Hit.MGB_score = re.findall("\d+\.?\d*", hit_sublist[3])
    Hit.cubit_score = re.findall("\d+\.?\d*", hit_sublist[4])

    # Collect hit gene locations
    Hit.location_columns = ["Locus", "Start", "Stop", "Strand", "Annotation", "Comment"]

    # Find end of location table:
    table_split = hit_sublist.index(
        'Table of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):')
    end_location_table = table_split
    start_blast_table = table_split + 1

    Hit.location_table = pd.read_table(
        io.StringIO(u'\n'.join([row.rstrip('\t') for row in hit_sublist[6:end_location_table]])),
                    names=Hit.location_columns)

    # Figure out which is the most common strand
    Hit.dominant_strand = Hit.location_table["Strand"].value_counts().idxmax()
    # Switch between strands
    if Hit.dominant_strand == '+':
        Hit.operon_start = min(Hit.location_table['Start'])
        Hit.start_locus = Hit.location_table['Locus'][Hit.location_table['Start'].idxmin()]

        Hit.operon_end = max(Hit.location_table['Stop'])
        Hit.end_locus = Hit.location_table['Locus'][Hit.location_table['Start'].idxmax()]

        Hit.operon_length = Hit.operon_end - Hit.operon_start

    elif Hit.dominant_strand == '-':
        Hit.operon_start = max(Hit.location_table['Start'])
        Hit.start_locus = Hit.location_table['Locus'][Hit.location_table['Start'].idxmax()]

        Hit.operon_end = min(Hit.location_table['Stop'])
        Hit.end_locus = Hit.location_table['Locus'][Hit.location_table['Start'].idxmin()]

        Hit.operon_length = Hit.operon_start - Hit.operon_end
    # Collect hit BLAST results

    Hit.blast_columns = ["Query", "Subject", "PercID", "Score", "PercCoverage", "E-Value"]
    Hit.blast_table = pd.read_table(
        io.StringIO(u'\n'.join([row.rstrip('\t') for row in hit_sublist[start_blast_table:]])),
                    names=Hit.blast_columns)

    return Hit


def parse_section(file, delim1, delim2):
    """Separate files in to sections according to delimiter pairs"""

    regex = '{}(.*?){}'.format(delim1, delim2)
    for result in re.findall(regex, file, re.S):
        # Throw away any blank lines remaining
        result = filter(None, result.split('\n'))

        return result


def main():
    """Call functions and parse results of MGB."""

    args = get_args()
    if args.references:
        display_refs()

    if args.verbose: print("Opening " + args.clusterfile + " for reading...")
    with open(args.clusterfile, 'r') as cfh:
        content = cfh.read()
        if args.verbose:
            print("Parsing Query details section...")
        header_section = parse_section(content, "^", "Significant hits:")
        if args.verbose:
            print("Parsing Significant hit section...")
        sighits_section = parse_section(content, "Significant hits:", "Details:")
        if args.verbose:
            print("Parsing Hit details section...")
        hit_section = parse_section(content, "Details:", "\Z")

    # Pass to simple parsing functions as each of these sections appears once
    Subject = create_header_class(header_section)
    SigHits = create_sighits_class(sighits_section)

    # Detailed hit section is more complicated as the section needs to
    # be broken up in to separate hits...
    if args.verbose:
        print("Splitting Hit details section...")
    all_hit_lists = []
    sublist = []
    for line in hit_section:
        if line == '>>':
            sublist = []
            all_hit_lists.append(sublist)
        else:
            sublist.append(line)

    # Create classes of each of the separated hits
    if args.verbose:
        print("Establishing Hit classes...")
    hit_classlist = []
    for entry in all_hit_lists:
        hit_classlist.append(create_hit_class(entry))

    # Write output files
    if args.outfile is None:
        args.outfile = args.clusterfile

    # Set up outfile filepaths
    query_outfile = os.path.join(os.path.dirname(args.clusterfile), args.outfile) + '_queryinfo.tsv'
    sighit_outfile = os.path.join(os.path.dirname(args.clusterfile), args.outfile) + '_sighitinfo.tsv'

    # Prepare header info
    if args.query is True:
        if args.verbose:
            print("Writing query details to file: " + query_outfile)
        print("Query sequence information:")
        print("===========================")
        print("Input file:" + Subject.filename)
        print(Subject.table)
        with open(query_outfile, 'w') as qfh:
            qfh.write("Details of query sequence:" + Subject.filename)
            Subject.table.to_csv(qfh, sep='\t')

    # Prepare Sighits info
    if args.sighits is True:
        if args.verbose:
            print("Writing Significant hit details to file: " + sighit_outfile)
        print("Significant Hit information:")
        print("============================")
        print(SigHits.Table.iloc[0:args.max_result])
        with open(sighit_outfile, 'w') as sfh:
            SigHits.Table.iloc[0:args.max_result].to_csv(sfh, sep='\t')

    # If fewer results than specified max, modify max for enumerate to display correctly
    if len(hit_classlist) < args.max_result:
        args.max_result = len(hit_classlist)

    # Prepare Details info on a class-by-class basis
    for i, Hit_instance in enumerate(hit_classlist[0:args.max_result]):
        location_outfile = os.path.join(
            os.path.dirname(args.clusterfile),
            args.outfile) + '_' + Hit_instance.hit_id + '_locationinfo.tsv'
        blast_outfile = os.path.join(
            os.path.dirname(args.clusterfile),
            args.outfile) + '_' + Hit_instance.hit_id + '_blastinfo.tsv'
        coords_outfile = os.path.join(
            os.path.dirname(args.clusterfile),
            args.outfile) + '_' + Hit_instance.hit_id + '_coords.tsv'

        if args.verbose:
            print("Writing Hit details for: {0}. {1} to {2} ({3} of {4})".format(
                Hit_instance.hit_no, Hit_instance.hit_id, location_outfile, i+1, args.max_result))
            print("Hit location information:")
            print("============================")
            print(Hit_instance.location_table)
        with open(location_outfile, 'w') as lfh:
            Hit_instance.location_table.to_csv(lfh, sep='\t')

        if args.blastfile is True:
            if args.verbose:
                print("Writing Hit BLAST details for: {0}. {1} to {2} ({3} of {4})".format(
                    Hit_instance.hit_no, Hit_instance.hit_id, blast_outfile, i+1, args.max_result))
                print("Hit location information:")
                print("============================")
                print(Hit_instance.location_table)
            with open(blast_outfile, 'w') as bfh:
                Hit_instance.blast_table.to_csv(bfh, sep='\t')

        coordlist = [Hit_instance.hit_no,
                     Hit_instance.hit_id,
                     Hit_instance.start_locus,
                     Hit_instance.end_locus,
                     Hit_instance.operon_start,
                     Hit_instance.operon_end,
                     Hit_instance.dominant_strand,
                     Hit_instance.source]
        coordstring = '\t'.join(map(str, coordlist))

        if args.verbose:
            print("Writing Hit coordinates for: {0}. {1} to {2} ({3} of {4})".format(
                Hit_instance.hit_no, Hit_instance.hit_id, coords_outfile, i+1, args.max_result))
            print("Hit coordinate information:")
            print("===========================")
            print('\t'.join(["Hit No", "ID",
                            "Start Locus",
                            "End Locus",
                            "Start Index",
                            "End Index",
                            "Main Strand",
                            "Source"]))
            print(coordstring)

        with open(coords_outfile, 'w') as cfh:
            cfh.write(coordstring)


if __name__ == "__main__":
    main()
