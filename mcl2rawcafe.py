import argparse
import sys
import os
import re
import copy

def mcl2rawcafe(mcl_dump, path_output_file, species_ids):
    """
    mcl_dump: mcl's dump file
    path_output_file: path to output .txt file
    species_ids: list of unique species identifiers (separated by white spaces, enclosed within double quotes)
                 (e.g., "ENSFCA ENSCAF ENSAME ENSMPU ENSBTA ENSSSC ENSECA ENSP00 ENSLAF")
    """
    sps_header = "\t".join(sp_id for sp_id in species_ids)
    species_ids_dict_0 = dict((sp_id, 0) for sp_id in species_ids) # initializing dictionary template

    with open(path_output_file, "w") as output_file:
        output_file.write("Desc\tFamily ID\t" + sps_header + "\n") # writing CAFE input header

        with open(mcl_dump, "r") as input_file:
            line_count = 1 # using a simple counter as Family ID's

            for line in input_file:
                species_gene_ids = line.rstrip().split() # list of species id in this line
                species_ids_dict = copy.deepcopy(species_ids_dict_0) # dictionary for this line

                for sp_gene_ids in species_gene_ids:
                    for sp_id in species_ids:
                        if sp_gene_ids.find(sp_id) != -1:
                            try:
                                species_ids_dict[sp_id] += 1

                            except:
                                species_ids_dict[sp_id] = 1

                output_file.write("\t".join(["(null)", str(line_count), ""])) # first 2 cols + "\t"

                for sp_id in species_ids:
                    output_file.write(str(species_ids_dict[sp_id]) + "\t") # gene family counts for all species

                output_file.write("\n")

                line_count += 1

if __name__ == "__main__":
    """
    Running standalone
    """
    parser = argparse.ArgumentParser(prog="mcl2rawcafe.py", description="Parses mcl's output dump file into a raw table CAFE can read")
    parser.add_argument("-i", "--input-file", action="store", dest="input_file", default=None, type=str, help="full path to mcl's output dump file")
    parser.add_argument("-o", "--output-file", action="store", dest="output_file", default=None, type=str, help="full path to file to be written")
    parser.add_argument("-sp", "--species-ids", action="store", dest="sps", default=None, type=str, help="list of unique species identifiers (separated by white spaces) present in the input file")

    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        exit("Could not find input file (mcl's output dump file). Exiting...\n")

    if len(args.sps.split(" ")) <= 1:
        exit("You specified a single species identifier. Exiting...\n")

    mcl2rawcafe(args.input_file, args.output_file, args.sps.split(" "))
