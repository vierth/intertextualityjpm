'''
This script takes each file in the results directory and compiles them
into a single file. It also filters out common results (and results that
are highly similar to the common quotes) if the user chooses to do so.
'''

import os
import sys
import time

import Levenshtein

from config.config import read_config_file

#**********************#
# FUNCTION DEFINITIONS #
#**********************#

# This function removes common quotes, and similar quotes as requested


def remove_common(quoteinfo, params):
    common_repetition_threshold = params["common_repetition_min"]
    short_quote_length = params["short_quote_length"]
    filter_similar = params["filter_similar"]
    similarity_threshold = params["similarity_threshold"]
    limit_similarty_check = params["limit_similarity_check"]
    limit_extent = params["limit_extent"]
    # Identify common quotes
    print(f"\nCounting short quotes")
    quote_scores = {}
    for line in quoteinfo:
        info = line.split("\t")
        # remove alignment information
        relevant_quote = info[6]
        # Calculate short quote frequency
        if len(relevant_quote) < short_quote_length:
            try:
                quote_scores[relevant_quote] += 1
            except:
                quote_scores[relevant_quote] = 1

    print(f"Identifying high incidence quotes")
    # return quotes that should be discarded
    to_discard = set([q for q, v in quote_scores.items()
                     if v > common_repetition_threshold])

    if filter_similar:
        print("Removing high incidence and similar results")
    else:
        print("Removing high incidence quotes")

    temptime = time.time()
    totaltime = 0
    save = []
    total_quotes = len(quoteinfo)

    # Iterate through the information and discard unwanted information
    for line in quoteinfo:
        add = True
        info = line.split("\t")
        relevant_quote = info[6]
        # Save the quote if it is longer than the cutoff
        if relevant_quote in to_discard:
            add = False
        # Otherwise, check its length. If it is shorter than the cuttoff
        # and filtersimilar is set to true, check to see if it is similar
        # to discarded quotes quotes.
        elif len(relevant_quote) < short_quote_length and filter_similar:
            for common in to_discard:
                if limit_similarty_check:
                    if common[0:limit_extent] == relevant_quote[0:limit_extent]:
                        if Levenshtein.ratio(relevant_quote, common) >= similarity_threshold:
                            add = False
                            break
                else:
                    if common[0] == relevant_quote[0]:
                        if Levenshtein.ratio(relevant_quote, common) >= similarity_threshold:
                            add = False
                            break
        if add:
            save.append(line)
        total_quotes -= 1
        if total_quotes % 10000 == 0:
            ft = time.time() - temptime
            totaltime += ft
            sys.stdout.write(
                f"{total_quotes} quotes remaining. {ft:.2f}/{totaltime:.2f} secs                \r")
            sys.stdout.flush()
            temptime = time.time()

    quoteinfo = save
    return quoteinfo

#*********************#
# START OF MAIN LOGIC #
#*********************#


def run(cfg):
    # Start timer
    gs = time.time()

    # extract paramaters
    params = cfg["params"]

    # Container for the resulting information
    results = []

    # Compile results into a single list
    for root, dirs, files in os.walk(cfg["paths"]["quote_result_corpus_folder"]):
        for i, f in enumerate(files):
            with open(f"{root}/{f}", "r", encoding="utf8") as rf:
                contents = rf.read().split("\n")
                for line in contents:
                    if line[:11] != "TargetTitle" and len(line) > 0:
                        results.append(f"{f[:-4]}\t{line}")
            sys.stdout.write(f"{i + 1} results of {len(files)} processed.\r")
            sys.stdout.flush()

    # Filter empty results (just in case)
    results = [r for r in results if len(r) > 0]

    
    # Filter out common quotes if desired
    if params["filter_common"]:
        results = remove_common(results, params)

    # Write results to file
    with open(cfg["files"]["output_file"], "w", encoding="utf8") as wf:
        wf.write(
            "SourceTitle\tTargetTitle\tLength\tRatio\tSourcePlace\tTargetPlace\tSourceText\tTargetText\n")
        wf.write("\n".join(results))

    ge = time.time()
    gt = ge-gs
    print(f"Global Operation completed in {gt:.2f} seconds")


if __name__ == "__main__":
    # read config file
    cfg = read_config_file("config/compile_config.json")
    run(cfg)
