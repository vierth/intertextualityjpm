'''
This script finds instances of intertextuality in a large corpus.
This assumes you are using the Anaconda distribution of Python,
otherwise you will have to install a number of dependencies.
The only non-Anaconda library is the Levenshtein library, which
can be installed with pip install python-Levenshtein
'''


import itertools
import json
import os
import pickle
import shutil
import sqlite3
import sys
import time
from itertools import repeat
from multiprocessing import Pool

import Levenshtein

from config.config import read_config_file
from utilities.textanalysis import get_seeds

#*************************#
#  data loading functions #
#*************************#


def load_data(data_file, text_length_file):
    """Loads data from input files and outputs a document with length information
    Titles are at index 0, texts at index 1
    Args:
        data_file (str): the file name containing corpus data created by prepare_corpus
        text_length_file (str): the file name for an output file with info on text lengths

    Returns:
        tuple: returns two lists and a dictionary, one with all titles, one with all texts, and a dictionary relating 
        titles to places in the index
    """
    data = pickle.load(open(data_file, "rb"))
    # Get metadata
    all_titles = data[0]
    # Get texts
    all_texts = data[1]

    # map titles to locations for ease of lookup
    title_to_index = {}
    text_lengths = []
    for i, title in enumerate(all_titles):
        title_to_index[title] = i
        text_lengths.append(f"{title}\t{len(all_texts[i])}")
    with open(text_length_file, "w", encoding="utf8") as wf:
        wf.write("\n".join(text_lengths))

    return all_titles, all_texts, title_to_index


def load_analysis_files(file_name):
    """Loads files with information structured as one item per line and returns a list with the information

    Args:
        file_name (str): path to file

    Returns:
        list: returns data in file or empty list
    """
    if os.path.isfile(file_name):
        with open(file_name, "r", encoding='utf8') as f:
            return f.read().split("\n")
    else:
        return []


def load_index(index_file, total_texts):
    """This function loads a precreated index as created by index_corpus.py

    Args:
        index_file ([type]): [description]
        total_texts ([type]): [description]

    Returns:
        [type]: [description]
    """
    # If a prepared index was provided, extract the data and save into memory
    if not os.path.isfile(index_file):
        print("No index of this name found")
        return False, False, False, False

    print(f"Loading index from {index_file}")
    # Connect to database
    connection = sqlite3.connect(index_file)
    c = connection.cursor()
    # Data containers
    text_dictionaries = []
    text_indices = []
    text_seeds = []
    # Extract each index
    for i in range(total_texts):
        # Parse json and append to data containers
        indexdata = json.loads(
            c.execute(f"SELECT data FROM info WHERE textid = '{i}'").fetchone()[0])
        text_dictionaries.append(indexdata[0])
        text_indices.append(indexdata[1])
        text_seeds.append(set(indexdata[2]))
        sys.stdout.write(f"{i} of {total_texts} indexes loaded\r")
        sys.stdout.flush()

    return text_dictionaries, text_indices, text_seeds, True


def clear_history(completed_file, result_directory):
    """Clears history from previous runs (result outputs and the file that ensures texts are not processed more than
    once)

    Args:
        completed_file (str): file with history of completed texts
        result_directory (str): directory with results stored in it
    """

    if os.path.isfile(completed_file):
        os.remove(completed_file)
    if os.path.exists(result_directory):
        shutil.rmtree(result_directory)
        os.mkdir(result_directory)


def set_files_to_ignore(completed, analysis_files, comparative_files):
    """Eliminate files that have already been completed from the analysis and comparative files

    Args:
        completed ([type]): [description]
        analysis_files ([type]): [description]
        comparative_files ([type]): [description]

    Returns:
        [type]: [description]
    """

    if len(completed) > 0:
        analysis_files = set(analysis_files) ^ set(completed)
        comparative_files = set(comparative_files) ^ set(completed)

    # In case there are empty strings (a possibility with user created file lists)
    # delete the empty strings from the list
    analysis_files = [a for a in analysis_files if a != ""]
    comparative_files = [c for c in comparative_files if c != ""]
    return analysis_files, comparative_files
# Prepare an index that links the titles of the files in the corpus with their place in
# the all texts list. Additionally, save the text lenghts for later data visualization


#****************#
# TEXT FUNCTIONS #
#****************#


def match_locations(source_locations, target_locations, source_dictionary, target_dictionary, matches):
    """Associate the locations of a seed to where it occurs in the second text
    Returns an ordered list of seeds from the source text, and a dictionary
    of the corrosponding locations in the target text.

    Args:
        source_locations ([type]): [description]
        target_locations ([type]): [description]
        source_dictionary ([type]): [description]
        target_dictionary ([type]): [description]
        matches ([type]): [description]

    Returns:
        [type]: [description]
    """
    all_source_locations = []
    locations_in_target = {}
    for match in matches:
        source_locs = source_locations[source_dictionary[match]]
        target_locs = target_locations[target_dictionary[match]]
        for source in source_locs:
            all_source_locations.append(source)
            locations_in_target[source] = sorted(target_locs)
    sorted_locations = sorted(all_source_locations)
    return sorted_locations, locations_in_target


def matchlocationsnonindexed(sourcedict, targetdict, matchingseeds):
    """The output is the same as match_location function but it runs on the index
    created by the get_seeds function.
    Args:
        sourcedict ([type]): [description]
        targetdict ([type]): [description]
        matchingseeds ([type]): [description]

    Returns:
        [type]: [description]
    """
    allsourcelocations = []
    locationsintarget = {}
    for seed in matchingseeds:
        sourcelocs = sourcedict[seed]
        targetlocs = targetdict[seed]
        for source in sourcelocs:
            allsourcelocations.append(source)
            locationsintarget[source] = targetlocs
    sortedlocations = sorted(allsourcelocations)
    return sortedlocations, locationsintarget

#


def extend(sourcetext, targettext, ss, ts, threshold, matchlength, maxlengthlev):
    """ Extend the two matching seeds until they fall below the set matching threshold.
    Returns their final length and the final similarity.

    Args:
        sourcetext ([type]): [description]
        targettext ([type]): [description]
        ss ([type]): [description]
        ts ([type]): [description]
        threshold ([type]): [description]
        matchlength ([type]): [description]
        maxlengthlev ([type]): [description]

    Returns:
        [type]: [description]
    """
    # determine the end of the strings
    se = ss + matchlength
    te = ts + matchlength

    # Make sure the end of the string does not extend past the end of the document in question
    if se >= len(sourcetext):
        se = len(sourcetext)
    if te >= len(targettext):
        te = len(targettext)

    # Get the string slices
    sourcestring = sourcetext[ss:se]
    targetstring = targettext[ts:te]

    # Measure initial similarity
    similarity = Levenshtein.ratio(sourcestring, targetstring)

    # Establish tracking information
    # How far has the quote extended?
    extender = 0
    # How many instances of straight decrease?
    straightdecrease = 0

    # Track the similarity in the last loop
    previoussimilarity = similarity

    # Save the final high similarity. I do this so I don't need to remeasure similarity after
    # backing the quote up.
    finalsimilarity = similarity

    # While similarity is above the matching threshold, and the end of the quotes are within the two texts
    # keep extending the matches
    while similarity >= threshold and se + extender <= len(sourcetext) and te + extender <= len(targettext):
        extender += 1

        # If the length is over a certain amount then limit the Lev measurement
        if (se + extender - ss >= maxlengthlev) and maxlengthlev:
            adjust = se + extender - ss - maxlengthlev
        else:
            adjust = 0

        # Check similarity of extended quote
        similarity = Levenshtein.ratio(
            sourcetext[ss+adjust:se+extender], targettext[ts+adjust:te+extender])

        # If the similarity is less than the previous instance, increment straight decrease
        # Otherwise, reset straight decrease to 0 and reset final similarity to similarity
        if similarity < previoussimilarity:
            straightdecrease += 1
        else:
            straightdecrease = 0
            finalsimilarity = similarity

        # Save similarity to previous similarity variable for use in the next iteration
        previoussimilarity = similarity
    # Back the length of the resulting quote up to the last time where its value began falling
    # and ended below the threshhold
    length = se+extender-straightdecrease - ss

    # return the length and final similarity
    return length, finalsimilarity


def alltextmatches(sourcelocations, targetlocationdict, sourcetext, targettext, matchlength, levthreshold, targettitle, maxcomp):
    """Find all of the matching quotes within two documents
    This function accepts the output of the match locations function
    And calls the extend function to extend and evaluate the quotes

    Args:
        sourcelocations ([type]): [description]
        targetlocationdict ([type]): [description]
        sourcetext ([type]): [description]
        targettext ([type]): [description]
        matchlength ([type]): [description]
        levthreshold ([type]): [description]
        targettitle ([type]): [description]
        maxcomp ([type]): [description]

    Returns:
        [type]: [description]
    """
    # container for results
    quoteinfo = []

    # To avoid duplcation of effort, track the extent of previously identified source quote
    # and end of the target quote
    sourcequoterange = [-1, -1]
    targetquoteend = -1

    # Iterate through each seed in the source lcation
    for source in sourcelocations:
        # If the current value of the source seed does not fall within the last identified quote range
        # reset the target quote end to search the whole text again.
        if not (source >= sourcequoterange[0] and source < sourcequoterange[1]):
            targetquoteend = -1

        # Get the corresponding locations
        targetlocations = targetlocationdict[source]

        # Iterate through all of the corresponding target locations
        for target in targetlocations:

            # Ensure the source start is not inside the last source quote and that the target seed
            # does not occur before the end of the last target quote.
            # If all of these things are true, this means the current seeds are internal to already
            # identified quotes and should be skipped
            if not (source >= sourcequoterange[0] and source < sourcequoterange[1] and target < targetquoteend):

                # Get the length and similarity of the strings begining at source and target
                length, similarity = extend(
                    sourcetext, targettext, source, target, levthreshold, matchlength, maxcomp)

                # Save the results if the length is above the minimum matching length and similarity
                # is above the set similarity threshold.
                if length >= matchlength and similarity >= levthreshold:
                    # Append the information to the data container as a tab seperated string
                    # This facilitate writing to file later
                    quoteinfo.append("\t".join([targettitle, str(length), str(similarity), str(source), str(
                        target), sourcetext[source:source+length], targettext[target:target+length]]))
                    # check if the start of the quote is inside the last source quote
                    # if it is not, reset the range to this quote. Otherwise, remain the same
                    if source >= sourcequoterange[1]:
                        sourcequoterange = [source, source+length]

                    # Set the end of the identified target quote
                    targetquoteend = target+length

    return quoteinfo


def comparetexts(sourcetext, sourcetitle, targetmeta, targettext, matchparams, sourceseeddict=None, index_data=None):
    """Compare two texts. This function is written this way to be fed in to the the
    multiprocessing starmap, which allows text comparisons to be conducted in paralle

    Args:
        sourcetext ([type]): [description]
        sourcetitle ([type]): [description]
        targetmeta ([type]): [description]
        targettext ([type]): [description]
        matchparams ([type]): [description]
        sourceseeddict ([type], optional): [description]. Defaults to None.
        index_data ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    # Construct the target title
    targettitle = targetmeta
    if index_data:
        # Retrieve global variables so the indices are accessible
        text_dictionaries, text_indices, text_seeds, title_to_index = index_data

        # get the index location the source and target document
        sourceindex = title_to_index[sourcetitle]
        targetindex = title_to_index[targettitle]

        # get the seed lists for both documents (stored as sets)
        sourceseedlist = text_seeds[sourceindex]
        targetseedlist = text_seeds[targetindex]

        # Find the interesection between the two sets
        matches = set.intersection(sourceseedlist, targetseedlist)

        # If the intersection between the two sets is more than zero,
        # get the seed dictionaries for both and provide this to the
        # match locations function
        if len(matches) > 0:
            # The source locations and target locations contain the index
            # locations for all the seeds
            source_locations = text_indices[sourceindex]
            target_locations = text_indices[targetindex]

            # These two dictionaries link the seeds themselves to the locations
            # in the above two lists
            source_dictionary = text_dictionaries[sourceindex]
            target_dictionary = text_dictionaries[targetindex]

            # get the matching locations
            sourcelocations, targetlocationdict = match_locations(
                source_locations, target_locations, source_dictionary, target_dictionary, matches)

    else:
        # Get the seeds from the source text dictionary (created outside the starmap process)
        # Starmap cannot accept a keys object, so it must be extracted here.
        sourceseedlist = sourceseeddict.keys()

        # Index the target text and get seeds
        targetseeddictionary = get_seeds(
            targettext, matchparams["seed_length"])
        targetseedlist = targetseeddictionary.keys()

        # Find the intersection between the sets and if there are matching seeds, gt locations
        matches = sourceseedlist & targetseedlist
        if len(matches) > 0:
            sourcelocations, targetlocationdict = matchlocationsnonindexed(
                sourceseeddict, targetseeddictionary, matches)

    # If there are matches, find all the text matches. Otherewise, return an empty list
    if len(matches) > 0:
        quoteinfo = alltextmatches(sourcelocations, targetlocationdict,
                                   sourcetext, targettext, matchparams["match_length"], matchparams["similarity_threshold"], targettitle, matchparams["max_comparisons"])
    else:
        quoteinfo = []
    return quoteinfo

#*********************#
# START OF MAIN LOGIC #
#*********************#


def run(cfg):

    # if clear history is set to true, clear the history:
    if cfg["run_settings"]["clear_history"]:
        clear_history(cfg["files"]["completed_files"],
                      cfg["paths"]["result_directory"])

    # load info
    all_titles, all_texts, title_to_index = load_data(
        cfg["files"]["data_file"],
        cfg["files"]["corpus_text_length_file"])

    # load analysis files
    # If there is a text to analyze file, load that in to memory
    if cfg["files"]["texts_to_analyze"]:
        analysisfiles = load_analysis_files(cfg["files"]["texts_to_analyze"])
    if not cfg["files"]["texts_to_analyze"] or len(analysisfiles) == 0:
        analysisfiles = all_titles

    # If there is a corpus composition file, load that into memory
    if cfg["files"]["corpus_composition"]:
        comparativefiles = load_analysis_files(
            cfg["files"]["corpus_composition"])
    if not cfg["files"]["corpus_composition"] or len(comparativefiles) == 0:
        comparativefiles = all_titles

    # check if any files have been completed already. these can be skipped
    jpm_eds = ['44975', '44974', '45765']
    completed = load_analysis_files(cfg["files"]["completed_files"])
    
    # remove files to be skipped
    analysisfiles, comparativefiles = set_files_to_ignore(
        completed, analysisfiles, comparativefiles)


    comparativefiles = [f for f in comparativefiles if f not in jpm_eds]
    
    # load index files if used
    if cfg["files"]["index_file"]:
        text_dictionaries, text_indices, text_seeds, prepped_index = load_index(
            cfg["files"]["index_file"], len(all_titles))
        index_data = (text_dictionaries, text_indices,
                      text_seeds, title_to_index)
    else:
        prepped_index = False

    # Create an empty list to keep track of how long the program is running
    runtimes = []

    # Initialize thread pool for parallel processing
    pool = Pool(60) #maxtasksperchild=cfg["run_settings"]["max_child_tasks"])

    # If the RESULT_DIRECTORY does not exist, create it
    if not os.path.exists(cfg["paths"]["result_directory"]):
        os.mkdir(cfg["paths"]["result_directory"])

    # Track how many files have been completed
    total_completed_files = len(completed)

    # If there is no prepared index, order the texts by size.
    # This speeds up calculation times by ensuring that long indices are not calculated
    # more than necessary. This can optionally be turned off.
    if not prepped_index and cfg["run_settings"]["front_load"]:
        analysislengths = {}
        for f in analysisfiles:
            textlocation = title_to_index[f]
            text = all_texts[textlocation]
            analysislengths[f] = len(text)
        analysisfiles = sorted(
            analysislengths, key=analysislengths.get, reverse=True)

        comparativelengths = {}
        for c in comparativefiles:
            textlocation = title_to_index[c]
            text = all_texts[textlocation]
            comparativelengths[c] = len(text)
        comparativefiles = sorted(
            comparativelengths, key=comparativelengths.get, reverse=True)

    # Iterate through every analysis file.
    for f in analysisfiles:
        # Record start time
        s = time.time()
        
        # Gather text information
        textlocation = title_to_index[f]
        text = all_texts[textlocation]
        title = all_titles[textlocation]

        # Get the comparative info. This needs to be recalculated each loop so
        # information about finished texts is not sent into the analysis function
        comparativemeta = []
        comparativetexts = []
        for cf in comparativefiles:
            if cf != f:
                location = title_to_index[cf]
                comparativemeta.append(all_titles[location])
                comparativetexts.append(all_texts[location])

        # Print current status
        print(
            f"Analyzing text {title}: {total_completed_files + 1} vs {len(comparativetexts)}: {f} from {all_titles[textlocation]} (length: {len(text)})")

        # If the index has been prepped, give the data to the comparetexts function via starmap.
        # Otherwise, create an index for the source text and then give the data
        if prepped_index:
            results = pool.starmap(comparetexts, zip(repeat(text),
                                                     repeat(title),
                                                     comparativemeta,
                                                     comparativetexts,
                                                     repeat(cfg["params"]),
                                                     index_data=repeat(index_data)))
        else:
            sourcedict = get_seeds(text, cfg["params"]["seed_length"])
            results = pool.starmap(comparetexts, zip(repeat(text),
                                                     repeat(title),
                                                     comparativemeta,
                                                     comparativetexts,
                                                     repeat(cfg["params"]),
                                                     repeat(sourcedict)))

        # Filter out results where nothing was returned and flatten the list
        # because results returns a list of lists
        filtered_results = [r for r in results if len(r) != 0]
        filtered_results = list(itertools.chain(*filtered_results))

        # write the results to file:
        with open(os.path.join(cfg["paths"]["result_directory"], f"{f}.txt"), "w", encoding="utf8") as wf:
            wf.write(
                "TargetTitle\tLength\tratio\tSource place\tTarget place\tAnalysis text\tTarget text\n")
            wf.write("\n".join(filtered_results))

        # Append the analysis filename to the completed lists and remove it from the comparative files
        completed.append(f)
        if f in comparativefiles:
            comparativefiles.remove(f)
        # if f in analysisfiles:
        #     analysisfiles.remove(f)

        # write the completed list to file:
        with open("completed_files.txt", "w", encoding="utf8") as wf:
            wf.write("\n".join(completed))
            wf.close()

        # Increment the completed files tracker
        total_completed_files += 1

        # Record the finished time, append to the runtimes list, and calculated time usage information
        e = time.time()
        t = e-s
        runtimes.append(t)

        total = sum(runtimes)
        average = total/len(runtimes)

        # Print information on time usage.
        print(
            f"Operation completed in {t:.2f} seconds (averaging {average:.2f}, in total {total:.2f})")


if __name__ == "__main__":
    # read config file
    cfg = read_config_file("config/intertextuality_config.json")
    # run the analysis
    run(cfg)
