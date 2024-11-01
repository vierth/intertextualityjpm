'''
This script takes every file in a folder named 'corpus' and transforms it into
a serialized object that can be used to study intertextuality at a large level.
This is the first step in the analysis workflow. It removes unwanted punctuation
and other artifacts from the texts, parses metadata from the file titles and saves
it into a "pickle" file. Note that all pickle files in this workflow are generated
by these scripts. If you do not trust the source of your pickle files,  you should
never open them, as they can contain arbitrary code! This assumes you are using the
Anaconda distribution of Python 3.

Assumed Filename format is this:
Title-Era-Author_TextDivison.txt

TextDivision should be zero if the text has not been divided. Otherwise, it can be
chapter number, volume number, etc.
'''

import os
import pickle
import sys

from config.config import read_config_file
from utilities.clean import clean


def run():
    with open("../metainfo/refined_meta.tsv", 'r', encoding='utf8') as rf:
        text = rf.read().split("\n")
        data = [l.split("\t") for l in text]
        use_t = [d[0] for d in data]
        

    # get configuration settings
    cfg = read_config_file("config/prepare_config.json")

    # containers for the data. The final file will be a list of lists. The first item
    # will be metadata, and the second item will be the texts themselves.
    all_filenames = []
    all_texts = []

    # Track the total length of the corpus.
    totalcharacters = 0

    # Retrieve the directory information
    for folder in cfg["paths"]["input_folders"]:
        for root, dirs, files in os.walk(folder):

            # remove license file if present
            files = [f for f in files if f != "LICENSE" and f[0] != "."]
            print(files)
            files = [f for f in files if f[:-4] in use_t]
            # Iterate through each file in the filelist
            for i, f in enumerate(files):

                # remove extension
                file_name = f[:-4]

                with open(f"{root}/{f}", 'r', encoding='utf8') as tf:
                    # clean text, append it to all texts, and increment total length tracker
                    text = clean(
                        tf.read(), cfg["params"]["characters_to_remove"], cfg["params"]["simplify"])
                    all_texts.append(text)
                    totalcharacters += len(text)

                    # Save metadata
                    all_filenames.append(file_name)

                # print tracking statement
                sys.stdout.write(f"{i + 1} documents of {len(files)} completed\r")
                sys.stdout.flush()

    # Print basic corpus information
    print(f"\n{totalcharacters} from {len(all_texts)} documents.")

    # Save data to pickle
    pickle.dump([all_filenames, all_texts], open(
        cfg["files"]["output_file"], "wb"))

    # write meta file
    with open(cfg["files"]["content_meta"], 'w', encoding='utf8') as wf:
        wf.write("\n".join(all_filenames))


if __name__ == "__main__":
    run()
