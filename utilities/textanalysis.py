"""
Utility functions for the text analysis
"""


def get_seeds(text, n):
    """
    create seed index for text of given size n

    Args:
        text (str): text to divide into seed
        n (int): length of seed

    Returns:
        dict: dictionary containing index
    """
    seed_dict = {}
    for i in range(len(text)-n):
        seed = text[i:i+n]
        if seed in seed_dict:
            seed_dict[seed].append(i)
        else:
            seed_dict[seed] = [i]
    return seed_dict
