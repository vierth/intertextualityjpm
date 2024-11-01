import re

from hanziconv import HanziConv


# clean the text. This will remove everything in the above list from the text
def clean(content, remove, simplified=True):
    # These two lines are useful for Chinese texts where there was no whitespace or punctuation
    # in the original documents
    content = re.sub('\s+', '', content)
    content = content.translate(str.maketrans('', '', "".join(remove)))

    if simplified:
        content = HanziConv.toSimplified(content)

    return content
