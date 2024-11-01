import json


def read_config_file(config_file, type_of_config="general"):
    with open(config_file, 'r', encoding="utf8") as rf:
        return json.load(rf)
