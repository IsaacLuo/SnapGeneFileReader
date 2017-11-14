import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation

from SnapGeneFileReader import SnapGeneFileReader

def parse_qualifiers(qualifiers_list):
    qualifiers_dict = {}
    for q in qualifiers_list:
        name = q['name']
        val = q['value']
        for key in val:
            try:
                qualifiers_dict[name] = val[key]
            except:
                qualifiers_dict[name] = 'unknown'
            break
    return qualifiers_dict

def snapgene_file_to_seq_record(file_name):
    obj_dict = SnapGeneFileReader().read_file(file_name)

    record = SeqRecord(Seq(obj_dict['seq'], DNAAlphabet() ))

    strand_dict = {'+':1, '-':-1, '.':0}

    for feature in obj_dict['features']:
        record.features.append(
            SeqFeature(
                FeatureLocation(
                    feature['start'], 
                    feature['end'],
                    strand=strand_dict[feature['strand']],
                    ),
                strand=strand_dict[feature['strand']],
                id=feature['text'],
                type=feature['type'],
                qualifiers=parse_qualifiers(feature['qualifiers'])
                )
        )

    record.name = record.id = obj_dict['name']

    return record

