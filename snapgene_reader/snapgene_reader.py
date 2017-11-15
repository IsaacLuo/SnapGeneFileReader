import struct
import xmltodict
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
import html2text


html_parser = html2text.HTML2Text()
html_parser.ignore_emphasis = True
html_parser.ignore_links = True
html_parser.body_width = 0

def parse(val):
    return html_parser.handle(val).strip() if isinstance(val, str) else val

def snapgene_file_to_dict(filepath=None, fileobject=None):
    """Return a dictionnary containing the data from a ``*.dna`` file.

    Parameters
    ----------
    filepath
      Path to a .dna file created with SnapGene
    fileobject
      On object-like pointing to the data of a .dna file created with SnapGene
    """

    if filepath is not None:
        fileobject = open(filepath, 'rb')

    if fileobject.read(1) != b'\t':
        raise ValueError('Wrong format for a SnapGene file !')

    def unpack(size, mode):
        return struct.unpack('>' + mode, fileobject.read(size))[0]

    # READ THE DOCUMENT PROPERTIES
    length = unpack(4, 'I')
    title = fileobject.read(8).decode('ascii')
    if length != 14 or title != 'SnapGene':
        raise ValueError('Wrong format for a SnapGene file !')

    data = dict(
        isDNA=unpack(2, 'H'),
        exportVersion=unpack(2, 'H'),
        importVersion=unpack(2, 'H'),
        features=[]
    )

    while True:
        # READ THE WHOLE FILE, BLOCK BY BLOCK, UNTIL THE END
        next_byte = fileobject.read(1)

        if next_byte == b'':
            # END OF FILE
            break

        block_size = unpack(4, 'I')

        if ord(next_byte) == 0:
            # READ THE SEQUENCE AND ITS PROPERTIES
            props = unpack(1, 'b')
            data["dna"] = dict(
                topology="circular" if props & 0x01 else "linear",
                strandedness="double" if props & 0x02 > 0 else "single",
                damMethylated=props & 0x04 > 0,
                dcmMethylated=props & 0x08 > 0,
                ecoKIMethylated=props & 0x10 > 0,
                length=block_size - 1
            )
            data["seq"] = fileobject.read(block_size - 1).decode('ascii')

        elif ord(next_byte) == 10:
            # READ THE FEATURES
            strand_dict = {"0": ".", "1": "+", "2": "-", "3": "="}
            format_dict = {'@text': parse, '@int': int}
            features_data = xmltodict.parse(fileobject.read(block_size))

            for f in features_data["Features"]["Feature"]:
                segments = f["Segment"]
                if not isinstance(segments, list):
                    segments = [segments]
                segments_ranges = [
                   sorted([int(e) for e in segment['@range'].split('-')])
                   for segment in segments
                ]
                qualifiers = f.get('Q', [])
                if not isinstance(qualifiers, list):
                    qualifiers = [qualifiers]
                parsed_qualifiers = {}
                for q in qualifiers:
                    if isinstance(q['V'], list):
                        if len(q['V'][0].items()) == 1:
                            parsed_qualifiers[q['@name']] = l = []
                            for e in q['V']:
                                fmt, value = e.popitem()
                                fmt = format_dict.get(fmt, parse)
                                l.append(fmt(value))
                        else:
                            parsed_qualifiers[q['@name']] = d = {}
                            for e in q['V']:
                                (fmt1, value1), (fmt2, value2) = e.items()
                                fmt = format_dict.get(fmt1, parse)
                                d[value2] = fmt(value1)

                    else:
                        fmt, value = q['V'].popitem()
                        fmt = format_dict.get(fmt, parse)
                        parsed_qualifiers[q['@name']] = fmt(value)

                if 'label' not in parsed_qualifiers:
                    parsed_qualifiers['label'] = f['@name']
                if 'note' not in parsed_qualifiers:
                    parsed_qualifiers['note'] = []
                if not isinstance(parsed_qualifiers['note'], list):
                    parsed_qualifiers['note'] = [parsed_qualifiers['note']]
                color = segments[0]['@color']
                parsed_qualifiers['note'].append("color: " + color)

                data["features"].append(dict(
                    start=min([start for (start, end) in segments_ranges]),
                    end=max([end for (start, end) in segments_ranges]),
                    strand=strand_dict[f.get('@directionality', "0")],
                    type=f['@type'],
                    name=f['@name'],
                    color=segments[0]['@color'],
                    textColor='black',
                    segments=segments,
                    row=0,
                    isOrf=False,
                    qualifiers=parsed_qualifiers
                ))

        else:
            # WE IGNORE THE WHOLE BLOCK
            fileobject.read(block_size)

    fileobject.close()

    return data

def snapgene_file_to_seqrecord(filepath=None, fileobject=None):
    """Return a BioPython SeqRecord from the data of a ``*.dna`` file.

    Parameters
    ----------
    filepath
      Path to a .dna file created with SnapGene
    fileobject
      On object-like pointing to the data of a .dna file created with SnapGene
    """
    data = snapgene_file_to_dict(filepath=filepath, fileobject=fileobject)
    strand_dict = {'+': 1, '-': -1, '.': 0}

    return SeqRecord(
        seq=Seq(data['seq'], alphabet=DNAAlphabet()),
        features=[
            SeqFeature(
                location=FeatureLocation(
                    start=feature['start'],
                    end=feature['end'],
                    strand=strand_dict[feature['strand']]
                ),
                strand=strand_dict[feature['strand']],
                type=feature['type'],
                qualifiers= feature['qualifiers']
            )
            for feature in data['features']
        ]
    )
