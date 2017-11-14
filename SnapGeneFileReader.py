import struct
import xmltodict
import json
import sys
import os
import os.path

class SnapGeneFileReader:
    def _open_file(self, file_name):
        self.sg_file = open(file_name,"rb") 
        self.data = {"name":os.path.splitext(os.path.basename(file_name))[0]}

    def _open_file_object(self, file_object, name='unknown'):
        self.sg_file = file_object
        self.data = {"name": name}

    def _read_first_byte(self):
        s = self.sg_file.read(1)
        if s == b'':
            return s
        b = ord(s)
        return b

    def _read_description_segment(self):
        length = struct.unpack('>I', self.sg_file.read(4))[0]
        title = self.sg_file.read(8).decode('ascii')
        if length != 14 or title != 'SnapGene':
            raise "not a snapgene file"

        self.data["isDNA"] = struct.unpack('>H', self.sg_file.read(2))[0]
        self.data["exportVersion"] = struct.unpack('>H', self.sg_file.read(2))[0]
        self.data["importVersion"] = struct.unpack('>H', self.sg_file.read(2))[0]

    def _read_dna_segment(self):
        self.data["dna"] = {}
        d = self.data["dna"]
        dataLength = struct.unpack('>I', self.sg_file.read(4))[0]
        pType = struct.unpack('>b', self.sg_file.read(1))[0]
        d["topology"] = "circular" if pType & 0x01 else "linear"
        d["strandedness"] = "double" if pType & 0x02 > 0 else "single"
        d["damMethylated"] = pType & 0x04 >0
        d["dcmMethylated"] = pType & 0x08 >0
        d["ecoKIMethylated"] = pType & 0x10 >0
        d["length"] = dataLength - 1
        self.data["seq"] = self.sg_file.read(d["length"]).decode('ascii')

    def _read_additional_sequence_properties(self):
        dataLength = struct.unpack('>I', self.sg_file.read(4))[0]
        data = self.sg_file.read(dataLength)

    def _read_features(self):
        dataLength = struct.unpack('>I', self.sg_file.read(4))[0]
        data = self.sg_file.read(dataLength)
        dct = xmltodict.parse(data)
        #self.data["features"] = dct["Features"]["Feature"]
        original = dct["Features"]["Feature"]
        self.data['features'] = []

        strandDict = {"0":".","1":"+","2":"-","3":"="}

        for f in original:
            d = {}
            d['text']=f['@name']
            d['textColor']='black'
            directionality = f['@directionality'] if '@directionality' in f else "0" 
            d['strand']= strandDict[directionality]
            d['type'] = f['@type']
            d['segments'] = []
            if type(f["Segment"]) == type([]):
                segments = f["Segment"]
            else:
                segments = [f["Segment"],]

            r = segments[0]["@range"].split('-')
            totalStart = int(r[0])
            totalEnd = int(r[1])

            for s in segments:
                seg = {}
                r = s["@range"].split('-')
                seg['start'] = int(r[0])
                seg['end'] = int(r[1])
                if seg['start']<totalStart:
                    totalStart = seg['start']
                if seg['end']>totalEnd:
                    totalEnd = seg['end'] 

                seg['color'] = s['@color']
                seg['name'] = s['@name'] if '@name' in s else ""
                d['segments'].append(seg)

            d['start'] = totalStart
            d['end'] = totalEnd
            d['color'] = segments[0]['@color']

            d['row'] = 0
            d['isORf'] = False

            d['qualifiers'] = []
            if "Q" in f:
                if type(f["Q"]) == type([]):
                    qualifiers = f["Q"]
                else:
                    qualifiers = [f["Q"],]
                for q in qualifiers:
                    qu = {}
                    qu['name'] = q['@name']
                    qu['value'] = q['V']
                    d['qualifiers'].append(qu)

            self.data['features'].append(d)

    def _pass_this_data_block(self):
        dataLength = struct.unpack('>I', self.sg_file.read(4))[0]
        data = self.sg_file.read(dataLength)
        #print "passed %d"%dataLength

    def read_file(self, file_name, file_object=None):
        """
        read the snapgene file into a python dict.

        Args:
            fileName: the file name of a file path to a snapgene dna file, or it
                can be just a string name if fileObject is provided.
            fileObject: if the file is already open, pass the file object into the
                function and give it a name, the fileObject can be a io.BytesIO object
        Retruns:
            A dict which has name, seq, features and other information read from the file.
        """
        err = sys.stderr
        if file_object:
            self._open_file_object(file_object, file_name)
        else:
            self._open_file(file_name)

        b = self._read_first_byte()
        if b != 9:
            err.write("not a snapgene file\n")
            raise RuntimeError('not a snapgene file')
        else:
            self._read_description_segment()

        try:
            while self.sg_file:
                b = self._read_first_byte()
                if b==b'':
                    break
                elif b==0:
                    self._read_dna_segment()
                elif b==8:
                    # err.write("additional sequence properties seqment\n")
                    self._read_additional_sequence_properties()
                elif b==5:
                    # err.write("primers\n")
                    self._pass_this_data_block()
                elif b==6:
                    # err.write("notes\n")
                    self._pass_this_data_block()
                elif b==17:
                    # err.write("alignable sequence\n")
                    self._pass_this_data_block()
                elif b==16:
                    # err.write("alignable sequence\n")
                    self._pass_this_data_block()
                elif b==7:
                    # err.write("history tree\n")
                    self._pass_this_data_block()
                elif b==11:
                    # err.write("history node\n")
                    self._pass_this_data_block()
                elif b==1:
                    # err.write("compressed DNA\n")
                    self._pass_this_data_block()
                elif b==10:
                    # err.write("features\n")
                    self._read_features()
                elif b==18:
                    # err.write("sequence trace\n")
                    self._pass_this_data_block()
                else:
                    self._pass_this_data_block()

            self.sg_file.close()
            
        except Exception as e: 
            err.write(str(e))
            err.write("\n")
            raise

        return self.data