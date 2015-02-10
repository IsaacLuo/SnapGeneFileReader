#SnapGene .dna files store information about the DNA, features, primers, notes, history, display options, and
#more. This document will not cover the file segments that are typically irrelevant for applications other than
#SnapGene. For example, rapid searching of a .dna file is facilitated by indexes, including a MICA index (MICA:
#desktop software for comprehensive searching of DNA databases. Stokes WA, Glick BS. BMC Bioinformatics.
#2006 Oct 3;7:427). Here, we will describe the following file segments:
#File Description
#DNA
#Additional Sequence Properties
#Features
#Primers
#Notes

import struct
import xmltodict
import json
import sys
import os
import os.path


def SnapGeneFileToJson(fileName):
	err = sys.stderr
	s = SnapGeneFileReader()
	s.openFile(fileName)
	b = s.readFirstByte()
	if b!=9:
		err.write("not a snapgene file\n")
	else:
		s.readDescriptionSegment()

	try:
		while s.sgFile:
			b = s.readFirstByte()
			if b=="":
				break
			elif b==0:
				s.readDnaSegment()
			elif b==8:
				err.write("additional sequence properties seqment\n")
				s.readAdditionalSequenceProperties()
			elif b==5:
				err.write("primers\n")
				s.passThisDataBlock()
			elif b==6:
				err.write("notes\n")
				s.passThisDataBlock()
			elif b==17:
				err.write("alignable sequence\n")
				s.passThisDataBlock()
			elif b==16:
				err.write("alignable sequence\n")
				s.passThisDataBlock()
			elif b==7:
				err.write("history tree\n")
				s.passThisDataBlock()
			elif b==11:
				err.write("history node\n")
				s.passThisDataBlock()
			elif b==1:
				err.write("compressed DNA\n")
				s.passThisDataBlock()
			elif b==10:
				err.write("features\n")
				s.readFeatures()
			elif b==18:
				err.write("sequence trace\n")
				s.passThisDataBlock()
			else:
				s.passThisDataBlock()
			
	except Exception as e: 
		#err.write(e+"\n")
		#print "=============***++++======="
		#print e
		#print "=============***++++======="
		err.write(e.message)
		err.write("\n")

	return s.data

		


class SnapGeneFileReader:
	def openFile(self, fileName):
		self.sgFile = open(fileName,"rb") 
		self.data = {"name":os.path.splitext(os.path.basename(fileName))[0]}


	def readFirstByte(self):
		s = self.sgFile.read(1)
		if s == "":
			return s
		b = ord(s)
		return b
		
	def readDescriptionSegment(self):
		length = struct.unpack('>I', self.sgFile.read(4))[0]
		title = self.sgFile.read(8)
		if length != 14 or title != "SnapGene":
			raise "not a snapgene file"

		self.data["isDNA"] = struct.unpack('>H', self.sgFile.read(2))[0]
		self.data["exportVersion"] = struct.unpack('>H', self.sgFile.read(2))[0]
		self.data["importVersion"] = struct.unpack('>H', self.sgFile.read(2))[0]
	
	def readDnaSegment(self):
		self.data["dna"] = {}
		d = self.data["dna"]
		dataLength = struct.unpack('>I', self.sgFile.read(4))[0]
		pType = struct.unpack('>b', self.sgFile.read(1))[0]
		d["topology"] = "circular" if pType & 0x01 else "linear"
		d["strandedness"] = "double" if pType & 0x02 > 0 else "single"
		d["damMethylated"] = pType & 0x04 >0
		d["dcmMethylated"] = pType & 0x08 >0
		d["ecoKIMethylated"] = pType & 0x10 >0
		d["length"] = dataLength - 1
		self.data["seq"] = self.sgFile.read(d["length"])

	
	def readAdditionalSequenceProperties(self):
		dataLength = struct.unpack('>I', self.sgFile.read(4))[0]
		data = self.sgFile.read(dataLength)

	def readFeatures(self):
		dataLength = struct.unpack('>I', self.sgFile.read(4))[0]
		data = self.sgFile.read(dataLength)
		dct = xmltodict.parse(data)
		#self.data["features"] = dct["Features"]["Feature"]
		original = dct["Features"]["Feature"]
		self.data['features'] = []

		strandDict = {"0":".","1":"+","2":"-","3":"="}

		for f in original:
			d = {}
			d['text']=f['@name']
			d['textColor']='black'
			directionality = f['@directionality'] if f.has_key('@directionality') else "0" 
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
				if r[0]<totalStart:
					totalStart = seg['start']
				if r[1]>totalEnd:
					totalEnd = seg['end'] 

				seg['color'] = s['@color']
				seg['name'] = s['@name'] if s.has_key('@name') else ""
				d['segments'].append(seg)

			d['start'] = totalStart
			d['end'] = totalEnd
			d['color'] = segments[0]['@color']

			d['row'] = 0
			d['isORf'] = False

			d['qualifiers'] = []
			if f.has_key("Q"):
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



	def passThisDataBlock(self):
		dataLength = struct.unpack('>I', self.sgFile.read(4))[0]
		data = self.sgFile.read(dataLength)
		#print "passed %d"%dataLength
	



if __name__ == '__main__':
	if len(sys.argv)>1:
		print json.dumps(SnapGeneFileToJson(sys.argv[1]))
	else:
		print json.dumps(SnapGeneFileToJson("pIB2-SEC13-mEGFP.dna"))
	


		
		

