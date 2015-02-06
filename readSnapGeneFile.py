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


def SnapGeneFileToJson(fileName):
	s = SnapGeneFileReader()
	s.openFile(fileName)
	b = s.readFirstByte()
	if b!=9:
		print "nota snapgene file"
	else:
		s.readDescriptionSegment()

	dct = {}

	try:
		while s.sgFile:
			b = s.readFirstByte()
			if b==0:
				s.readDnaSegment()
			elif b==8:
				print "additional sequence properties seqment"
				s.readAdditionalSequenceProperties()
			elif b==5:
				print "primers"
				s.passThisDataBlock()
			elif b==6:
				print "notes"
				s.passThisDataBlock()
			elif b==17:
				print "alignable sequence"
				s.passThisDataBlock()
			elif b==16:
				print "alignable sequence"
				s.passThisDataBlock()
			elif b==7:
				print "history tree"
				s.passThisDataBlock()
			elif b==11:
				print "history node"
				s.passThisDataBlock()
			elif b==1:
				print "compressed DNA"
				s.passThisDataBlock()
			elif b==10:
				print "features"
				dct = s.readFeatures()
			elif b==18:
				print "sequence trace"
				s.passThisDataBlock()
			else:
				print b
				s.passThisDataBlock()
			
	except:
		print "end"
	return s.data

		


class SnapGeneFileReader:
	def openFile(self, fileName):
		self.sgFile = open(fileName,"rb") 
		self.data = {}


	def readFirstByte(self):
		s = self.sgFile.read(1)
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
		print self.data
	
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
		d["sequence"] = self.sgFile.read(d["length"])
		print d

	
	def readAdditionalSequenceProperties(self):
		dataLength = struct.unpack('>I', self.sgFile.read(4))[0]
		data = self.sgFile.read(dataLength)
		print data

	def readFeatures(self):
		dataLength = struct.unpack('>I', self.sgFile.read(4))[0]
		data = self.sgFile.read(dataLength)
		dct = xmltodict.parse(data)
		self.data["features"] = dct["Features"]["Feature"]
		print dct
		return dct

	def passThisDataBlock(self):
		dataLength = struct.unpack('>I', self.sgFile.read(4))[0]
		data = self.sgFile.read(dataLength)
		#print "passed %d"%dataLength
	



if __name__ == '__main__':
	print SnapGeneFileToJson("pIB2-SEC13-mEGFP.dna")
	


		
		

