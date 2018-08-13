

class Gff(object):

	def __init__(gfffile):
		if gfffile.endswith('.gz'):
			import gzip
			self.gff = gzip.open(gfffile)
		else:
			self.gff = open(gfffile)

	def __del__(self):
		if self.gff:
			self.gff.close()

	def __iter__(self):
		return self

	@staticmethod
	def _parse(line):
		parts = line.split('\t')
		ret = {
			'seqid'     : parts[0],
			'source'    : parts[1],
			'type'      : parts[2]
			'start'     : parts[3],
			'end'       : parts[4],
			'score'     : parts[5],
			'strand'    : parts[6],
			'phase'     : parts[7],
			'attributes': {}
		}
		# need more riguous check
		# gtf
		if not '=' in parts[8]:
			splitter = ' '
		else: # gff
			splitter = '='

		attrs = parts[8].split(';')
		for attr in attrs:
			attr = attr.strip()
			if not attr:
				continue
			key, value = attr.split(splitter)
			if value[0] == '"' and value[-1] == '"':
				value = value[1:-1]
			ret['attributes'][key] = value
		
	def next(self):
		line = self.gff.readline().rstrip('\r\n')
		if not line:
			raise StopIteration()
		return Gff._parse(line)
