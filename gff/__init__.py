

class Gff(object):

	def __init__(self, gfffile):
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
	def toBed(gfffile, name = '{attributes[gene_id]}'):
		gff = Gff(gfffile)
		for g in gff:
			yield {
				'CHR'   : g['seqid'],
				'START' : g['start'],
				'END'   : g['end'],
				'NAME'  : name(g) if callable(g) else name.format(**g),
				'SCORE' : g['score'],
				'STRAND': g['strand'],
			}

	@staticmethod
	def toBedfile(gfffile, bedfile, name = '{g[attributes][gene_id]}'):
		with open(bedfile, 'w') as f:
			for g in Gff.toBed(gfffile, name):
				f.write("{CHR}\t{START}\t{END}\t{NAME}\t{SCORE}\t{STRAND}\n".format(**g))

	@staticmethod
	def _parse(line):
		parts = line.split('\t')
		ret = {
			'seqid'     : parts[0],
			'source'    : parts[1],
			'type'      : parts[2],
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
		return ret


	def next(self):
		line = self.gff.readline().rstrip('\r\n')
		if not line:
			raise StopIteration()
		return Gff._parse(line)

	__next__ = next