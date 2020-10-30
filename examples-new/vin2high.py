'''
v.in transform between low and high precision
'''

import os
import shutil


def vin2high(vin_file):
	"""to high precision"""
	vin_high = vin_file + '.high'
	with open(vin_file, 'r') as fi, open(vin_high, 'w') as fo:
		line = fi.readline()
		idx = 0
		while line:
			tokens = line.strip().split()
			tokens = list(map(float, tokens))

			if (idx + 1) % 3 == 0:
				fo.write(' '*3)
				fo.write(''.join(['%8d' % t for t in tokens]))
				fo.write('\n')
			else:
				fo.write('%2d ' % tokens[0])
				fo.write(''.join(['%8.3f' % t for t in tokens[1:]]))
				fo.write('\n')

			line = fi.readline()
			idx += 1


def vin2low(vin_file):
	"""to low precision"""
	vin_high = vin_file + '.low'
	with open(vin_file, 'r') as fi, open(vin_high, 'w') as fo:
		line = fi.readline()
		idx = 0
		while line:
			tokens = line.strip().split()
			tokens = list(map(float, tokens))

			if (idx + 1) % 3 == 0:
				fo.write(' '*3)
				fo.write(''.join(['%7d' % t for t in tokens]))
				fo.write('\n')
			else:
				fo.write('%2d ' % tokens[0])
				fo.write(''.join(['%7.2f' % t for t in tokens[1:]]))
				fo.write('\n')

			line = fi.readline()
			idx += 1


if __name__ == '__main__':
	CUR_DIR = os.path.dirname(os.path.abspath(__file__))
	examples = ['e' + str(i) for i in range(1, 9)]
	for e in examples:
		vin = os.path.join(CUR_DIR, e, 'v.in')
		vin2high(vin)
		vin2low(vin)
		# shutil.copy(vin, vin + '.low')
		shutil.copy(vin + '.high', vin)
