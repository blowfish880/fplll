import os
import subprocess

cmd = "./test_bkz"

result_dir = "results"
raw_out_files = "raw_{0}_{1}_{2}.txt"

dims = [80]
blocks = [35]
N = 2

if not os.path.isdir(result_dir):
	os.makedirs(result_dir)

for dim in dims:
	for block in blocks:
		for n in range(N):
			with open(os.path.join(result_dir, raw_out_files.format(dim, block, n)), 'w') as f:
				subprocess.call([cmd , str(dim), str(block), str(n)], stdout=f, stderr=f)

