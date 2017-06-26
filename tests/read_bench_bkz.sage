import os
import numpy as np

variables = ["dims", "blocks", "N", "result_dir", "raw_out_files"]

test_script = "bench_bkz.sage"
with open(test_script) as f:
	for line in f:
		for var in variables:
			if line.startswith(var):
				exec(line)

def get_key_value(key, line):
	for l in line.strip().split(','):
		if l.strip().startswith(key):
			return l.split('=')[-1].strip()

def get_data(fn, keys):
	data = []
	with open(fn) as f:
		d = {}
		for line in f:
			if line.startswith("End of BKZ loop"):
				for key in keys:
					d[key] = get_key_value(key, line)
			if line.startswith("End of BKZ: success"):
				data.append(copy(d))
				d = {}
			if "error" in line:
				data.append({})
	return data

keys = ["time", "r_0", "log2(nodes)"]

data = {}
avg_data = {}
for dim in dims:
	data[dim] = {}
	avg_data[dim] = {}
	for block in blocks:
		data[dim][block] = {}
		avg_data[dim][block] = {}
		for n in range(N):
			dat = get_data(os.path.join(result_dir, raw_out_files.format(dim,block,n)), keys)
			for d in dat:
				if d:
					d["time"] = float(d["time"][:-1])
					d["r_0"] = float(d["r_0"])
					d["log2(nodes)"] = float(d["log2(nodes)"])
			data[dim][block][n] = dat
		for k in keys:
			avg_data[dim][block][k] = np.average([data[dim][blocks][i][k] for i in range(N)])

