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
			value = l.split('=')[-1].strip()
			return float(value[:-1]) if key == "time" else float(value)

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
				for key in keys:
					d[key] = NaN
				data.append(copy(d))
				d = {}
	return data

def read_dir(fld, fns, keys, dims, blocks, N):
	data = {}
	avg_data = {}
	for dim in dims:
		data[dim] = {}
		avg_data[dim] = {}
		for block in blocks:
			data[dim][block] = {}
			avg_data[dim][block] = {}
			for n in range(N):
				dat = get_data(os.path.join(fld, fns.format(dim,block,n)), keys)
				#~ for d in dat:
					#~ if d:
						#~ d["time"] = float(d["time"][:-1])
						#~ d["r_0"] = float(d["r_0"])
						#~ d["log2(nodes)"] = float(d["log2(nodes)"])
				data[dim][block][n] = dat
			for k in keys:
				avg_data[dim][block][k] = [np.average(filter(lambda x: x is not NaN, [data[dim][block][n][flag][k] for n in range(N)])) for flag in range(4)]
	return data, avg_data

keys = ["time", "r_0", "log2(nodes)"]
flags = ["no", "blll", "nolll", "b_no_lll"]
colors = ["blue", "red", "black", "green"]
plt_dir = result_dir + "_plts"
plt_name = "plt_{0}_{1}.pdf"

if not os.path.isdir(plt_dir):
	os.makedirs(plt_dir)

data, avg_data = read_dir(result_dir, raw_out_files, keys, dims, blocks, N)
for dim in dims:
	for k in keys:
		plt = list_plot([])
		for f in range(3): 
			plt += line2d([(block, avg_data[dim][block][k][f]) for block in blocks], color=colors[f], legend_label=flags[f])
		d3 = filter(lambda x: x[1] is not NaN, [(block, avg_data[dim][block][k][3]) for block in blocks])
		plt += points(d3, color=colors[3], legend_label=flags[3])
		plt.axes_labels(["block", k])
		plt.save(os.path.join(plt_dir, plt_name.format(dim,k)))
