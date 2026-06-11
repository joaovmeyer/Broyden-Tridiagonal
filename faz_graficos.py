import matplotlib.pyplot as plt
import subprocess
import os


Ns = [32, 64, 128, 256, 512, 1000, 2000, 4000, 8000, 9000, 10000, 20000]


resultados = {}
for path in os.listdir("./"):
	if not path.startswith("versao") or not os.path.isdir(path):
		continue
		
	if path != "versao_base":
		continue
	
	try:
		subprocess.run(["make", "benchmark"], cwd=path, capture_output=True, text=True, check=True)
	except subprocess.CalledProcessError as e:
		print(f"'make benchmark' falhou para {path}: {e.stderr}")
		continue
	
	for N in Ns:
		for group in ["L3", "L2CACHE", "FLOPS_DP", "FLOPS_AVX"]:
			print(f"Medindo {group} para a {path}")
			res = subprocess.run(["likwid-perfctr", "-O", "-C", "3", "-g", group, "-m", "./broyden"], input=f"{N} -1 0 25", cwd=path, capture_output=True, text=True, check=True)
			
			
			print(res.stdout)
			input()
