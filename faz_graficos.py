import matplotlib.pyplot as plt
import subprocess
import os

plt.style.use('ggplot')


Ns = [32, 64, 128, 256, 512, 1000, 2000, 4000, 8000, 9000, 10000, 20000]
groups = ["L3", "L2CACHE", "FLOPS_AVX", "FLOPS_DP"]
targets = ["L3 bandwidth [MBytes/s]", "L2 miss ratio", "Packed DP [MFLOP/s]", "DP [MFLOP/s]"]
regions = ["Total", "Jacobiana", "Sistema Linear"]


graficos = {}
for group in groups + ["tempo"]:
	graficos[group] = {}
	for region in regions:
		graficos[group][region] = {}

for path in os.listdir("./"):
	if not path.startswith("versao") or not os.path.isdir(path):
		continue
	
	try:
		subprocess.run(["make", "benchmark"], cwd=path, capture_output=True, text=True, check=True)
	except subprocess.CalledProcessError as e:
		print(f"'make benchmark' falhou para {path}: {e.stderr}")
		continue
	
	for group in groups + ["tempo"]:
		for region in regions:
			graficos[group][region][path] = []
	
	for N in Ns:

		# pega a média dos tempos para todas as execuções (uma por grupo)
		tempos = { "Total": 0, "Jacobiana": 0, "Sistema Linear": 0 }
		vezes_calculado = { "Total": 0, "Jacobiana": 0, "Sistema Linear": 0 }
			
		for group, target in zip(groups, targets):
				
			# print(f"Medindo {group} para a {path}")
			
			res = subprocess.run(["likwid-perfctr", "-O", "-C", "3", "-g", group, "-m", "./broyden"], input=f"{N} -1 0 25", cwd=path, capture_output=True, text=True, check=True)
			
			region = None
			for line in res.stdout.split("\n"):
				
				if line.startswith("# Tempo "):
					regiao = line.removeprefix("# Tempo ").split(": ")[0]
					regiao = regiao if regiao != "SL" else "Sistema Linear" # arruma o nome pra ficar que nem nos gráficos
					tempo = float(line.removeprefix("# Tempo ").split(": ")[1])
					tempos[regiao] += tempo
					vezes_calculado[regiao] += 1
					
					continue
				
				parts = line.split(",")
				if len(parts) < 2:
					continue
					
				if parts[0] == "TABLE":
					
					if parts[2].endswith("Raw"): # pega só a métrica não o dado "cru"
						region = None
					else:
						region = parts[1].removeprefix("Region ")
						continue
						
				elif parts[0] == "STRUCT":
					region = None
				
				if region is None:
					continue
				
				if parts[0] == target or parts[0] == target[0]:
					graficos[group][region][path].append(float(parts[1]))
					# print(f"{region}: {parts[0]} -> {parts[1]}")
			
			
			# input()
			
		for key in tempos.keys():
			if vezes_calculado[key] == 0:
				graficos["tempo"][key][path].append(0)
			else:
				graficos["tempo"][key][path].append(tempos[key] / vezes_calculado[key])
	
	for group in groups + ["tempo"]:
		for region in regions:
			if len(graficos[group][region][path]) == 0:
				del graficos[group][region][path]



for group, target in zip(groups + ["tempo"], targets + ["Tempo"]):
	for region in regions:
		for k, v in graficos[group][region].items():
			plt.plot(Ns, v, label=k)
		
		plt.title(f"{region} - {group}")
		plt.legend()
		plt.xscale("log")
		plt.xlabel("N (log)")
		plt.ylabel(target)
		if target == "Tempo":
			plt.yscale("log")
			plt.ylabel(target + " (log)")
		
		os.makedirs(f"resultados/{region}/", exist_ok=True)
		plt.savefig(f"resultados/{region}/{group}.png", bbox_inches="tight")
		plt.close()





















