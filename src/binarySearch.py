import multiprocessing
import networkx as nx
import random
from concurrent.futures import ThreadPoolExecutor

import definition
import simulation


GRAPH_NUM = 16


def binary_search(design: definition.BLIF_design):
	initial_config = [False] * len(design.lut62)
	e_th = design.e_th
	core_used = int(multiprocessing.cpu_count() * 2 / 4)

	# First Iteration: Initial Multi-thread Arguments
	jobs = []
	jobs_ID = []
	# Iterate over each design option.
	for i in range(len(initial_config)):
		# Generate a new design by changing the i-th design option.
		new_config = initial_config[:]
		new_config[i] = True
		# Add new config to the job list
		jobs.append(new_config)
		jobs_ID.append(i)
	# Simulate the accuracy of the new design in parallel
	with ThreadPoolExecutor(max_workers=core_used) as executor:
		results = executor.map(simulation.accuracy_sim, [design]*len(jobs), jobs, jobs_ID)
	# Sort the designs by accuracy.
	results_lst = list(results)
	sorted_results = sorted(results_lst, key=lambda x: x[1], reverse=False)
	# Update best implementation
	simulation.find_best_implement(results_lst, design)

	# Binary Search Process
	left, right = 0, len(sorted_results) - 1
	iteration_ID = 1
	best_design = None
	while True:
		mid = (left + right) // 2
		# Cannot be further searched
		if mid == left:
			break
		pair_ID = [rst[2] for rst in sorted_results[:mid + 1]]
		# Create a graph
		G = nx.Graph()
		edges = []
		for ID in pair_ID:
			edges.append((design.lut62[ID].lut1, design.lut62[ID].lut2, ID))
		for u, v, tag in edges:
			G.add_edge(u, v, tag=tag)
		# For each iteration, find several maximum matching
		matchings = []
		edge_to_remove = None
		for _ in range(GRAPH_NUM):
			# Find one maximum matching
			G_copy = G.copy()
			if edge_to_remove is not None:
				G_copy.remove_edge(*edge_to_remove)
			matching = nx.algorithms.matching.max_weight_matching(G_copy, maxcardinality=True)
			if len(matching) == 0:
				break
			edge_to_remove = random.choice(list(matching))
			if matching not in matchings:
				matchings.append(matching)
		# Initial Multi-thread Arguments
		jobs = []
		jobs_ID = []
		for matching in matchings:
			# Generate new config according to the maximum matching
			new_config = initial_config[:]
			for u, v in matching:
				new_config[G[u][v]['tag']] = True
			# Add new config to the job list
			jobs.append(new_config)
			jobs_ID.append(iteration_ID)
			iteration_ID += 1
		# Simulate all maximum matchings
		with ThreadPoolExecutor(max_workers=core_used) as executor:
			results = executor.map(simulation.accuracy_sim, [design]*len(jobs), jobs, jobs_ID)
		# Sort the designs by accuracy.
		new_rst_list = list(results)
		new_sorted_results = sorted(new_rst_list, key=lambda x: x[1], reverse=False)
		if len(new_sorted_results) == 0:
			break
		if new_sorted_results[0][1] <= e_th:
			best_design = new_sorted_results[0]
			left = mid
		else:
			right = mid

	# If no improvement, return the original design
	if best_design is None:
		return 0, simulation.synthesis(design, initial_config)
	return best_design[1], simulation.synthesis(design, best_design[0])
