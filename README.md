# QUC3M – Quantum Virtual Screening (Hackathon Project)

This repository contains our solution for the **Quantum Madrid × UPM Hackathon**, where the challenge was to detect molecular similarity using the **Fujitsu Digital Annealer**.  
Our team, **QUC3M**, received the **Honorary Mention for Innovative Applicability** for demonstrating that the same QUBO-based pipeline generalizes far beyond chemistry.

## About the Project

The hackathon problem asked us to formalize molecular similarity as a **Maximum Weighted Independent Set (MWIS)** problem on a conflict graph.  
We implemented:

- RDKit-based **molecular graph construction**
- Automatic extraction of **pharmacophoric features**
- Pairwise **conflict graph** generation
- **MWIS-via-QUBO formulation** suitable for quantum or classical solvers
- Experiments on multiple QUBO SDKs (`dimod`, `pyqubo`, `dadk`)

We later showed the same pipeline applies to a completely different domain:  
**stabilizing unstable electric grids using action-conflict graphs**, which earned us the Innovative Applicability award.

---

## 📓 Hackathon Problem & Solution Notebooks

➡️ **The Jupyter notebooks in this repo include the full hackathon problem and our final solution.**  
They document the reasoning, the QUBO formulation, intermediate experiments, and final results.

If you're here to understand the challenge or reproduce our solution, **start with the notebooks**.

---

## Repository Structure

- **RetoMolecularComparison/modules/**  
  Core implementation: molecular graphs, conflict graphs, feature logic, node definitions.
  
- **dadk-documentation/**  
  Minimal docs + benchmarking assets.

- **Documentation Viewer.ipynb**  
  Helper to spin up a tiny local documentation server.

- **dadk_benchmark.py**  
  Benchmark script comparing QUBO SDKs on a TSP instance.

---

## Installation

Install dependencies (RDKit + scientific packages):

`pip install rdkit-pypi networkx matplotlib dimod pyqubo memory_profiler`

--- 

## Minimal Usage Example

```python
from modules.molecular_graph import MolecularGraphDA
from modules.conflict_graph import ConflictGraph
from rdkit import Chem

mol1 = MolecularGraphDA(Chem.MolFromSmiles("CCO"))
mol2 = MolecularGraphDA(Chem.MolFromSmiles("CCN"))

cg = ConflictGraph(mol1, mol2)
print("Nodes:", len(cg.graph.nodes))
print("Edges:", len(cg.graph.edges))
```
