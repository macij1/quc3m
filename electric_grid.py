"""
This script demonstrates the full workflow for finding an
Optimal Corrective Action plan for a power grid.

It follows this process:
1. Model the physical grid (the "initial graph").
   **This version uses a more realistic 14-bus system.**
2. Define a set of candidate actions.
3. Simulate pairs of actions to find conflicts.
4. Build a "conflict graph" based on these incompatibilities.
5. **Plot both the physical grid and the conflict graph.**

"""

import networkx as nx
import matplotlib.pyplot as plt
import itertools
from dataclasses import dataclass, field
from typing import List, Dict, Any

# --- Step 1: Model the Physical Grid ---

def create_physical_grid():
    """
    Creates a 'networkx' graph model of a power grid.
    This is based on the standard IEEE 14-BUS TEST CASE,
    representing a more complex and realistic network.
    """
    grid = nx.Graph()

    # Node Types:
    # - SLACK: The "swing" bus, main generator that balances the system.
    # - PV: "Generator" bus (voltage-controlled).
    # - PQ: "Load" bus (consumes power).
    # - SUB: A simple substation, just for connecting lines.

    # Add 14 buses (nodes)
    grid.add_node("B1_SLACK", type="generator", output_mw=232.4)
    grid.add_node("B2_PV", type="generator", output_mw=40.0)
    grid.add_node("B3_PV", type="generator", output_mw=0.0) # Synchronous condenser
    grid.add_node("B4_SUB", type="substation")
    grid.add_node("B5_SUB", type="substation")
    grid.add_node("B6_PV", type="generator", output_mw=0.0) # Synchronous condenser
    grid.add_node("B7_SUB", type="substation")
    grid.add_node("B8_PV", type="generator", output_mw=0.0) # Synchronous condenser
    grid.add_node("B9_SUB", type="substation")
    grid.add_node("B10_SUB", type="substation")
    grid.add_node("B11_SUB", type="substation")
    grid.add_node("B12_SUB", type="substation")
    grid.add_node("B13_SUB", type="substation")
    grid.add_node("B14_SUB", type="substation")
    
    # Add loads to the buses (in a real model, these would be
    # separate nodes, but we can attach them as attributes)
    grid.nodes["B2_PV"]["demand_mw"] = 21.7
    grid.nodes["B3_PV"]["demand_mw"] = 94.2
    grid.nodes["B4_SUB"]["demand_mw"] = 47.8
    grid.nodes["B5_SUB"]["demand_mw"] = 7.6
    grid.nodes["B6_PV"]["demand_mw"] = 11.2
    grid.nodes["B9_SUB"]["demand_mw"] = 29.5
    grid.nodes["B10_SUB"]["demand_mw"] = 9.0
    grid.nodes["B11_SUB"]["demand_mw"] = 3.5
    grid.nodes["B12_SUB"]["demand_mw"] = 6.1
    grid.nodes["B13_SUB"]["demand_mw"] = 13.5
    grid.nodes["B14_SUB"]["demand_mw"] = 14.9
    
    # Add 20 transmission lines & transformers (edges)
    # (Simplified: not showing real impedance, just capacity)
    grid.add_edge("B1_SLACK", "B2_PV", line_id="L1-2", capacity_mw=120)
    grid.add_edge("B1_SLACK", "B5_SUB", line_id="L1-5", capacity_mw=120)
    grid.add_edge("B2_PV", "B3_PV", line_id="L2-3", capacity_mw=120)
    grid.add_edge("B2_PV", "B4_SUB", line_id="L2-4", capacity_mw=120)
    grid.add_edge("B2_PV", "B5_SUB", line_id="L2-5", capacity_mw=120)
    grid.add_edge("B3_PV", "B4_SUB", line_id="L3-4", capacity_mw=120)
    grid.add_edge("B4_SUB", "B5_SUB", line_id="L4-5", capacity_mw=60) # Lower capacity line
    grid.add_edge("B4_SUB", "B7_SUB", line_gide="L4-7", capacity_mw=120)
    grid.add_edge("B4_SUB", "B9_SUB", line_id="L4-9", capacity_mw=120)
    grid.add_edge("B5_SUB", "B6_PV", line_id="L5-6", capacity_mw=120)
    grid.add_edge("B6_PV", "B11_SUB", line_id="L6-11", capacity_mw=60)
    grid.add_edge("B6_PV", "B12_SUB", line_id="L6-12", capacity_mw=60)
    grid.add_edge("B6_PV", "B13_SUB", line_id="L6-13", capacity_mw=120)
    grid.add_edge("B7_SUB", "B8_PV", line_id="L7-8", capacity_mw=60)
    grid.add_edge("B7_SUB", "B9_SUB", line_id="L7-9", capacity_mw=60)
    grid.add_edge("B9_SUB", "B10_SUB", line_id="L9-10", capacity_mw=60)
    grid.add_edge("B9_SUB", "B14_SUB", line_id="L9-14", capacity_mw=60)
    grid.add_edge("B10_SUB", "B11_SUB", line_id="L10-11", capacity_mw=60)
    grid.add_edge("B12_SUB", "B13_SUB", line_id="L12-13", capacity_mw=60)
    grid.add_edge("B13_SUB", "B14_SUB", line_id="L13-14", capacity_mw=60)
    
    # We assume this grid is in a strained "N-1" state.
    # A real simulation would calculate flows, but for this
    # model, just having the topology is enough.

    print(f"Physical Grid (IEEE 14-Bus): {grid.number_of_nodes()} nodes and {grid.number_of_edges()} edges.")
    return grid

# --- Step 2: Define Candidate Actions ---

@dataclass
class CorrectiveAction:
    action_id: str
    description: str
    value: float  # The "goodness" of this action (its weight)
    params: Dict[str, Any] = field(default_factory=dict)

def get_candidate_actions() -> List[CorrectiveAction]:
    """
    Returns a list of potential actions to solve the grid instability.
    These are just examples; a real system would generate them
    based on the 14-bus topology.
    """
    return [
        CorrectiveAction(
            action_id="A1_REROUTE_G1",
            description="Ramp down B1_SLACK by 50MW",
            value=8.0,
            params={"gen": "B1_SLACK", "change": -50}
        ),
        CorrectiveAction(
            action_id="A2_REROUTE_G2",
            description="Ramp up B2_PV by 30MW",
            value=7.0,
            params={"gen": "B2_PV", "change": 30}
        ),
        CorrectiveAction(
            action_id="A3_SHED_LOAD_B9",
            description="Shed 20MW of non-critical load at B9_SUB",
            value=4.0, # Lower value
            params={"load": "B9_SUB", "change": -20}
        ),
        CorrectiveAction(
            action_id="A4_TOPOLOGY_B4",
            description="Split busbar at Substation B4_SUB",
            value=9.0, # High value
            params={"substation": "B4_SUB", "action": "split_bus"}
        ),
        CorrectiveAction(
            action_id="A5_MUTUALLY_EXCLUSIVE_G1",
            description="Boost B1_SLACK stability (exclusive with A1)",
            value=6.0,
            params={"gen": "B1_SLACK", "action": "boost"}
        ),
        CorrectiveAction(
            action_id="A6_SHED_LOAD_B14",
            description="Shed 10MW of non-critical load at B14_SUB",
            value=2.0,
            params={"load": "B14_SUB", "change": -10}
        )
    ]

# --- Step 3: The Simulator (The "Hard Part") ---

def run_simulation(
    physical_grid: nx.Graph,
    actions: List[CorrectiveAction]
) -> str:
    """
    A **DUMMY SIMULATOR** that predicts the grid state after
    applying a set of actions.

    The logic here is hard-coded for this demo, but in a real
    system it would be running a power flow simulation on the
    14-bus grid model.
    """
    action_ids = {a.action_id for a in actions}

    # CONFLICT 1: Mutual Exclusion (e.g., can't ramp down AND boost B1)
    if "A1_REROUTE_G1" in action_ids and "A5_MUTUALLY_EXCLUSIVE_G1" in action_ids:
        return "UNSTABLE_MUTUAL_EXCLUSION_B1"

    # CONFLICT 2: Negative Interaction (e.g., two actions overload L4-9)
    if "A1_REROUTE_G1" in action_ids and "A2_REROUTE_G2" in action_ids:
        return "OVERLOAD_L4-9_UNSTABLE"

    # CONFLICT 3: Resource Contention (e.g., two actions try to use B4)
    if "A2_REROUTE_G2" in action_ids and "A4_TOPOLOGY_B4" in action_ids:
        return "UNSTABLE_RESOURCE_CONFLICT_B4"
        
    # CONFLICT 4: New conflict for our new actions
    if "A3_SHED_LOAD_B9" in action_ids and "A4_TOPOLOGY_B4" in action_ids:
        # Let's pretend shedding load at B9 while splitting B4
        # isolates B10 and causes a voltage collapse.
        return "VOLTAGE_COLLAPSE_B10"
    
    if "A3_SHED_LOAD_B9" in action_ids and "A2_REROUTE_G2" in action_ids:
            return "VOLTAGE_COLLAPSE_B9"
    
    if "A6_SHED_LOAD_B14" in action_ids and "A2_REROUTE_G2" in action_ids:
            return "VOLTAGE_COLLAPSE_B9"


    # If no conflicts are triggered, we assume the state is stable
    return "STABLE"

# --- Step 4: Build the Conflict Graph ---

def build_conflict_graph(
    candidate_actions: List[CorrectiveAction],
    physical_grid: nx.Graph
) -> nx.Graph:
    """
    Builds the conflict graph by simulating all pairs of actions.
    """
    print("Building conflict graph...")
    conflict_graph = nx.Graph()

    # 1. Add all candidate actions as NODES.
    for action in candidate_actions:
        conflict_graph.add_node(action.action_id, value=action.value)
        print(f"  + Added Node: {action.action_id} (Value: {action.value})")

    # 2. Test all pairs of actions to find EDGES (conflicts).
    print("Testing action pairs for conflicts...")
    for action_a, action_b in itertools.combinations(candidate_actions, 2):
        state = run_simulation(physical_grid, [action_a, action_b])

        if state != "STABLE":
            print(f"  ! CONFLICT: {action_a.action_id} <> {action_b.action_id} (Reason: {state})")
            
            conflict_graph.add_edge(
                action_a.action_id,
                action_b.action_id,
                reason=state 
            )
        
    return conflict_graph

def plot_physical_grid(grid: nx.Graph):
    """
    Uses Matplotlib to draw the physical power grid.
    """
    plt.figure(1, figsize=(12, 10)) 
    plt.title("Physical Grid (IEEE 14-Bus Model)")

    # 1. Create a color map based on node type
    node_types = nx.get_node_attributes(grid, 'type')
    color_map = []
    for node in grid.nodes():
        if node_types.get(node) == 'generator':
            color_map.append('skyblue')
        elif grid.nodes[node].get("demand_mw", 0) > 0:
            color_map.append('salmon')
        else:
            color_map.append('lightgray')

    # 2. Define layout
    pos = nx.kamada_kawai_layout(grid) 

    # 3. Draw nodes, edges, and labels
    nx.draw_networkx_nodes(grid, pos, node_color=color_map, node_size=2000)
    nx.draw_networkx_edges(grid, pos, width=1.5, alpha=0.7)
    nx.draw_networkx_labels(grid, pos, font_size=8, font_weight='bold')

    # 4. Draw edge labels (line IDs)
    edge_labels = nx.get_edge_attributes(grid, 'line_id')
    
    plt.axis('off')


def plot_conflict_graph(graph: nx.Graph):
    """
    Uses Matplotlib to draw the abstract conflict graph.
    """
    plt.figure(2, figsize=(12, 10))
    plt.title("Conflict Graph (Optimization Problem)")

    node_values = nx.get_node_attributes(graph, 'value')
    labels = {}
    for node_id, value in node_values.items():
        short_id = node_id.split('_')[0]
        labels[node_id] = f"{short_id}\n(Val: {value})"

    pos = nx.spring_layout(graph, seed=42) # Ensure consistent layout

    nx.draw_networkx(
        graph,
        pos,
        labels=labels,
        with_labels=True,
        node_color='coral',
        node_size=3000,
        font_size=10,
        font_weight='bold'
    )
    
    edge_labels = {}
    for u, v in graph.edges():
        short_u = u.split('_')[0]
        short_v = v.split('_')[0]
        edge_labels[(u, v)] = f"$\lambda_{{{short_u},{short_v}}}$"

    label_offsets = {
        ("A1_REROUTE_G1", "A2_REROUTE_G2"): (-0.02, 0.05), # Push A1-A2 label up a bit
        ("A2_REROUTE_G2", "A6_SHED_LOAD_B14"): (0.02, 0.05), # Push A2-A6 label up and right
        ("A2_REROUTE_G2", "A3_SHED_LOAD_B9"): (-0.03, -0.05), # Push A2-A3 label down and left
        ("A2_REROUTE_G2", "A4_TOPOLOGY_B4"): (0.05, -0.01), # Push A2-A4 label right
        ("A3_SHED_LOAD_B9", "A4_TOPOLOGY_B4"): (0.00, -0.05), # Push A3-A4 label down
        ("A1_REROUTE_G1", "A5_MUTUALLY_EXCLUSIVE_G1"): (0.02, 0.05), # Push A1-A5 label up and right
    }

    # Draw edge labels, applying manual offsets
    for (u, v), label in edge_labels.items():
        # Get the midpoint of the edge
        x_mid = (pos[u][0] + pos[v][0]) / 2
        y_mid = (pos[u][1] + pos[v][1]) / 2

        # Apply custom offset if available
        dx, dy = label_offsets.get((u, v), (0, 0))
        # Also check for (v,u) in case the edge was added in reverse order
        dx_rev, dy_rev = label_offsets.get((v, u), (0, 0))
        dx += dx_rev
        dy += dy_rev

        # Add the label as an annotation for fine-grained control
        plt.annotate(
            label,
            xy=(x_mid + dx, y_mid + dy),
            xytext=(x_mid + dx, y_mid + dy), # Same as xy for text position
            fontsize=12,
            color='black',
            ha='center', # Horizontal alignment
            va='center', # Vertical alignment
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.6) # Light background for readability
        )
    
    plt.axis('off')
    plt.show() # Display the plot

# --- Main Execution ---

if __name__ == "__main__":
    print("--- Power Grid Corrective Action Optimizer ---")

    # 1. Setup the physical system
    print("\n[Phase 1: Setup]")
    grid = create_physical_grid()
    actions = get_candidate_actions()
    action_map = {a.action_id: a for a in actions}

    # 2. Build the optimization model
    print("\n[Phase 2: Building Conflict Graph]")
    conflict_model = build_conflict_graph(actions, grid)

    # 3. Plot the graphs
    print("\n[Phase 3: Plotting Graphs]")
    try:
        plot_physical_grid(grid)
        plot_conflict_graph(conflict_model)
        
        print("... Displaying plot windows. Close the plot windows to exit.")
        plt.show()

    except ImportError:
        print("\n[Plotting Skipped]")
        print("Please install 'matplotlib' to visualize the graphs.")
        print("Run: pip install matplotlib")

    print("\n[Phase 4: Conflict Model Ready]")
    print(f"Nodes (Actions): {list(conflict_model.nodes(data=True))}")
    print(f"Edges (Conflicts): {list(conflict_model.edges())}")
    print("\nYou can now pass the 'conflict_model' object to your MWIS solver.")