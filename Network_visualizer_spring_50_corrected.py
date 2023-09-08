
import networkx as nx
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description="Generate a protein-protein similarity network and extract all the nodes a certain number of edges away from any member of a set of interest.")
parser.add_argument('--maxdist', type=int, default=2, help='Maximum distance from any member of a set of interest.')
parser.add_argument('--input_fasta', type=str, required=True, help='Input protein FASTA file.')
parser.add_argument('--output_fasta', type=str, required=True, help='Output subset FASTA file.')
parser.add_argument('--png', type=str, default='network_diagram.png', help='Output PNG file.')
args = parser.parse_args()

# Read the data
with open('alignments.tsv', 'r') as file:
    data = file.readlines()

# Create a new directed graph
G = nx.DiGraph()

# Store colors of the nodes
node_colors = {}

for line in data:
    # Split the line into columns
    columns = line.split("\t")

    # Get the query protein, subject protein, and percent identity
    query_protein = columns[0]
    subject_protein = columns[1]
    percent_identity = float(columns[2])

    # Assign node colors based on source of the sequence
    for protein in (query_protein, subject_protein):
        if 'Giant_Proteins' in protein:
            node_colors[protein] = 'purple'
        elif 'Archaea' in protein:
            node_colors[protein] = 'green'
        elif 'Uniprot' in protein:
            node_colors[protein] = 'blue'
        elif 'Omni' in protein:
            node_colors[protein] = 'red'
        elif 'Bacteria__' in protein:
            node_colors[protein] = 'gold'
        else:
            #Pink is bad! You should not see pink!
            node_colors[protein] = 'pink'

    # Add an edge between the query and subject proteins with weight >= 40%
    if query_protein != subject_protein and percent_identity >= 40:
        G.add_edge(query_protein, subject_protein, weight=percent_identity)

# Define a set for nodes with distance less than or equal to maxdist from any purple node
close_nodes = set()

# Define a function to find all nodes within maxdist from a given node
def get_nodes_within_maxdist(G, start_node, maxdist):
    nodes_within_maxdist = set()
    nodes_to_visit = [(start_node, 0)]
    while nodes_to_visit:
        current_node, current_dist = nodes_to_visit.pop()
        if current_node not in nodes_within_maxdist and current_dist <= maxdist:
            nodes_within_maxdist.add(current_node)
            nodes_to_visit.extend((n, current_dist + 1) for n in nx.neighbors(G, current_node))
    return nodes_within_maxdist

# For each purple node, we look for nodes within maxdist and add them to the set
for node in G.nodes():
    if node_colors[node] == 'purple':
        close_nodes.update(get_nodes_within_maxdist(G, node, args.maxdist))

# Now we filter the graph to consider only the nodes in close_nodes
G = G.subgraph(close_nodes)

# Read the sequences from the input FASTA file, keeping only the first sequence for each identifier
sequences = {}
for record in SeqIO.parse(args.input_fasta, "fasta"):
    if record.id not in sequences:
        sequences[record.id] = record

# Filter the sequences to include only those that are in the graph
filtered_sequences = [seq for seq_id, seq in sequences.items() if seq_id in G.nodes]

# Write the filtered sequences to the output FASTA file
SeqIO.write(filtered_sequences, args.output_fasta, "fasta")

# Define size of the figure (you may need to adjust this)
size = 30
plt.figure(figsize=(size, size))

# Set the title of the plot
plt.title('C39 Peptidase Identity Network', fontsize=40)

# Define size of the nodes and width of the edges (you may need to adjust these)
node_size = 100
edge_width = 0.1

# Define position of the nodes using the spring layout
pos = nx.spring_layout(G, k=2)

# Draw the nodes and edges, but not the nodes of interest
other_nodes = [node for node in G.nodes() if node_colors[node] != 'purple']
other_colors = [node_colors[node] for node in other_nodes]
nx.draw_networkx_nodes(G, pos, nodelist=other_nodes, node_color=other_colors, node_size=node_size)

nx.draw_networkx_edges(G, pos, width=edge_width)

# Draw the nodes of interest
interest_nodes = [node for node in G.nodes() if node_colors[node] == 'purple']
nx.draw_networkx_nodes(G, pos, nodelist=interest_nodes, node_color='purple', node_size=node_size)

# Create a legend
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label='Giant Protein', markersize=15, markerfacecolor='purple'),
                   plt.Line2D([0], [0], marker='o', color='w', label='Omnitrophota', markersize=15, markerfacecolor='red'),
                   plt.Line2D([0], [0], marker='o', color='w', label='Archaea', markersize=15, markerfacecolor='green'),
                   plt.Line2D([0], [0], marker='o', color='w', label='Uniprot', markersize=15, markerfacecolor='blue'),
                   plt.Line2D([0], [0], marker='o', color='w', label='GTDB Bacteria', markersize=15, markerfacecolor='gold')]
plt.legend(handles=legend_elements, loc='upper left')

# Save the plot as a PNG file
plt.savefig(args.png)
