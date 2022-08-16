
# CONDA ENV: 

from graphviz import Graph, Digraph, Source
import json
import os

#def create_digraph_system(model):
#    dot = Digraph(comment='System')
#    for edge in model['edges']:
#        dot.edge(edge[0], edge[1])
#    dot.render('system.gv', view=True)
    
models_path = 'gj_json_files'
model_template = os.path.join(
    models_path, 'gj_{region}_section.json')

fill_colors = {
    'storage': 'blue',
    'catchment': 'lightblue',
    'input': 'lightblue',
    'output': 'pink',
}

font_colors = {
    'storage': 'white'
}

node_shapes = {
    'storage': 'rect'
}

font_sizes = {
    'storage': 14
}

subgraph_attrs = {
    'style': 'filled',
    'color': '#f0f0f0'
}


def create_digraph(nodes, edges, node_lookup, up_nodes, down_nodes, 
    focal_system, filename, focal_names=None, output_format='png'):

    print('creating ' + focal_system)
    all_nodes = []
    all_edges = []    
    focal_nodes = None
    all_focal_names = []
    
    subgraphs = {}
    
    parent = Digraph(name=focal_system, format=output_format)
    
    if focal_names:
        focal_nodes = [node for node in nodes if node['name'] in focal_names]
        all_focal_names = focal_names[:]
        for focal_name in focal_names:
    #         print(focal_name)
            all_focal_names.extend([n['name'] for n in up_nodes.get(focal_name, []) if n['name'] not in all_focal_names])
            all_focal_names.extend([n['name'] for n in down_nodes.get(focal_name, []) if n['name'] not in all_focal_names])
    #         up_node_names = [n['name'] for n in up_nodes.get(node['name'], [])]
    #         down_node_names = [n['name'] for n in down_nodes.get(node['name'], [])]

    def add_node(node, subgraph, color='black'):
        if node['name'] in all_nodes:
            return

        all_nodes.append(node['name'])
        label = node['name'].replace('_', ' ')

        fontsize = str(font_sizes.get(node['type'], 10))
        fillcolor = fill_colors.get(node['type'], 'lightgray')
        fontcolor = font_colors.get(node['type'], 'black')
        shape = node_shapes.get(node['type'], 'ellipse')

        if node['type'] == 'storage':
            if node['name'].find('_') >= 0:
                label = label.title()
            else: # this is a WTP
                fillcolor = 'grey'
                fontcolor = 'black'

        if node['name'].find('evap') == 0:
            fillcolor = 'white'
        elif node['name'].find('ifr') == 0:
            fillcolor = '#99ffcc'

        if node['name'].find('hetch_hetchy') != 0:
            label = label.replace('hetch hetchy', '\nhetch hetchy')
        label = label.replace('reservoir', 'reservoir\n')
        
        graphprops = dict(
            shape=shape,
            color=color,
            style='filled',
            fontcolor=fontcolor,
            fillcolor=fillcolor,
            label=label,
            fontsize=fontsize
        )

        if all_focal_names and node['name'] not in all_focal_names:
            for prop in ['style', 'fillcolor', 'fontcolor', 'color']:
                graphprops.pop(prop, None)
#             graphprops['fillcolor'] = 'lightgrey'

        subgraph.node(node['name'], **graphprops)
        
    for node in nodes:
        
        region = node.get('region')
        if region not in subgraphs:
            # note: subgraph name must start with "cluster_"
            subgraphs[region] = Digraph(name='cluster_{}'.format(region), graph_attr=subgraph_attrs)
        subgraph = subgraphs[region]
        
        _down_nodes = down_nodes.get(node['name'], [])
        _up_nodes = up_nodes.get(node['name'], [])
        if _up_nodes or _down_nodes:
            add_node(node, subgraph)

        for down_node in _down_nodes:
            edge = (node['name'], down_node['name'])
            if edge in all_edges:
                continue
            all_edges.append(edge)
            add_node(down_node, subgraph)
            color = 'black'
            if all_focal_names and not (edge[0] in all_focal_names and edge[1] in all_focal_names):
                color = 'grey'
            graph = subgraph if node['region'] == down_node['region'] else parent
            graph.edge(edge[0], edge[1], color=color)

        for up_node in _up_nodes:
            edge = (up_node['name'], node['name'])
            if edge in all_edges:
                continue
            all_edges.append(edge)
            add_node(up_node, subgraph)

            color = 'black'
            if all_focal_names and not (edge[0] in all_focal_names and edge[1] in all_focal_names):
                color = 'grey'
            
            graph = subgraph if node['region'] == down_node['region'] else parent
            graph.edge(edge[0], edge[1], color=color)
    
    for child in subgraphs.values():
        parent.subgraph(child)
    parent.render(filename, view=False)
    del parent

# setup model
nodes = []
edges = []
node_lookup = {}
up_nodes = {}
down_nodes = {}

exclude = ['dummy_output']

for region in ['allende', 'alzate','ameche','arandas','corrales','markazuza','modulo2','ocampo','outlet','pericos','purissima','ramirez','salamanca','salvatierra','solis','tepetitlan','tepuxtepec','turbio','yuriria']:
    with open(model_template.format(region=region)) as f:
        model = json.load(f)
        edges.extend([e for e in model['edges'] if not(e[0] in exclude or e[1] in exclude)])
        for node in model['nodes']:
            if node['name'] in exclude:
                continue
            node['region'] = region
            nodes.append(node)
            
print('Total nodes: {}'.format(len(nodes)))

node_lookup = {node['name']: node for node in nodes}

for edge in edges:
    up_node = node_lookup.get(edge[0])
    down_node = node_lookup.get(edge[1])
    if up_node and down_node:
        up_nodes[down_node['name']] = up_nodes.get(down_node['name'], []) + [up_node]
        down_nodes[up_node['name']] = down_nodes.get(up_node['name'], []) + [down_node]
    
focal_systems = {
    'hetch_hetchy': ['hetch_hetchy_reservoir'],
    'cherry_eleanor': ['cherry_reservoir', 'eleanor_reservoir'],
    'don_pedro': ['don_pedro_reservoir'],
    'east_bay': ['san_antonio_reservoir', 'calaveras_reservoir', 'recapture_pit'],
    'pilarcitos': ['pilarcitos_reservoir'],
    'east_peninsula': ['san_andreas_reservoir', 'crystal_springs_reservoir'],
    'peninsula': ['san_andreas_reservoir', 'crystal_springs_reservoir', 'pilarcitos_reservoir']
}

OUTPUT_FORMAT_LIST = ['png', 'pdf']
for OUTPUT_FORMAT in OUTPUT_FORMAT_LIST:
    create_digraph(nodes, edges, node_lookup, up_nodes, down_nodes,
        'system', 'schematics/system.gv', output_format=OUTPUT_FORMAT)
    for focal_system in focal_systems:
    #     create_digraph(focal_system, 'schematics/{}.gv'.format(focal_system))
        create_digraph(nodes, edges, node_lookup, up_nodes, down_nodes, focal_system,
                       'schematics/{}.gv'.format(focal_system), 
                       focal_names=focal_systems[focal_system], 
                       output_format=OUTPUT_FORMAT)