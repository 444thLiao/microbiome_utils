from __future__ import division
import json
import math
import numpy as np

def visualize(complex, color_function="", path_html="mapper_visualization_output.html", title="My Data",
              graph_link_distance=50, graph_gravity=0.1, graph_charge=-50, custom_tooltips=None, width_html=0,
              height_html=0, show_tooltips=True, show_title=True, show_meta=True):
    # Turns the dictionary 'complex' in a html file with d3.js
    #
    # Input:      complex. Dictionary (output from calling .map())
    # Output:      a HTML page saved as a file in 'path_html'.
    #
    # parameters
    # ----------
    # color_function    	string. Not fully implemented. Default: "" (distance to origin)
    # path_html        		file path as string. Where to save the HTML page.
    # title          		string. HTML page document title and first heading.
    # graph_link_distance  	int. Edge length.
    # graph_gravity     	float. "Gravity" to center of layout.
    # graph_charge      	int. charge between nodes.
    # custom_tooltips   	None or Numpy Array. You could use "y"-label array for this.
    # width_html        	int. Width of canvas. Default: 0 (full width)
    # height_html       	int. Height of canvas. Default: 0 (full height)
    # show_tooltips     	bool. default:True
    # show_title      		bool. default:True
    # show_meta        		bool. default:True

    ####################################################################################################################
    # Format JSON for D3 graph
    json_s = complex
    #json_s["nodes"] = []
    #json_s["links"] = []
    #k2e = {}  # a key to incremental int dict, used for id's when linking

    node_size_scale_factor = 100
    node_sizes = [len(it['samples']) for it in json_s['nodes']]
    node_size_max, node_size_min = max(node_sizes), min(node_sizes)

    link_weight = [it['weight'] for it in json_s['links']]
    link_weight_max, link_weight_min = max(link_weight), min(link_weight)

    for it in json_s['nodes']:
        it['name'] = str(it['id'])
        it['group'] = math.sqrt(((len(it['samples'])-node_size_min)*1.0/node_size_max) * node_size_scale_factor + 50)
        if '-N' in it['name']:
            it['color'] = 9
        else:
            it['color'] = 21
    for it in json_s['links']:
        it['weight'] = ((it['weight']-link_weight_min)*1.0/link_weight_max) * 20 + 2

    new_json = {}
    new_json['nodes'] = json_s['nodes']
    new_json['links'] = json_s['links']
    json_s = new_json
        #it['color'] = np.mean(np.asarray([custom_tooltips[idx] for idx in it['samples']]))

    # Width and height of graph in HTML output
    if width_html == 0:
        width_css = "100%"
        width_js = 'document.getElementById("holder").offsetWidth-20'
    else:
        width_css = "%spx" % width_html
        width_js = "%s" % width_html
    if height_html == 0:
        height_css = "100%"
        height_js = 'document.getElementById("holder").offsetHeight-20'
    else:
        height_css = "%spx" % height_html
        height_js = "%s" % height_html

    # Whether to show certain UI elements or not
    if show_tooltips == False:
        tooltips_display = "display: none;"
    else:
        tooltips_display = ""

    if show_meta == False:
        meta_display = "display: none;"
    else:
        meta_display = ""

    if show_title == False:
        title_display = "display: none;"
    else:
        title_display = ""

    with open(path_html, "wb") as outfile:
        html = open('./template.html').read() % (
        title, width_css, height_css, title_display, meta_display, tooltips_display,
        width_js, height_js, graph_charge, graph_gravity, json.dumps(json_s))
        outfile.write(html.encode("utf-8"))
