def ReadLayerSets():
    ##################################################################
    ## STORE LAYER SETS ##############################################
    ##################################################################
    ## ELEMENTS on layers ############################################
    layer_sets = {}
    layer_elements_list = [
    ]
    layer_sets['Layer0'] = layer_elements_list
    layer_elements_list = [
    1 ,
    ]
    layer_sets['down'] = layer_elements_list
    layer_elements_list = [
    2 ,
    ]
    layer_sets['up'] = layer_elements_list
    layer_elements_list = [
    ]
    layer_sets['master'] = layer_elements_list
    layer_elements_list = [
    ]
    layer_sets['slave'] = layer_elements_list
    return layer_sets

def ReadInnerBoundaryElements():
    ## ELEMENTS on inner boundaries ##################################
    inner_boundary_elements = [
    ]
    return inner_boundary_elements

def ReadLayerNodesSets():
    ## NODES on layers ###############################################
    layer_nodes_sets = {}
    layer_nodes_list = [
    ]
    layer_nodes_sets['Layer0'] = layer_nodes_list
    layer_nodes_list = [
    13 ,
    14 ,
    15 ,
    16 ,
    ]
    layer_nodes_sets['down'] = layer_nodes_list
    layer_nodes_list = [
    1 ,
    4 ,
    5 ,
    10 ,
    ]
    layer_nodes_sets['up'] = layer_nodes_list
    layer_nodes_list = [
    2 ,
    6 ,
    7 ,
    11 ,
    ]
    layer_nodes_sets['master'] = layer_nodes_list
    layer_nodes_list = [
    3 ,
    8 ,
    9 ,
    12 ,
    ]
    layer_nodes_sets['slave'] = layer_nodes_list
    return layer_nodes_sets

def ReadContactMasterNodes():
    ## CONTACT MASTER NODES ##########################################
    contact_master_nodes = [
    ]
    return contact_master_nodes

def ReadContactSlaveNodes():
    ## CONTACT SLAVE NODES ###########################################
    contact_slave_nodes = [
    ]
    return contact_slave_nodes

def ReadInnerBoundaryNodes():
    ## INNER BOUNDARY NODES ##########################################
    inner_boundary_nodes = [
    ]
    return inner_boundary_nodes

def ReadTopSurfaceNodes():
    ##################################################################
    ## STORE NODES ON GROUND SURFACE #################################
    ##################################################################
    top_surface_nodes = []
    return top_surface_nodes

def ReadNodeGroups():
    ##################################################################
    ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
    ##################################################################
    node_groups = {}
    if( not (' master ' in node_groups) ):
        node_groups['master'] = []
    node_groups['master'].append( 2 )
    if( not ('slave' in node_groups) ):
        node_groups['slave'] = []
    node_groups['slave'].append( 3 )
    if( not ('master' in node_groups) ):
        node_groups['master'] = []
    node_groups['master'].append( 6 )
    if( not ('master' in node_groups) ):
        node_groups['master'] = []
    node_groups['master'].append( 7 )
    if( not ('slave' in node_groups) ):
        node_groups['slave'] = []
    node_groups['slave'].append( 8 )
    if( not ('slave' in node_groups) ):
        node_groups['slave'] = []
    node_groups['slave'].append( 9 )
    if( not ('master' in node_groups) ):
        node_groups['master'] = []
    node_groups['master'].append( 11 )
    if( not ('slave' in node_groups) ):
        node_groups['slave'] = []
    node_groups['slave'].append( 12 )
    return node_groups

