def ReadLayerSets():
    ##################################################################
    ## STORE LAYER SETS ##############################################
    ##################################################################
    ## ELEMENTS on layers ############################################
    layer_sets = {}
    layer_elements_list = [
    1 ,
    2 ,
    3 ,
    4 ,
    5 ,
    6 ,
    7 ,
    8 ,
    9 ,
    10 ,
    11 ,
    12 ,
    13 ,
    14 ,
    15 ,
    16 ,
    17 ,
    18 ,
    19 ,
    20 ,
    21 ,
    22 ,
    23 ,
    24 ,
    25 ,
    26 ,
    27 ,
    28 ,
    29 ,
    30 ,
    31 ,
    32 ,
    33 ,
    34 ,
    35 ,
    ]
    layer_sets['Layer0'] = layer_elements_list
    layer_elements_list = [
    ]
    layer_sets['base'] = layer_elements_list
    layer_elements_list = [
    ]
    layer_sets['load'] = layer_elements_list
    return layer_sets

def ReadLayerNodesSets():
    ## NODES on layers ###############################################
    layer_nodes_sets = {}
    layer_nodes_list = [
    4 ,
    6 ,
    7 ,
    8 ,
    11 ,
    12 ,
    15 ,
    16 ,
    17 ,
    18 ,
    19 ,
    20 ,
    21 ,
    22 ,
    23 ,
    24 ,
    25 ,
    26 ,
    27 ,
    28 ,
    30 ,
    31 ,
    32 ,
    33 ,
    34 ,
    35 ,
    36 ,
    37 ,
    38 ,
    39 ,
    40 ,
    41 ,
    42 ,
    43 ,
    44 ,
    45 ,
    46 ,
    47 ,
    48 ,
    49 ,
    50 ,
    51 ,
    52 ,
    53 ,
    54 ,
    55 ,
    56 ,
    57 ,
    58 ,
    59 ,
    60 ,
    61 ,
    62 ,
    63 ,
    64 ,
    65 ,
    66 ,
    67 ,
    68 ,
    69 ,
    70 ,
    71 ,
    72 ,
    76 ,
    77 ,
    83 ,
    ]
    layer_nodes_sets['Layer0'] = layer_nodes_list
    layer_nodes_list = [
    1 ,
    2 ,
    3 ,
    5 ,
    9 ,
    10 ,
    13 ,
    14 ,
    29 ,
    ]
    layer_nodes_sets['base'] = layer_nodes_list
    layer_nodes_list = [
    73 ,
    74 ,
    75 ,
    78 ,
    79 ,
    80 ,
    81 ,
    82 ,
    84 ,
    85 ,
    86 ,
    87 ,
    88 ,
    89 ,
    90 ,
    91 ,
    ]
    layer_nodes_sets['load'] = layer_nodes_list
    return layer_nodes_sets

def ReadBoundaryNodes():
    ## BOUNDARY NODES ##########################################
    boundary_nodes = {}
    boundary_nodes['IsBoundary'] = []
    boundary_nodes['IsBottom'] = []
    boundary_nodes['IsSide'] = []
    boundary_nodes['IsTop'] = []

    return boundary_nodes

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

def ReadLiningEndNodes():
    ## LINING END NODES ###########################################
    lining_end_nodes = [
    ]
    return lining_end_nodes

def ReadTBMEndNodes():
    ## TBM END NODES ###########################################
    tbm_end_nodes = [
    ]
    return tbm_end_nodes

def ReadInnerBoundaryNodes():
    ## INNER BOUNDARY NODES ##########################################
    inner_boundary_nodes = [
    ]
    return inner_boundary_nodes

def ReadTopSurfaceNodes():
    ##################################################################
    ## STORE NODES ON GROUND SURFACE #################################
    ##################################################################
    top_surface_nodes = [
    ]
    return top_surface_nodes

def ReadNodeGroups():
    ##################################################################
    ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
    ##################################################################
    node_groups = {}

    return node_groups

