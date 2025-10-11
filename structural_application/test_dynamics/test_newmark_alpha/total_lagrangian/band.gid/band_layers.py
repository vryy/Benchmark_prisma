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
    ]
    layer_sets['Layer0'] = layer_elements_list
    return layer_sets

def ReadLayerNodesSets():
    ## NODES on layers ###############################################
    layer_nodes_sets = {}
    layer_nodes_list = [
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
    73 ,
    74 ,
    75 ,
    76 ,
    ]
    layer_nodes_sets['Layer0'] = layer_nodes_list
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
    node_groups['Interface'] = []

    #  Interface 
    node_groups['Interface'].append( 1 )
    # Interface
    node_groups['Interface'].append( 1 )
    # Interface
    node_groups['Interface'].append( 2 )
    # Interface
    node_groups['Interface'].append( 3 )
    # Interface
    node_groups['Interface'].append( 4 )
    # Interface
    node_groups['Interface'].append( 5 )
    # Interface
    node_groups['Interface'].append( 5 )
    # Interface
    node_groups['Interface'].append( 8 )
    # Interface
    node_groups['Interface'].append( 9 )
    # Interface
    node_groups['Interface'].append( 12 )
    # Interface
    node_groups['Interface'].append( 13 )
    # Interface
    node_groups['Interface'].append( 16 )
    # Interface
    node_groups['Interface'].append( 17 )
    # Interface
    node_groups['Interface'].append( 20 )
    # Interface
    node_groups['Interface'].append( 21 )
    # Interface
    node_groups['Interface'].append( 24 )
    # Interface
    node_groups['Interface'].append( 25 )
    # Interface
    node_groups['Interface'].append( 28 )
    # Interface
    node_groups['Interface'].append( 29 )
    # Interface
    node_groups['Interface'].append( 32 )
    # Interface
    node_groups['Interface'].append( 33 )
    # Interface
    node_groups['Interface'].append( 36 )
    # Interface
    node_groups['Interface'].append( 37 )
    # Interface
    node_groups['Interface'].append( 40 )
    # Interface
    node_groups['Interface'].append( 41 )
    # Interface
    node_groups['Interface'].append( 44 )
    # Interface
    node_groups['Interface'].append( 45 )
    # Interface
    node_groups['Interface'].append( 48 )
    # Interface
    node_groups['Interface'].append( 49 )
    # Interface
    node_groups['Interface'].append( 52 )
    # Interface
    node_groups['Interface'].append( 53 )
    # Interface
    node_groups['Interface'].append( 56 )
    # Interface
    node_groups['Interface'].append( 57 )
    # Interface
    node_groups['Interface'].append( 60 )
    # Interface
    node_groups['Interface'].append( 61 )
    # Interface
    node_groups['Interface'].append( 64 )
    # Interface
    node_groups['Interface'].append( 65 )
    # Interface
    node_groups['Interface'].append( 68 )
    # Interface
    node_groups['Interface'].append( 69 )
    # Interface
    node_groups['Interface'].append( 72 )
    # Interface
    node_groups['Interface'].append( 73 )
    # Interface
    node_groups['Interface'].append( 76 )
    return node_groups

