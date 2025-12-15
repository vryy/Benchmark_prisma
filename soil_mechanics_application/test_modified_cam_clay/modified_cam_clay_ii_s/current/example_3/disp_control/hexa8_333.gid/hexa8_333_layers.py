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
    node_groups['loady'] = []
    node_groups['loadz'] = []
    node_groups['loadx'] = []

    #  loady 
    node_groups['loady'].append( 1 )
    # loadz
    node_groups['loadz'].append( 1 )
    # loady
    node_groups['loady'].append( 2 )
    # loady
    node_groups['loady'].append( 3 )
    # loadz
    node_groups['loadz'].append( 3 )
    # loadz
    node_groups['loadz'].append( 4 )
    # loadz
    node_groups['loadz'].append( 5 )
    # loady
    node_groups['loady'].append( 6 )
    # loady
    node_groups['loady'].append( 9 )
    # loady
    node_groups['loady'].append( 10 )
    # loadz
    node_groups['loadz'].append( 10 )
    # loadz
    node_groups['loadz'].append( 11 )
    # loadz
    node_groups['loadz'].append( 12 )
    # loady
    node_groups['loady'].append( 14 )
    # loadz
    node_groups['loadz'].append( 15 )
    # loady
    node_groups['loady'].append( 17 )
    # loady
    node_groups['loady'].append( 21 )
    # loadz
    node_groups['loadz'].append( 22 )
    # loady
    node_groups['loady'].append( 26 )
    # loadz
    node_groups['loadz'].append( 27 )
    # loadx
    node_groups['loadx'].append( 28 )
    # loady
    node_groups['loady'].append( 28 )
    # loadz
    node_groups['loadz'].append( 28 )
    # loady
    node_groups['loady'].append( 30 )
    # loadz
    node_groups['loadz'].append( 33 )
    # loadx
    node_groups['loadx'].append( 34 )
    # loady
    node_groups['loady'].append( 34 )
    # loadx
    node_groups['loadx'].append( 35 )
    # loadz
    node_groups['loadz'].append( 35 )
    # loadx
    node_groups['loadx'].append( 38 )
    # loadx
    node_groups['loadx'].append( 41 )
    # loadz
    node_groups['loadz'].append( 41 )
    # loady
    node_groups['loady'].append( 42 )
    # loadz
    node_groups['loadz'].append( 44 )
    # loadx
    node_groups['loadx'].append( 45 )
    # loady
    node_groups['loady'].append( 45 )
    # loadx
    node_groups['loadx'].append( 47 )
    # loadx
    node_groups['loadx'].append( 48 )
    # loadx
    node_groups['loadx'].append( 54 )
    # loadx
    node_groups['loadx'].append( 56 )
    # loady
    node_groups['loady'].append( 56 )
    # loadx
    node_groups['loadx'].append( 57 )
    # loadz
    node_groups['loadz'].append( 57 )
    # loadx
    node_groups['loadx'].append( 59 )
    # loadx
    node_groups['loadx'].append( 60 )
    # loadx
    node_groups['loadx'].append( 62 )
    # loadx
    node_groups['loadx'].append( 63 )
    # loadx
    node_groups['loadx'].append( 64 )
    return node_groups

