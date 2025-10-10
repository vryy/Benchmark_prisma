def ReadLayerSets():
    ##################################################################
    ## STORE LAYER SETS ##############################################
    ##################################################################
    ## ELEMENTS on layers ############################################
    layer_sets = {}
    layer_elements_list = [
    1 ,
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

    return node_groups

