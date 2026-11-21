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
    ]
    layer_sets['Layer0'] = layer_elements_list
    layer_elements_list = [
    ]
    layer_sets['inner'] = layer_elements_list
    layer_elements_list = [
    ]
    layer_sets['outer'] = layer_elements_list
    return layer_sets

def ReadLayerNodesSets():
    ## NODES on layers ###############################################
    layer_nodes_sets = {}
    layer_nodes_list = [
    4 ,
    5 ,
    8 ,
    9 ,
    10 ,
    11 ,
    12 ,
    13 ,
    14 ,
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
    49 ,
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
    77 ,
    79 ,
    80 ,
    89 ,
    90 ,
    91 ,
    92 ,
    95 ,
    96 ,
    97 ,
    98 ,
    99 ,
    100 ,
    101 ,
    102 ,
    105 ,
    106 ,
    108 ,
    109 ,
    113 ,
    114 ,
    115 ,
    118 ,
    119 ,
    120 ,
    121 ,
    124 ,
    126 ,
    127 ,
    128 ,
    133 ,
    134 ,
    135 ,
    136 ,
    137 ,
    138 ,
    139 ,
    140 ,
    141 ,
    142 ,
    143 ,
    144 ,
    145 ,
    146 ,
    147 ,
    151 ,
    152 ,
    156 ,
    157 ,
    158 ,
    159 ,
    162 ,
    163 ,
    167 ,
    168 ,
    169 ,
    170 ,
    171 ,
    172 ,
    173 ,
    174 ,
    175 ,
    177 ,
    179 ,
    180 ,
    183 ,
    184 ,
    185 ,
    186 ,
    187 ,
    188 ,
    189 ,
    190 ,
    191 ,
    198 ,
    199 ,
    200 ,
    201 ,
    202 ,
    203 ,
    204 ,
    205 ,
    206 ,
    209 ,
    210 ,
    211 ,
    212 ,
    213 ,
    214 ,
    215 ,
    221 ,
    222 ,
    223 ,
    224 ,
    225 ,
    226 ,
    227 ,
    230 ,
    231 ,
    232 ,
    233 ,
    234 ,
    235 ,
    236 ,
    240 ,
    242 ,
    245 ,
    246 ,
    247 ,
    248 ,
    249 ,
    250 ,
    251 ,
    252 ,
    253 ,
    254 ,
    255 ,
    256 ,
    257 ,
    258 ,
    259 ,
    260 ,
    261 ,
    264 ,
    265 ,
    266 ,
    267 ,
    268 ,
    269 ,
    270 ,
    271 ,
    272 ,
    276 ,
    277 ,
    278 ,
    279 ,
    280 ,
    283 ,
    284 ,
    ]
    layer_nodes_sets['Layer0'] = layer_nodes_list
    layer_nodes_list = [
    81 ,
    82 ,
    83 ,
    84 ,
    85 ,
    86 ,
    87 ,
    88 ,
    93 ,
    94 ,
    103 ,
    104 ,
    107 ,
    116 ,
    117 ,
    122 ,
    123 ,
    125 ,
    131 ,
    132 ,
    148 ,
    149 ,
    150 ,
    160 ,
    161 ,
    164 ,
    165 ,
    166 ,
    181 ,
    182 ,
    192 ,
    193 ,
    194 ,
    207 ,
    208 ,
    216 ,
    217 ,
    219 ,
    228 ,
    229 ,
    241 ,
    243 ,
    244 ,
    ]
    layer_nodes_sets['inner'] = layer_nodes_list
    layer_nodes_list = [
    1 ,
    2 ,
    3 ,
    6 ,
    7 ,
    15 ,
    16 ,
    17 ,
    30 ,
    31 ,
    48 ,
    50 ,
    51 ,
    76 ,
    78 ,
    110 ,
    111 ,
    112 ,
    129 ,
    130 ,
    153 ,
    154 ,
    155 ,
    176 ,
    178 ,
    195 ,
    196 ,
    197 ,
    218 ,
    220 ,
    237 ,
    238 ,
    239 ,
    262 ,
    263 ,
    273 ,
    274 ,
    275 ,
    281 ,
    282 ,
    285 ,
    286 ,
    287 ,
    ]
    layer_nodes_sets['outer'] = layer_nodes_list
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
    node_groups['outer'] = []
    node_groups['inner'] = []

    #  outer 
    node_groups['outer'].append( 1 )
    # outer
    node_groups['outer'].append( 2 )
    # outer
    node_groups['outer'].append( 3 )
    # outer
    node_groups['outer'].append( 6 )
    # outer
    node_groups['outer'].append( 7 )
    # outer
    node_groups['outer'].append( 15 )
    # outer
    node_groups['outer'].append( 16 )
    # outer
    node_groups['outer'].append( 17 )
    # outer
    node_groups['outer'].append( 30 )
    # outer
    node_groups['outer'].append( 31 )
    # outer
    node_groups['outer'].append( 48 )
    # outer
    node_groups['outer'].append( 50 )
    # outer
    node_groups['outer'].append( 51 )
    # outer
    node_groups['outer'].append( 76 )
    # outer
    node_groups['outer'].append( 78 )
    # inner
    node_groups['inner'].append( 81 )
    # inner
    node_groups['inner'].append( 82 )
    # inner
    node_groups['inner'].append( 83 )
    # inner
    node_groups['inner'].append( 84 )
    # inner
    node_groups['inner'].append( 85 )
    # inner
    node_groups['inner'].append( 86 )
    # inner
    node_groups['inner'].append( 87 )
    # inner
    node_groups['inner'].append( 88 )
    # inner
    node_groups['inner'].append( 93 )
    # inner
    node_groups['inner'].append( 94 )
    # inner
    node_groups['inner'].append( 103 )
    # inner
    node_groups['inner'].append( 104 )
    # inner
    node_groups['inner'].append( 107 )
    # outer
    node_groups['outer'].append( 110 )
    # outer
    node_groups['outer'].append( 111 )
    # outer
    node_groups['outer'].append( 112 )
    # inner
    node_groups['inner'].append( 116 )
    # inner
    node_groups['inner'].append( 117 )
    # inner
    node_groups['inner'].append( 122 )
    # inner
    node_groups['inner'].append( 123 )
    # inner
    node_groups['inner'].append( 125 )
    # outer
    node_groups['outer'].append( 129 )
    # outer
    node_groups['outer'].append( 130 )
    # inner
    node_groups['inner'].append( 131 )
    # inner
    node_groups['inner'].append( 132 )
    # inner
    node_groups['inner'].append( 148 )
    # inner
    node_groups['inner'].append( 149 )
    # inner
    node_groups['inner'].append( 150 )
    # outer
    node_groups['outer'].append( 153 )
    # outer
    node_groups['outer'].append( 154 )
    # outer
    node_groups['outer'].append( 155 )
    # inner
    node_groups['inner'].append( 160 )
    # inner
    node_groups['inner'].append( 161 )
    # inner
    node_groups['inner'].append( 164 )
    # inner
    node_groups['inner'].append( 165 )
    # inner
    node_groups['inner'].append( 166 )
    # outer
    node_groups['outer'].append( 176 )
    # outer
    node_groups['outer'].append( 178 )
    # inner
    node_groups['inner'].append( 181 )
    # inner
    node_groups['inner'].append( 182 )
    # inner
    node_groups['inner'].append( 192 )
    # inner
    node_groups['inner'].append( 193 )
    # inner
    node_groups['inner'].append( 194 )
    # outer
    node_groups['outer'].append( 195 )
    # outer
    node_groups['outer'].append( 196 )
    # outer
    node_groups['outer'].append( 197 )
    # inner
    node_groups['inner'].append( 207 )
    # inner
    node_groups['inner'].append( 208 )
    # inner
    node_groups['inner'].append( 216 )
    # inner
    node_groups['inner'].append( 217 )
    # outer
    node_groups['outer'].append( 218 )
    # inner
    node_groups['inner'].append( 219 )
    # outer
    node_groups['outer'].append( 220 )
    # inner
    node_groups['inner'].append( 228 )
    # inner
    node_groups['inner'].append( 229 )
    # outer
    node_groups['outer'].append( 237 )
    # outer
    node_groups['outer'].append( 238 )
    # outer
    node_groups['outer'].append( 239 )
    # inner
    node_groups['inner'].append( 241 )
    # inner
    node_groups['inner'].append( 243 )
    # inner
    node_groups['inner'].append( 244 )
    # outer
    node_groups['outer'].append( 262 )
    # outer
    node_groups['outer'].append( 263 )
    # outer
    node_groups['outer'].append( 273 )
    # outer
    node_groups['outer'].append( 274 )
    # outer
    node_groups['outer'].append( 275 )
    # outer
    node_groups['outer'].append( 281 )
    # outer
    node_groups['outer'].append( 282 )
    # outer
    node_groups['outer'].append( 285 )
    # outer
    node_groups['outer'].append( 286 )
    # outer
    node_groups['outer'].append( 287 )
    return node_groups

