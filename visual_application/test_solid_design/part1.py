dd = designDoc
D1 = 30
D2 = 20
Hc = 50
L = 50
W = 50
H = 10
dd.addBox(0, 0, 0, 0, L, W, H)
dd.addCylinder(0, L/2, W/2, Hc/10, D1/2, Hc)
dd.subtractCylinder(0, L/2, W/2, -Hc/10, D2/2, Hc*1.2)

elem_size = 1
dd.setMeshSize(elem_size)
mesh_params = {'edge_size': elem_size, \
    'facet_angle': 25, \
    'facet_size': elem_size, \
    'facet_distance': elem_size/10, \
    'cell_radius_edge_ratio': 3, \
    'cell_size': elem_size, \
    'refine': 2}
dd.createMeshWithFeatures(mesh_params)

dd.setEdgeMode(1)
dd.zoomToFit()
