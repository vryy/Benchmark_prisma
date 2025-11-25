import math

dd = designDoc
D1 = 80
D2 = 100
D3 = 8
Hc = 10

dd.pauseRendering()

elem_size = 1
dd.setMeshSize(0, elem_size)

dd.addCylinder(0, 0, 0, 0, D2/2, Hc)
dd.subtractCylinder(0, 0, 0, -Hc/10, D1/2, Hc*1.2)

nb = 10
for i in range(0, nb):
    a = 2*math.pi/nb
    c = D2/2*math.cos(i*a)
    s = D2/2*math.sin(i*a)
    dd.addSphere(0, c, s, Hc/2, D3/2)

# mesh_params = {'edge_size': elem_size, \
#     'facet_angle': 25, \
#     'facet_size': elem_size, \
#     'facet_distance': elem_size/10, \
#     'cell_radius_edge_ratio': 3, \
#     'cell_size': elem_size, \
#     'refine': 2}
# dd.createMeshWithFeatures(mesh_params)

dd.buildDataSet()
dd.render()

dd.setEdgeMode(1)
dd.zoomToFit()
