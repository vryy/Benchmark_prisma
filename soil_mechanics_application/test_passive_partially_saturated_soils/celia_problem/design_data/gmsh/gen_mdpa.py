import gmshparser

from KratosMultiphysics import *
from KratosMultiphysics.LayerApplication import *

import gmsh_mdpa_interface

# name = "ground_h27"

# mesh = gmshparser.parse(name + ".msh")

# dim = 3
# mdpa_writer = gmsh_mdpa_interface.GmshMDPAWriter(dim, mesh)
# mdpa_writer.Print()

# mdpa_writer.SetElementTag(1, gmsh_mdpa_interface.VOLUME, "PassivePartiallySaturatedSoilsKinematicLinear")

# mdpa_writer.SetLayerTag(1, gmsh_mdpa_interface.VOLUME, "Ground")

# mdpa_writer.WriteMDPA(name)
# mdpa_writer.WriteLayers(name)

######

name = "ground_l2"

mesh = gmshparser.parse(name + ".msh")

dim = 1
mdpa_writer = gmsh_mdpa_interface.GmshMDPAWriter(dim, mesh)
mdpa_writer.Print()

mdpa_writer.SetElementTag(1, gmsh_mdpa_interface.EDGE, "RichardsElement")

mdpa_writer.SetLayerTag(1, gmsh_mdpa_interface.EDGE, "Ground")

mdpa_writer.WriteMDPA(name)
mdpa_writer.WriteLayers(name)
