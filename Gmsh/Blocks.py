import gmsh

gmsh.initialize()

gmsh.model.add('Stationary-Block')

meshSize: float = 20

# Add nodes
gmsh.model.geo.addPoint(    0, 0,  0, meshSize, 1)
gmsh.model.geo.addPoint(    0, 0, 50, meshSize, 2)
gmsh.model.geo.addPoint( -250, 0, 50, meshSize, 3)
gmsh.model.geo.addPoint( -250, 0,  0, meshSize, 4)


gmsh.model.geo.addPoint(    0, 1500,  0, meshSize, 5)
gmsh.model.geo.addPoint(    0, 1500, 50, meshSize, 6)
gmsh.model.geo.addPoint( -250, 1500, 50, meshSize, 7)
gmsh.model.geo.addPoint( -250, 1500,  0, meshSize, 8)


# Add lines
gmsh.model.geo.addLine(1, 2, 1) 
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)


gmsh.model.geo.addLine(5, 6, 5) 
gmsh.model.geo.addLine(6, 7, 6)
gmsh.model.geo.addLine(7, 8, 7)
gmsh.model.geo.addLine(8, 5, 8)

gmsh.model.geo.addLine(1, 5,  9) 
gmsh.model.geo.addLine(2, 6, 10)
gmsh.model.geo.addLine(3, 7, 11)
gmsh.model.geo.addLine(4, 8, 12)


# Add Curve Loop
gmsh.model.geo.addCurveLoop([ 1,  2,  3,   4], 1)
gmsh.model.geo.addCurveLoop([ 5,  6,  7,   8], 2)
gmsh.model.geo.addCurveLoop([ 2, 11, -6, -10], 3)
gmsh.model.geo.addCurveLoop([ 1, 10, -5,  -9], 4)
gmsh.model.geo.addCurveLoop([-4, 12,  8,  -9], 5)
gmsh.model.geo.addCurveLoop([-3, 11,  7, -12], 6)

# Add Surfaces
gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.geo.addPlaneSurface([2], 2)
gmsh.model.geo.addPlaneSurface([3], 3)
gmsh.model.geo.addPlaneSurface([4], 4)
gmsh.model.geo.addPlaneSurface([5], 5)
gmsh.model.geo.addPlaneSurface([6], 6)


# Add Volume
gmsh.model.geo.addSurfaceLoop([1, 2, 3, 4, 5, 6], 1)
gmsh.model.geo.addVolume([1])



gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("Stationary-Block.msh")

gmsh.fltk.run()
gmsh.finalize()