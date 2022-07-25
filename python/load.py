import bpy

#vertices = [
#  [0,0,0],
#  [1,0.57735026,0],
#  [0.33333334,0.57735026,1],
#  [0.6666666,0,1],
#]
#faces = [
#  [0, 1, 2, ],
#  [1, 0, 3, ],
#]

exec(open("/tmp/geom.py").read(), globals(), locals())

name = "maze"
mesh = bpy.data.meshes.new(name)
mesh.from_pydata(vertices, [] , faces)

obj = bpy.data.objects.get(name)
if obj is None:
  obj = bpy.data.objects.new(name, mesh)
else:
  obj.data = mesh

try:
  bpy.context.scene.collection.objects.link(obj)
except:
  pass

