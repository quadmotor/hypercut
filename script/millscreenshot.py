import polyscope as ps
import meshio


dir = "C:/tmp/hc/"

# Adjust some screenshot default settings if you'd like
ps.set_screenshot_extension(".jpg");
ps.init()
# ps.set_screenshot_extension(".jpg");

jsonstr = '{"farClipRatio":20.0,"fov":45.0,"nearClipRatio":0.005,"projectionMode":"Perspective","viewMat":[1.0,0.0,0.0,-500.0,0.0,0.997785151004791,-0.0665190070867538,-465.633117675781,0.0,0.0665190070867538,0.997785151004791,-2256.9814453125,0.0,0.0,0.0,1.0],"windowHeight":818,"windowWidth":1184}'
ps.set_view_from_json(jsonstr)

for i in range(100, 17900, 100):
    fn16500 = dir + "solid{}.obj".format(i)

    mesh = meshio.read(fn16500)

    crd3 = mesh.points
    print(crd3, mesh.cells_dict)
    tri3 = mesh.cells_dict["quad"]

    ps_mesh = ps.register_surface_mesh("my mesh", crd3, tri3)
    ps_mesh.set_edge_color((0, 0, 0)) 
    ps_mesh.set_edge_width(1.0)
    # ps.show()



    # Take a screenshot
    # It will be written to your current directory as screenshot_000000.jpg, etc
    ps.screenshot()

    ps_mesh.remove()
    
