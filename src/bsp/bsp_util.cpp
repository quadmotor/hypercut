#include "bsp_util.h"





void random_sphere(std::string outputfn, const int nofplanes)
{
    std::vector<int> randomx = { 21,1,13,18,4,2,16,6,9,20,7,26,26,1,16,17,15,24,27,7,3,13,3,12,24,19,23,18,5,9,12,1,7,6,18,2,14,0,28,2,15,16,10,13,5,17,22,9,18,6,17,0,3,9,22,12,6,0,14,6,9,1,12,7,23,10,20,20,23,19,12,10,0,15,2,25,12,2,18,28,26,7,16,29,0,13,1,14,15,12,15,3,6,28,14,29,4,11,2, 0 };


    randomx.clear();
    for (int i = 0; i < nofplanes; i++)
    {
        randomx.push_back(rand() / 30);
    }

    Mesh sh;


    std::vector<pmp::Vertex> vertices;



    std::vector< BspPlane> sphereplanes;

    for (int i = 0; i < randomx.size() / 3; i++)
    {
        int* p = &(randomx[i * 3]);
        Eigen::Vector3d pt;
        pt << p[0], p[1], p[2];

        pt.normalize();
        pt *= 900;

        BspPlane bsplane = { int(pt[0]), int(pt[1]), int(pt[2]), -900 * 20 };

        addrawquad(sh, bsplane, vertices);

        sphereplanes.push_back(bsplane);
    }

    if (0)
        sh.write("randomsphere.obj");


    {
        int failcnt = 0;
        Bsp bsp;
        Mesh sh;
        build_cube(bsp, sh, 40);

        for (int i = 0; i < sphereplanes.size(); i++)
        {

            auto pl = sphereplanes[i];
            if (0 != convexmesh_cut(bsp, sh, pl))
                failcnt++;

            if (0)
                dumpmesh(i, sh);

            if (0 == i % 100) {
                // printf("%d %d\n", i, failcnt);
            }
        }

        // printf("failed count %d %d", failcnt, int(sphereplanes.size()));

        if (!outputfn.empty()) {
            dumpmesh(sh, outputfn, 0, true);
        }
    }
}